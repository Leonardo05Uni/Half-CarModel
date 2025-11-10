# -*- coding: utf-8 -*-

"""
Half-car (4-DOF) suspension model driving over a measured road profile.

DOFs:
    - Body heave:      z       [m]   (positive upwards)
    - Body pitch:      theta   [rad] (small-angle approximation)
    - Front wheel:     z_wf    [m]
    - Rear wheel:      z_wr    [m]

Includes:
    - Suspension springs and dampers (front & rear)
    - Tyre vertical stiffness (front & rear)
    - Road profile from CSV, smoothed with a spline
    - Vehicle travelling at constant speed along the road

Main outputs:
    - Time histories of z, theta, wheel motions
    - Vertical acceleration at passenger location
    - RMS acceleration metrics for comfort analysis
"""

from dataclasses import dataclass
from typing import Callable, Tuple

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import UnivariateSpline


# ---------------------------------------------------------------------------
# Data classes for parameters and simulation options
# ---------------------------------------------------------------------------

@dataclass
class CarParams:
    """Physical parameters of the 4-DOF half-car model."""

    # Sprung body properties
    body_M: float          # Sprung mass [kg]
    body_inertia: float    # Pitch inertia about CG [kg·m²]
    body_a: float          # Distance CG → front axle [m]  (positive forward)
    body_b: float          # Distance CG → rear axle [m]   (positive backward)

    # Suspension (between body and wheel) – front & rear
    FWS_k: float           # Front suspension spring stiffness [N/m]
    FWD_c: float           # Front suspension damping coefficient [N·s/m]
    RWS_k: float           # Rear suspension spring stiffness [N/m]
    RWD_c: float           # Rear suspension damping coefficient [N·s/m]

    # Geometric pickup locations (not used in current equations but kept for extension)
    FWP_theta: float       # Front pickup pitch coordinate (if needed later)
    FWP_z: float           # Front pickup vertical offset (if needed later)
    RWP_theta: float       # Rear pickup pitch coordinate (if needed later)
    RWP_z: float           # Rear pickup vertical offset (if needed later)

    # Unsprung masses (wheels)
    m_wf: float            # Front unsprung mass [kg]
    m_wr: float            # Rear unsprung mass [kg]

    # Tyre stiffness (vertical)
    k_tf: float            # Front tyre stiffness [N/m]
    k_tr: float            # Rear tyre stiffness [N/m]


@dataclass
class SimulationOptions:
    """Integrator options and initial conditions."""
    t_span: tuple          # (t0, t_final) [s]
    y_0: list              # Initial state vector (length 8)
    r_tolerance: float = 1e-7
    a_tolerance: float = 1e-9
    dense: bool = True     # Whether to store a dense interpolant sol.sol(t)


# ---------------------------------------------------------------------------
# Road profile: load CSV, build spline, then define a BaseInput function
# ---------------------------------------------------------------------------

# Type alias for a road base function:
#   base(t) -> (y_f, y_r, y_f_dot, y_r_dot)
BaseInput = Callable[[float], Tuple[float, float, float, float]]

# Load road profile (distance vs height) from CSV
df = pd.read_csv("bumpy_road_cords.csv")
x = df["distance"].values  # [m] along road
y = df["height"].values    # [m] road elevation

# Smooth interpolation of road profile
# s is a smoothing factor: larger s = smoother curve (less true to data)
spline = UnivariateSpline(x, y, s=0.4)
dsdx = spline.derivative()   # derivative dy/dx for computing dy/dt


def make_road_base(p: CarParams,
                   spline: UnivariateSpline,
                   dsdx: UnivariateSpline,
                   v: float = 8.0,
                   x0: float = 0.0,
                   clamp: bool = True) -> BaseInput:
    """
    Build a road input function for a vehicle travelling at constant speed v.

    Parameters
    ----------
    p : CarParams
        Car geometry (body_a, body_b) used to locate axles relative to CG.
    spline : UnivariateSpline
        Road height y(x).
    dsdx : UnivariateSpline
        Spatial derivative dy/dx of the road profile.
    v : float
        Vehicle speed [m/s].
    x0 : float
        CG longitudinal position at t = 0 [m].
    clamp : bool
        If True, clamp query points to profile domain [xmin, xmax].

    Returns
    -------
    base : function
        base(t) → (y_f, y_r, y_f_dot, y_r_dot) [m, m, m/s, m/s]
        Road height and vertical velocity under front and rear wheels.
    """
    xmin, xmax = float(np.min(x)), float(np.max(x))

    def base(t: float) -> Tuple[float, float, float, float]:
        # Longitudinal position of CG
        x_cg = x0 + v * t

        # Front & rear contact patch positions (front ahead by a, rear behind by b)
        x_f = x_cg + p.body_a
        x_r = x_cg - p.body_b

        if clamp:
            # Stop queries going outside the measured road profile
            x_fq = np.clip(x_f, xmin, xmax)
            x_rq = np.clip(x_r, xmin, xmax)
        else:
            x_fq, x_rq = x_f, x_r

        # Road height under front and rear wheels
        y_f = float(spline(x_fq))
        y_r = float(spline(x_rq))

        # Vertical road velocity: dy/dt = (dy/dx) * dx/dt = (dy/dx) * v
        y_f_dot = float(dsdx(x_fq)) * v
        y_r_dot = float(dsdx(x_rq)) * v

        return y_f, y_r, y_f_dot, y_r_dot

    return base


# ---------------------------------------------------------------------------
# Linearised body-only mass, damping, stiffness matrices (2-DOF body model)
# ---------------------------------------------------------------------------

def build_matrices_mck(p: CarParams):
    """
    Assemble M, C, K matrices for the 2-DOF sprung body (heave & pitch only).

    DOFs in order: [z, theta]

    This ignores unsprung masses and tyres and is mainly used for
    sanity-checking undamped natural frequencies of the body.
    """
    m, I = p.body_M, p.body_inertia
    a, b = p.body_a, p.body_b
    k_f, k_r = p.FWS_k, p.RWS_k
    c_f, c_r = p.FWD_c, p.RWD_c

    M = np.array([[m, 0.0],
                  [0.0, I]], dtype=float)

    C = np.array([[c_f + c_r,     a * c_f - b * c_r],
                  [a * c_f - b * c_r,  a * a * c_f + b * b * c_r]], dtype=float)

    K = np.array([[k_f + k_r,     a * k_f - b * k_r],
                  [a * k_f - b * k_r,  a * a * k_f + b * b * k_r]], dtype=float)

    return M, C, K


def undamped_naturals(p: CarParams):
    """
    Undamped natural frequencies (Hz) of the 2-DOF body-only model.

    Solves K φ = λ M φ, λ = ω².

    Note: This is an approximate check for the sprung body only and does not
    include the unsprung masses or tyre stiffnesses.
    """
    M, _, K = build_matrices_mck(p)
    lam, _ = np.linalg.eig(np.linalg.solve(M, K))
    # Clip possible small negative numerical eigenvalues to zero before sqrt
    wn = np.sqrt(np.clip(lam, 0.0, None))     # rad/s
    return np.sort(wn / (2 * np.pi))          # Hz


# ---------------------------------------------------------------------------
# Equations of motion for the 4-DOF half-car model
# ---------------------------------------------------------------------------

def rhs_car(t: float, x: np.ndarray, p: CarParams, base: BaseInput) -> np.ndarray:
    """
    Right-hand side of the ODE: x_dot = f(t, x).

    State vector:
        x = [ z, theta, z_wf, z_wr,  z_dot, theta_dot, z_wf_dot, z_wr_dot ]

    where:
        z           : body heave [m]
        theta       : body pitch [rad]
        z_wf, z_wr  : front and rear wheel vertical positions [m]
        z_dot       : body heave velocity [m/s]
        theta_dot   : body pitch rate [rad/s]
        z_wf_dot    : front wheel vertical velocity [m/s]
        z_wr_dot    : rear wheel vertical velocity [m/s]

    The forces are derived from:
        - Suspension springs/dampers between body and wheels
        - Tyre springs between wheels and road
        - Small-angle approximation for mapping pitch to pickup vertical motion
    """
    # Unpack state
    z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot = x

    # Road height and vertical velocity under front & rear tyres
    y_f, y_r, y_f_dot, y_r_dot = base(t)

    # Shorthand for parameters
    a, b = p.body_a, p.body_b
    k_f, k_r = p.FWS_k, p.RWS_k
    c_f, c_r = p.FWD_c, p.RWD_c
    M, I = p.body_M, p.body_inertia
    m_wf, m_wr = p.m_wf, p.m_wr
    k_tf, k_tr = p.k_tf, p.k_tr

    # ------------------------------------------------------------------
    # Suspension deflections & velocities (body relative to wheel)
    # Small-angle: front pickup at z + a*theta, rear pickup at z - b*theta
    # ------------------------------------------------------------------
    dL_f = (z + a * theta) - z_wf
    dL_f_dot = (z_dot + a * theta_dot) - z_wf_dot

    dL_r = (z - b * theta) - z_wr
    dL_r_dot = (z_dot - b * theta_dot) - z_wr_dot

    # Suspension forces on body (positive if suspension is extended upwards)
    F_s_f = k_f * dL_f + c_f * dL_f_dot
    F_s_r = k_r * dL_r + c_r * dL_r_dot

    # Tyre forces on wheel (road vs wheel vertical displacement)
    # Positive when tyre is compressed (road higher than wheel)
    F_t_f = k_tf * (y_f - z_wf)
    F_t_r = k_tr * (y_r - z_wr)

    # ------------------------------------------------------------------
    # Equations of motion
    # ------------------------------------------------------------------
    # Body heave:   M * z̈ = -(F_s_f + F_s_r)
    # Negative because a positive suspension extension (body higher than wheel)
    # produces a restoring force pulling the body back down.
    z_dot_dot = -(F_s_f + F_s_r) / M

    # Body pitch about CG:
    #   I * thetä = -(a * F_s_f - b * F_s_r)
    # Sign convention:
    #   - a*F_s_f tends to pitch the nose down if front suspension pulls down,
    #   - b*F_s_r tends to pitch the nose up or down depending on sign of F_s_r.
    theta_dot_dot = -(a * F_s_f - b * F_s_r) / I

    # Front unsprung mass (wheel):
    #   m_wf * z̈_wf = F_s_f + F_t_f
    z_wf_dot_dot = (F_s_f + F_t_f) / m_wf

    # Rear unsprung mass (wheel):
    #   m_wr * z̈_wr = F_s_r + F_t_r
    z_wr_dot_dot = (F_s_r + F_t_r) / m_wr

    # Pack derivatives
    return np.array([
        z_dot,            # d/dt z
        theta_dot,        # d/dt theta
        z_wf_dot,         # d/dt z_wf
        z_wr_dot,         # d/dt z_wr
        z_dot_dot,        # d/dt z_dot
        theta_dot_dot,    # d/dt theta_dot
        z_wf_dot_dot,     # d/dt z_wf_dot
        z_wr_dot_dot      # d/dt z_wr_dot
    ], dtype=float)


# ---------------------------------------------------------------------------
# Simulation helpers
# ---------------------------------------------------------------------------

def run_simulation(p: CarParams, base: BaseInput, opts: SimulationOptions):
    """
    Integrate the 4-DOF ODE system using solve_ivp.

    Returns
    -------
    sol : OdeSolution
        Object with sol.t and sol.y, and dense interpolant sol.sol(t) if enabled.
    """
    def fun(t, x):
        return rhs_car(t, x, p, base)

    sol = solve_ivp(fun,
                    opts.t_span,
                    opts.y_0,
                    rtol=opts.r_tolerance,
                    atol=opts.a_tolerance,
                    dense_output=opts.dense)
    return sol


def sample_states(sol, t_span, n=2000):
    """
    Sample the continuous solution on a uniform time grid.

    Parameters
    ----------
    sol : OdeSolution
        Output of run_simulation (solve_ivp with dense_output=True).
    t_span : tuple
        (t0, tf) [s]
    n : int
        Number of time points.

    Returns
    -------
    ts : ndarray
        Time vector [s]
    z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot : ndarrays
        State histories sampled at ts.
    """
    ts = np.linspace(t_span[0], t_span[1], n)
    z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot = sol.sol(ts)
    return ts, z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot


def accelerations_from_rhs(t_arr,
                           z_arr, theta_arr, z_wf_arr, z_wr_arr,
                           z_dot_arr, theta_dot_arr, z_wf_dot_arr, z_wr_dot_arr,
                           p: CarParams,
                           base: BaseInput):
    """
    Compute accelerations [z̈, θ̈, z̈_wf, z̈_wr] from the RHS for a given time history.

    This re-evaluates rhs_car at each time step using the stored state.
    It is numerically consistent with the ODE definition.
    """
    z_dot_dot = []
    theta_dot_dot = []
    zwf_dot_dot = []
    zwr_dot_dot = []

    for ti, zi, thi, zwfi, zwri, zdi, thdi, zwfdi, zwrdi in zip(
            t_arr,
            z_arr, theta_arr, z_wf_arr, z_wr_arr,
            z_dot_arr, theta_dot_arr, z_wf_dot_arr, z_wr_dot_arr):

        x_state = np.array([zi, thi, zwfi, zwri,
                            zdi, thdi, zwfdi, zwrdi], dtype=float)

        xdot = rhs_car(ti, x_state, p, base)

        # xdot indices:
        #   4 → z̈
        #   5 → θ̈
        #   6 → z̈_wf
        #   7 → z̈_wr
        z_dot_dot.append(xdot[4])
        theta_dot_dot.append(xdot[5])
        zwf_dot_dot.append(xdot[6])
        zwr_dot_dot.append(xdot[7])

    return (np.array(z_dot_dot),
            np.array(theta_dot_dot),
            np.array(zwf_dot_dot),
            np.array(zwr_dot_dot))


def occupant_vertical_accel(z_dot_dot, theta_dot_dot, x_from_CG: float = 0.0):
    """
    Vertical acceleration at a point on the body located x_from_CG [m] ahead of CG.

    Small-angle approximation:
        a_vert(x) = z̈ + x * θ̈
    """
    return z_dot_dot + x_from_CG * theta_dot_dot


def rms(a):
    """Root-mean-square of a 1D signal."""
    return float(np.sqrt(np.mean(a ** 2)))


# ---------------------------------------------------------------------------
# Example: define parameters, run simulation, and plot
# ---------------------------------------------------------------------------

# Nominal car parameters (order-of-magnitude realistic values)
p = CarParams(
    body_M=1200.0,          # Sprung mass [kg]
    body_inertia=2200.0,    # Pitch inertia [kg·m²]
    body_a=1.2,             # CG to front axle [m]
    body_b=1.3,             # CG to rear axle [m]

    FWS_k=35e3,             # Front suspension spring [N/m]
    FWD_c=500.0,            # Front damper [N·s/m]
    RWS_k=30e3,             # Rear suspension spring [N/m]
    RWD_c=500.0,            # Rear damper [N·s/m]

    m_wf=40.0,              # Front unsprung mass [kg]
    m_wr=35.0,              # Rear unsprung mass [kg]
    k_tf=2.0e5,             # Front tyre stiffness [N/m]
    k_tr=1.8e5,             # Rear tyre stiffness [N/m]

    FWP_theta=0.0,
    FWP_z=0.0,
    RWP_theta=0.0,
    RWP_z=0.0
)

# Build road input: car travelling at 8 m/s starting at x0 = 0
road_base = make_road_base(p, spline, dsdx, v=8.0, x0=0.0)

# Initial conditions:
# wheels start on the road, body initially level with zero velocity
y_f0, y_r0, _, _ = road_base(0.0)
opts = SimulationOptions(
    t_span=(0.0, 20.0),
    y_0=[
        0.0,     # z          (body heave)
        0.0,     # theta       (body pitch)
        y_f0,    # z_wf        (front wheel on road at t=0)
        y_r0,    # z_wr        (rear wheel on road at t=0)
        0.0,     # z_dot
        0.0,     # theta_dot
        0.0,     # z_wf_dot
        0.0      # z_wr_dot
    ]
)

# Run ODE simulation
sol = run_simulation(p, road_base, opts)

# Sample states uniformly in time
ts, z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot = sample_states(
    sol, opts.t_span, n=2000
)

# Compute accelerations using the same RHS as the ODE
z_dot_dot, theta_dot_dot, z_wf_dot_dot, z_wr_dot_dot = accelerations_from_rhs(
    ts,
    z, theta, z_wf, z_wr,
    z_dot, theta_dot, z_wf_dot, z_wr_dot,
    p,
    road_base
)

# Passenger location (ahead of CG) for comfort metric
seat_x = 0.8  # [m] forward of CG
a_pass = occupant_vertical_accel(z_dot_dot, theta_dot_dot, x_from_CG=seat_x)

# RMS comfort metrics
g = 9.81
rms_z_dot_dot = rms(z_dot_dot)
rms_pass = rms(a_pass)

print(f"Undamped body-only natural frequencies (approx) [Hz]: "
      f"{undamped_naturals(p)}")
print(f"RMS heave accel at CG: {rms_z_dot_dot:.4f} m/s^2  "
      f"({rms_z_dot_dot / g:.4f} g)")
print(f"RMS vertical accel at passenger (x={seat_x} m): "
      f"{rms_pass:.4f} m/s^2  ({rms_pass / g:.4f} g)")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

# (A) Body heave displacement
plt.figure()
plt.plot(ts, z)
plt.xlabel("Time [s]")
plt.ylabel("Heave z [m]")
plt.title("Body Heave vs Time")
plt.grid(True)

# (B) Body pitch angle in degrees
plt.figure()
plt.plot(ts, np.degrees(theta))
plt.xlabel("Time [s]")
plt.ylabel("Pitch θ [deg]")
plt.title("Body Pitch vs Time")
plt.grid(True)

# (C) Passenger vertical acceleration
plt.figure()
plt.plot(ts, a_pass)
plt.xlabel("Time [s]")
plt.ylabel("Vertical acceleration [m/s²]")
plt.title(f"Passenger Vertical Accel (x = {seat_x} m from CG)")
plt.grid(True)

# (D) Unsprung (wheel) accelerations
plt.figure()
plt.plot(ts, z_wf_dot_dot, label="Front Wheel")
plt.plot(ts, z_wr_dot_dot, label="Rear Wheel")
plt.xlabel("Time [s]")
plt.ylabel("Wheel acceleration [m/s²]")
plt.title("Unsprung (Wheel) Accelerations")
plt.grid(True)
plt.legend()

plt.show()
