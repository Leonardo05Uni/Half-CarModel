# -*- coding: utf-8 -*-

from dataclasses import dataclass
from typing import Callable, Tuple
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import UnivariateSpline

@dataclass
class CarParams:
    """All physical parameters of the 4-DOF car body model."""
    body_M: float
    body_inertia: float
    body_a: float
    body_b: float
    FWS_k: float
    FWD_c: float
    RWS_k: float
    RWD_c: float
    FWP_theta: float
    FWP_z: float
    RWP_theta: float
    RWP_z: float
    m_wf: float            # front unsprung mass [kg]
    m_wr: float            # rear unsprung mass [kg]
    k_tf: float            # front tyre stiffness [N/m]
    k_tr: float            # rear tyre stiffness [N/m]


@dataclass
class SimulationOptions:
    """Integrator options and initial conditions."""
    t_span: tuple  # (t0, t_final) seconds
    y_0: list      # [z_0, theta_0, z_dot_0, theta_dot_0]
    r_tolerance: float = 1e-7
    a_tolerance: float = 1e-9
    dense: bool = True  # return continuous (dense) solution


# Road input signature:
#  returns (y_f, y_r, y_f_dot, y_r_dot) at time t
BaseInput = Callable[[float], Tuple[float, float, float, float]]

## ======== Start of CSV Loading and interpretation==============
# loading the CSV file and defining the columns
df = pd.read_csv("bumpy_road_cords.csv")
x = df['distance'].values
y = df['height'].values
# defining the spline action with smoothing factor
spline = UnivariateSpline(x, y, s = 0.4)  # adjust s as needed, high s is smoother but less true
#y_smooth is the set of smoothed y values
y_smooth = spline(x)
def zero_base(_: float) -> Tuple[float, float, float, float]:
    x_r = 39.6
    x_f = 39.6 + 2.5
    y_r, y_f, y_r_dot, y_f_dot = spline(x_r) , spline(x_f), spline.derivative()(x_r), spline.derivative()(x_f)
    return y_f, y_r, y_f_dot, y_r_dot # (y_f, y_r, y_f_dot, y_r_dot)

## ========= End of CSV Loading and Interpretation===========

def build_matrices_mck(p: CarParams):
    """
    Assemble mass (M), damping (C), stiffness (K) matrices for the 2-DOF body.
    DOFs are [z, theta].
    """
    m, I = p.body_M, p.body_inertia
    a, b = p.body_a, p.body_b
    k_f, k_r = p.FWS_k, p.RWS_k
    c_f, c_r = p.FWD_c, p.RWD_c

    M = np.array([[m, 0.0],
                  [0.0, I]], dtype = float)

    C = np.array([[c_f + c_r, a*c_f - b*c_r],
                  [a*c_f - b*c_r, a*a*c_f + b*b*c_r]], dtype = float)

    K = np.array([[k_f + k_r, a*k_f - b*k_r],
                  [a*k_f - b*k_r, a*a*k_f + b*b*k_r]], dtype = float)

    return M, C, K

# This function is defining the right hand side of our ODE x_dot = f(t,x)
# We are computing [ż, θ, z̈, θ̈] to feed into our ODE
def rhs_car(t, x, p: CarParams, base: BaseInput):
    """
    State x = [ z, theta, z_wf, z_wr,  z_dot, theta_dot, z_wf_dot, z_wr_dot ].
    Equations from force/moment balance on body with front/rear spring-dampers.
    """
    z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot = x
    y_f, y_r, y_f_dot, y_r_dot = base(t)

    a, b = p.body_a, p.body_b
    k_f, k_r = p.FWS_k, p.RWS_k
    c_f, c_r = p.FWD_c, p.RWD_c
    M, I = p.body_M, p.body_inertia
    m_wf, m_wr = p.m_wf, p.m_wr
    k_tf, k_tr = p.k_tf, p.k_tr

    # Relative spring/damper deflections at front/rear (small-angle)
    dL_f = (z + a*theta) - z_wf
    dL_f_dot = (z_dot + a*theta_dot) - z_wf_dot
    dL_r = (z - b*theta) - z_wr
    dL_r_dot = (z_dot - b*theta_dot) - z_wr_dot

    # Spring-damper forces on body (taking positive "up" in dL, then restoring sign in EOM)
    F_s_f = k_f * dL_f + c_f * dL_f_dot
    F_s_r = k_r * dL_r + c_r * dL_r_dot

    F_t_f = k_tf * (y_f - z_wf)
    F_t_r = k_tr * (y_r - z_wr)

    # Rigid-body equations (heave & pitch)
    z_dot_dot = -(F_s_f + F_s_r) / M # Negative sign due to restoring force driving body down
    theta_dot_dot = -(a*F_s_f - b*F_s_r) / I # Again restoring moments causes the negative sign

    z_wf_dot_dot  = (F_s_f + F_t_f) / m_wf
    z_wr_dot_dot  = (F_s_r + F_t_r) / m_wr

    return np.array([
        z_dot,          # d/dt z
        theta_dot,      # d/dt theta
        z_wf_dot,       # d/dt z_wf
        z_wr_dot,       # d/dt z_wr
        z_dot_dot,         # d/dt z_dot
        theta_dot_dot,     # d/dt theta_dot
        z_wf_dot_dot,      # d/dt z_wf_dot
        z_wr_dot_dot       # d/dt z_wr_dot
    ])


# Key Modelling assumptions for now:
# Small angles so vertical motion at pickups is z + a*theta
# Linear springs/dampers, connected directly to the body (no wheel unsprung mass yet)
# Linearisation about static equilibrium (gravity cancelled by static spring pre-load)

def run_simulation(p: CarParams, base: BaseInput, opts: SimulationOptions):
    """Integrate the ODE using solve_ivp and return the solution object."""
    def fun(t, x):
        return rhs_car(t, x, p, base)

    # rtol / atol are local error control. Roughly, solver tries to keep the local truncation error below atol + rtol * abs(y)
    # dense_output = True builds a continuous time interval
    sol = solve_ivp(fun, opts.t_span, opts.y_0,
                    rtol = opts.r_tolerance, atol = opts.a_tolerance,
                    dense_output = opts.dense)
    return sol

# This function was done by GPT
# This is getting the Eigenvalues of K*phi = lambda*M*phi and converting in Hz to find the resonant frequency of the system as a sanity check
def undamped_naturals(p: CarParams):
    """Return undamped natural frequencies (Hz) of the 2-DOF body model."""
    M, _, K = build_matrices_mck(p)
    lam, _ = np.linalg.eig(np.linalg.solve(M, K))
    wn = np.sqrt(np.clip(lam, 0.0, None))  # rad/s, clips negatives to 0 before square rooting them
    return np.sort(wn / (2*np.pi))         # Hz

# Spacing out uniformly the time and evaluates the dense interpolant sol.sol
def sample_states(sol, t_span, n=2000):
    t = np.linspace(t_span[0], t_span[1], n)
    z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot = sol.sol(t)
    return t, z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot

# In this function we are computing the accelerations using the values from the RHS.
def accelerations_from_rhs(t_arr,
                            z_arr, theta_arr, z_wf_arr, z_wr_arr,
                            z_dot_arr, theta_dot_arr, z_wf_dot_arr, z_wr_dot_arr,
                            p, base):
    """
    Compute accelerations [z̈, θ̈, z̈_wf, z̈_wr] for the 4-DOF half-car model.

    Inputs are the time histories of positions and velocities for body and wheels.
    """

    z_dot_dot = []
    theta_dot_dot = []
    zwf_dot_dot = []
    zwr_dot_dot = []

    for ti, zi, thi, zwfi, zwri, zdi, thdi, zwfdi, zwrdi in zip(
            t_arr,
            z_arr, theta_arr, z_wf_arr, z_wr_arr,
            z_dot_arr, theta_dot_arr, z_wf_dot_arr, z_wr_dot_arr
        ):
        x_state = np.array([
            zi, thi, zwfi, zwri,
            zdi, thdi, zwfdi, zwrdi
        ], dtype=float)

        xdot = rhs_car(ti, x_state, p, base)
        # xdot indices:
        # 4 -> z̈
        # 5 -> θ̈
        # 6 -> z̈_wf
        # 7 -> z̈_wr

        z_dot_dot.append(xdot[4])
        theta_dot_dot.append(xdot[5])
        zwf_dot_dot.append(xdot[6])
        zwr_dot_dot.append(xdot[7])

    return (np.array(z_dot_dot),
            np.array(theta_dot_dot),
            np.array(zwf_dot_dot),
            np.array(zwr_dot_dot))


# Computing accelerations using the sol.sol function
def accelerations_from_sol(sol, p: CarParams, base: BaseInput, t_grid):
    """Same as above but accepts a solution + time grid."""
    z_dot_dot, theta_dot_dot = [], []
    for ti in t_grid:
        z, theta, z_dot, theta_dot = sol.sol(ti)
        _, _, z_dot_doti, theta_dot_doti = rhs_car(ti, [z, theta, z_dot, theta_dot], p, base)
        z_dot_dot.append(z_dot_doti); theta_dot_dot.append(theta_dot_doti)
    return np.array(z_dot_dot), np.array(theta_dot_dot)


def occupant_vertical_accel(z_dot_dot, theta_dot_dot, x_from_CG = 0.0):
    """Vertical acceleration at a point x ahead of CG (x > 0 forward)."""
    # For small angles (small pitch), the vertical component of acceleration becomes z_dot_dot + x*theta_dot_dot
    return z_dot_dot + x_from_CG * theta_dot_dot

# Computes comfort metric, RMS of vertical accel at CG and at seat location
def rms(a):
    """Root-mean-square of a signal."""
    return float(np.sqrt(np.mean(a**2)))

#Parameters
# Nominal average car values
p = CarParams(
    body_M = 1200.0,          # Sprung mass (kg)
    body_inertia = 2200.0,    # Pitch inertia (kg·m²)
    body_a = 1.2,             # CG to front axle (m)
    body_b = 1.3,             # CG to rear axle (m)

    FWS_k = 35e3,             # Front suspension spring (N/m)
    FWD_c = 3.0e3,            # Front damper (N·s/m)
    RWS_k = 30e3,             # Rear suspension spring (N/m)
    RWD_c = 3.0e3,            # Rear damper (N·s/m)

    m_wf = 40.0,              # Front unsprung mass (kg)
    m_wr = 35.0,              # Rear unsprung mass (kg)
    k_tf = 2.0e5,             # Front tyre stiffness (N/m)
    k_tr = 1.8e5,             # Rear tyre stiffness (N/m)

    FWP_theta = 0.0,
    FWP_z = 0.0,
    RWP_theta = 0.0,
    RWP_z = 0.0
)


# Simulation options
# 5 second test, release the car from a 2cm free fall (wheels are touching the ground, spring is just stretched)
opts = SimulationOptions(
    t_span = (0.0, 5.0),      # 5 s
    y_0 = [
    0.02,        # z
    0.0,         # theta
    0.0,         # z_wf
    0.0,         # z_wr
    0.0,         # z_dot
    0.0,         # theta_dot
    0.0,         # z_wf_dot
    0.0          # z_wr_dot
]
)

# Run
# Integrates system with flat road
sol = run_simulation(p, zero_base, opts)

# Sample states uniformly
# ts is the continuous times, whereas the rest of the values are sampled from the sol.sol(ts) function which are the model's solutions from our equations
ts, z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot = sample_states(sol, opts.t_span, n = 2000)

# Accelerations (via RHS)
# Finding the accelerations via our formula (the physics)
z_dot_dot, theta_dot_dot, z_wf_dot_dot, z_wr_dot_dot = accelerations_from_rhs(
    ts,
    z, theta, z_wf, z_wr,
    z_dot, theta_dot, z_wf_dot, z_wr_dot,
    p,
    zero_base
)

# Occupant acceleration at a chosen x from CG
seat_x = 0.8
a_pass = occupant_vertical_accel(z_dot_dot, theta_dot_dot, x_from_CG = seat_x)

# Metrics
g = 9.81
rms_z_dot_dot = rms(z_dot_dot)
rms_pass = rms(a_pass)
print(f"RMS heave accel at CG: {rms_z_dot_dot:.4f} m/s^2  ({rms_z_dot_dot/g:.4f} g)")
print(f"RMS vertical accel at passenger (x={seat_x} m): {rms_pass:.4f} m/s^2  ({rms_pass/g:.4f} g)")

# --- Plots (one per figure; no specific colors) ---
# (A) Heave displacement
plt.figure()
plt.plot(ts, z)
plt.xlabel("Time [s]")
plt.ylabel("Heave z [m]")
plt.title("Body Heave vs Time")
plt.grid(True)

# (B) Pitch angle (deg)
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
plt.title(f"Passenger Vertical Accel (x={seat_x} m from CG)")
plt.grid(True)

# (D) Wheel accel
plt.figure()
plt.plot(ts, z_wf_dot_dot, label = 'Front Wheel')
plt.plot(ts, z_wr_dot_dot, label = 'Rear Wheel')
plt.xlabel("Time [s]")
plt.ylabel("Wheel acceleration [m/s²]")
plt.title("Unsprung (wheel) accelerations")
plt.grid(True)
plt.legend()


plt.show()


