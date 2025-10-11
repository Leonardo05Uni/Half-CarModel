# -*- coding: utf-8 -*-

# ----------------------------
# Imports
# ----------------------------
from dataclasses import dataclass
from typing import Callable, Tuple
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# ----------------------------
# Data containers
# ----------------------------
@dataclass
class CarParams:
    """All physical parameters of the 2-DOF car body model."""
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


@dataclass
class SimulationOptions:
    """Integrator options and initial conditions."""
    t_span: tuple  # (t0, t_final) seconds
    y_0: list      # [z_0, theta_0, z_dot_0, theta_dot_0]
    r_tolerance: float = 1e-7
    a_tolerance: float = 1e-9
    dense: bool = True  # return continuous (dense) solution


# ----------------------------
# Types & base (road) input
# ----------------------------
# Road input signature:
#  returns (y_f, y_r, y_f_dot, y_r_dot) at time t
BaseInput = Callable[[float], Tuple[float, float, float, float]]

def zero_base(_: float) -> Tuple[float, float, float, float]:
    """Flat road: zero displacement/velocity at both wheels."""
    return 0.0, 0.0, 0.0, 0.0 # (y_f, y_r, y_f_dot, y_r_dot)


# ----------------------------
# System matrices (for modal checks)
# ----------------------------
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
                  [0.0, I]], dtype=float)

    C = np.array([[c_f + c_r, a*c_f - b*c_r],
                  [a*c_f - b*c_r, a*a*c_f + b*b*c_r]], dtype=float)

    K = np.array([[k_f + k_r, a*k_f - b*k_r],
                  [a*k_f - b*k_r, a*a*k_f + b*b*k_r]], dtype=float)

    return M, C, K


# ----------------------------
# Dynamics (RHS)
# ----------------------------
# This function is defining the right hand side of our ODE x_dot = f(t,x)
# We are computing [ż, θ, z̈, θ̈] to feed into our ODE
def rhs_car(t, x, p: CarParams, base: BaseInput):
    """
    State x = [z, theta, z_dot, theta_dot].
    Equations from force/moment balance on body with front/rear spring-dampers.
    """
    z, theta, z_dot, theta_dot = x
    y_f, y_r, y_f_dot, y_r_dot = base(t)

    a, b = p.body_a, p.body_b
    k_f, k_r = p.FWS_k, p.RWS_k
    c_f, c_r = p.FWD_c, p.RWD_c

    # Relative spring/damper deflections at front/rear (small-angle)
    dL_f = (z + a*theta) - y_f
    dL_r = (z - b*theta) - y_r
    dL_f_dot = (z_dot + a*theta_dot) - y_f_dot # The time derivatives are needed for viscous damper forces later
    dL_r_dot = (z_dot - b*theta_dot) - y_r_dot

    # Spring-damper forces on body (taking positive "up" in dL, then restoring sign in EOM)
    F_f = k_f*dL_f + c_f*dL_f_dot
    F_r = k_r*dL_r + c_r*dL_r_dot

    # Rigid-body equations (heave & pitch)
    z_dot_dot      = -(F_f + F_r) / p.body_M # Negative sign due to restoring force driving body down
    theta_dot_dot  = -(a*F_f - b*F_r) / p.body_inertia # Again restoring moments causes the negative sign

    return [z_dot, theta_dot, z_dot_dot, theta_dot_dot]  # Intentionally keep order [ż, θ, z̈, θ̈] to later put in ODE solver solve_ivp


# Key Modelling assumptions for now:
# Small angles so vertical motion at pickups is z + a*theta
# Linear springs/dampers, connected directly to the body (no wheel unsprung mass yet)
# Linearisation about static equilibrium (gravity cancelled by static spring pre-load)

# ----------------------------
# Integrator wrapper
# ----------------------------
def run_simulation(p: CarParams, base: BaseInput, opts: SimulationOptions):
    """Integrate the ODE using solve_ivp and return the solution object."""
    def fun(t, x):
        return rhs_car(t, x, p, base)

    # rtol / atol are local error control. Roughly, solver tries to keep the local truncation error below atol + rtol * abs(y)
    # dense_output = True builds a continuous time interval
    sol = solve_ivp(fun, opts.t_span, opts.y_0,
                    rtol=opts.r_tolerance, atol=opts.a_tolerance,
                    dense_output=opts.dense)
    return sol


# ----------------------------
# Analysis helpers
# ----------------------------

# This function was done by GPT
# This is getting the Eigenvalues of K*phi = lambda*M*phi and converting in Hz to find the resonant frequency of the system as a sanity check
def undamped_naturals(p: CarParams):
    """Return undamped natural frequencies (Hz) of the 2-DOF body model."""
    M, _, K = build_matrices_mck(p)
    lam, _ = np.linalg.eig(np.linalg.solve(M, K))
    wn = np.sqrt(np.clip(lam, 0.0, None))  # rad/s, clips negatives to 0 before square rooting them
    return np.sort(wn / (2*np.pi))         # Hz

# Spacing out uniformly the time and evaluates the dense interpolant sol.sol
def sample_states(sol, t_span, n = 2000):
    """Uniformly sample the dense solution (requires dense_output=True)."""
    t = np.linspace(t_span[0], t_span[1], n)
    z, theta, z_dot, theta_dot = sol.sol(t)
    
    # Returns state components as arrays of length n
    return t, z, theta, z_dot, theta_dot

# In this function we are computing the accelerations using the values from the RHS.
def accelerations_from_rhs(t, z, theta, z_dot, theta_dot, p, base):
    """Compute [z̈, θ̈] by calling the RHS (noise-free vs numerical differentiation)."""
    z_dot_dot, theta_dot_dot = [], []
    for ti, zi, thetai, z_doti, theta_doti in zip(t, z, theta, z_dot, theta_dot):
        _, _, z_dot_doti, theta_dot_doti = rhs_car(ti, [zi, thetai, z_doti, theta_doti], p, base)
        z_dot_dot.append(z_dot_doti); theta_dot_dot.append(theta_dot_doti)
    return np.asarray(z_dot_dot), np.asarray(theta_dot_dot)

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


# ----------------------------
# Script entry point (example use)
# ----------------------------
#Parameters
# Nominal average car values
p = CarParams(
    body_M = 1200.0, body_inertia = 2200.0,
    body_a = 1.2, body_b = 1.3,
    FWS_k = 35e3, FWD_c = 3.0e3,
    RWS_k = 30e3, RWD_c = 3.0e3,
    FWP_theta = 0.0, FWP_z = 0.0,
    RWP_theta = 0.0, RWP_z = 0.0
)

# Simulation options
# 5 second test, release the car from a 2cm free fall (wheels are touching the ground, spring is just stretched)
opts = SimulationOptions(
    t_span = (0.0, 5.0),             # 5 s
    y_0 = [0.02, 0.0, 0.0, 0.0]      # 2 cm heave release, zero pitch & rates
)

# Run
# Integrates system with flat road
sol = run_simulation(p, zero_base, opts)

# Sample states uniformly
# ts is the continuous times, whereas the rest of the values are sampled from the sol.sol(ts) function which are the model's solutions from our equations
ts, z, theta, z_dot, theta_dot = sample_states(sol, opts.t_span, n = 2000)

# Accelerations (via RHS)
# Finding the accelerations via our formula (the physics)
z_dot_dot, theta_dot_dot = accelerations_from_rhs(ts, z, theta, z_dot, theta_dot, p, zero_base)

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

plt.show()
