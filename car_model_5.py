# -*- coding: utf-8 -*-

# 2D car (4-DOF) suspension model:
# - Sprung mass (body): heave z and pitch theta
# - Two unsprung masses: front and rear wheels (z_wf, z_wr are the vertical positions)
# - Linear springs/dampers in suspension, linear tyre stiffness
# - Car drives along a measured road profile (from CSV)

from dataclasses import dataclass # Makes classes easier to read and write removing _init_ and parameter confucion
from typing import Callable, Tuple

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import UnivariateSpline

@dataclass
class CarParams:

    # Car Body
    body_M: float # Mass (kg)
    body_inertia: float # Inertia about CG (kgm^2)
    body_a: float # Distance CG to front suspension (m)
    body_b: float # Distance CG to back suspension (m)

    # Body Suspension
    FWS_k: float # Front spring (N/m)
    FWD_c: float # Front damper (Ns/m)
    RWS_k: float # Rear spring (N/m)
    RWD_c: float # Rear damper (Ns/m)

    # Wheel Masses
    m_wf: float # Front wheel mass (kg)
    m_wr: float # Rear wheel mass (kg)

    # Tyre stiffness
    k_tf: float # Front tyre (N/m)
    k_tr: float # Rear tyre (N/m)


@dataclass
class SimulationOptions:
    t_span: tuple # t0 to t_end (s) start to end of simulation (length of road in our case)
    y_0: list # Initial state, just the position of car at beginning - will be 0
    r_tolerance: float = 1e-7 # Tolerances required for solving later
    a_tolerance: float = 1e-9
    dense: bool = True # Dense solution again required for solving later


# THIS NEXT SECTION OF THE CODE IS EXPLAINED IN THE "IMPORTING THE ROAD PROFILE" SHEET OF NOTES


# Road profile, the CSV has column "distance" (m) and "height" (m)
BaseInput = Callable[[float], Tuple[float, float, float, float]]

# make_road_base function outputs the base(t) function which returns (yf, yr, yfdot, yrdot) as per the next bits of code
def make_road_base(p: CarParams,
                   spline: UnivariateSpline,
                   dsdx: UnivariateSpline,
                   x: np.ndarray,
                   v: float = 8.0,
                   x0: float = 0.0,
                   clamp: bool = True) -> BaseInput: # The clamp clips the car to act within the maximum and minimum values of the road
    
    # Minimum and maximum values for the road
    xmin, xmax = float(np.min(x)), float(np.max(x))

    def base(t: float) -> Tuple[float, float, float, float]:
        # Body CG position along the road
        x_cg = x0 + v * t

        # Front and rear suspension locations
        x_f = x_cg + p.body_a
        x_r = x_cg - p.body_b

        # If clamp is true (it is) the front wheel is bounded to stay within xmax, the rear wheel is bounded to stay within xmin
        if clamp:
            x_fq = np.clip(x_f, xmin, xmax)
            x_rq = np.clip(x_r, xmin, xmax)
        else:
            x_fq, x_rq = x_f, x_r

        # Road height under each wheel
        y_f = float(spline(x_fq))
        y_r = float(spline(x_rq))

        # dy/dt = (dy/dx) * v
        y_f_dot = float(dsdx(x_fq)) * v
        y_r_dot = float(dsdx(x_rq)) * v

        return y_f, y_r, y_f_dot, y_r_dot

    return base


# THIS NEXT SECTION OF THE CODE IS EXPLAINED IN THE "BUILDING THE MATRICES AND THE UNDAMPED NATURALS" SHEET OF NOTES


# These matrices are built around the original simplied 2-DOF model to later check for resonant frequencies as a sanity check
def build_matrices_mck(p: CarParams):

    # Necessary variables to build matrices imported from the CarParams class
    m, I = p.body_M, p.body_inertia
    a, b = p.body_a, p.body_b
    k_f, k_r = p.FWS_k, p.RWS_k
    c_f, c_r = p.FWD_c, p.RWD_c

    # Mass matrix
    M = np.array([[m, 0.0],
                  [0.0, I]], dtype = float)

    # Damping matrix
    C = np.array([[c_f + c_r, a * c_f - b * c_r],
                  [a * c_f - b * c_r, a * a * c_f + b * b * c_r]], dtype = float)

    # Stiffness matrix
    K = np.array([[k_f + k_r, a * k_f - b * k_r],
                  [a * k_f - b * k_r, a * a * k_f + b * b * k_r]], dtype = float)

    return M, C, K


# THIS NEXT SECTION OF THE CODE IS EXPLAINED IN THE "DAMPING RATIOS FOR 4 DOFS" SHEET OF NOTES


# This function returns the frequencies in hz for the 2 modes (bounce and pitch)
def modal_damping_ratios(p: CarParams):

    # Mass, damping and stiffness matrices
    M, C, K = build_matrices_mck(p)

    # Solve the eigenproblem from before
    # By forming A = M^{-1} K and finding its eigenvalues and eigenvectors
    A = np.linalg.solve(M, K)
    eigvals, eigvecs = np.linalg.eig(A)

    # Natural frequencies in rad/s
    wn = np.sqrt(np.clip(eigvals, 0.0, None))

    # Sort from lowest to highest frequency
    order = np.argsort(wn)
    wn = wn[order]
    eigvecs = eigvecs[:, order]  # Columns: bounce, pitch

    # Create array for damping ratios
    zetas = np.zeros_like(wn, dtype = float)

    # Iterate over the 2 modes (bounce and pitch)
    for i in range(len(wn)):
        phi = eigvecs[:, i]  # Mode shape for this mode (length-2 array)

        # Modal mass: m_i = phi^T M phi
        m_i = float(phi.T @ M @ phi)

        # Modal stiffness: k_i = phi^T K phi
        k_i = float(phi.T @ K @ phi)

        # Modal damping: c_i = phi^T C phi
        c_i = float(phi.T @ C @ phi)

        # Damping ratio for this mode: dr_i = zeta_i / (2 * sqrt(m_i * k_i))
        zetas[i] = c_i / (2.0 * np.sqrt(m_i * k_i))

    # Convert frequencies to Hz for reporting
    freqs_hz = wn / (2.0 * np.pi)

    return freqs_hz, zetas


# THIS NEXT SECTION OF THE CODE IS EXPLAINED IN THE "THE EQUATIONS OF MOTION" SHEET OF NOTES


def rhs_car(t: float, x: np.ndarray, p: CarParams, base: BaseInput) -> np.ndarray:
    # Assign initial state and input variables
    z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot = x

    # Road input under the tyres, generated by CSV and converted to y(t) by make_road_base function
    y_f, y_r, y_f_dot, y_r_dot = base(t)

    # Conversion of necessary paramaters from CarParams into readable ones
    a, b = p.body_a, p.body_b
    k_f, k_r = p.FWS_k, p.RWS_k
    c_f, c_r = p.FWD_c, p.RWD_c
    M, I = p.body_M, p.body_inertia
    m_wf, m_wr = p.m_wf, p.m_wr
    k_tf, k_tr = p.k_tf, p.k_tr

    # Body suspension deflections– small angle approx:
    # Front pickup: z + a*theta, Rear pickup: z - b*theta
    dL_f = (z + a * theta) - z_wf
    dL_f_dot = (z_dot + a * theta_dot) - z_wf_dot

    dL_r = (z - b * theta) - z_wr
    dL_r_dot = (z_dot - b * theta_dot) - z_wr_dot

    # Suspension forces as per maths
    F_s_f = k_f * dL_f + c_f * dL_f_dot
    F_s_r = k_r * dL_r + c_r * dL_r_dot

    # Tyre forces as per maths
    F_t_f = k_tf * (y_f - z_wf)
    F_t_r = k_tr * (y_r - z_wr)

    # Body equations of motion
    z_dot_dot = -(F_s_f + F_s_r) / M
    theta_dot_dot = -(a * F_s_f - b * F_s_r) / I

    # Wheel equations of motion
    z_wf_dot_dot = (F_s_f + F_t_f) / m_wf
    z_wr_dot_dot = (F_s_r + F_t_r) / m_wr

    return np.array([
        z_dot,
        theta_dot,
        z_wf_dot,
        z_wr_dot,
        z_dot_dot,
        theta_dot_dot,
        z_wf_dot_dot,
        z_wr_dot_dot
    ], dtype = float)


# THIS NEXT SECTION OF THE CODE IS EXPLAINED IN THE "SIMULATION HELPERS" SHEET OF NOTES


# This function is wrapping the solver so that it can be called later
def run_simulation(p: CarParams, base: BaseInput, opts: SimulationOptions):
    def fun(t, x):
        return rhs_car(t, x, p, base)

    sol = solve_ivp(fun,
                    opts.t_span,
                    opts.y_0,
                    rtol = opts.r_tolerance,
                    atol = opts.a_tolerance,
                    dense_output = opts.dense)
    return sol


# This function distributes the x values along a uniform time grid
def sample_states(sol, t_span, n = 2000):
    ts = np.linspace(t_span[0], t_span[1], n)
    z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot = sol.sol(ts)
    return ts, z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot

# Calls rhs_car (equations of motion), to achieve accelerations from stored positions / velocities
def accelerations_from_rhs(t_arr,
                           z_arr, theta_arr, z_wf_arr, z_wr_arr,
                           z_dot_arr, theta_dot_arr, z_wf_dot_arr, z_wr_dot_arr,
                           p: CarParams,
                           base: BaseInput): # Function takes in time arrays for each input

    n = len(t_arr)

    # Pre-creating arrays with correct number of slots based on the length of time in t_arr
    z_dot_dot = np.zeros(n)
    theta_dot_dot = np.zeros(n)
    z_wf_dot_dot = np.zeros(n)
    z_wr_dot_dot = np.zeros(n)

    # Iterating through time
    for i in range(n):
        ti = t_arr[i]

        # The state of the input positions and velocities at this time step
        x_state = np.array([
            z_arr[i],
            theta_arr[i],
            z_wf_arr[i],
            z_wr_arr[i],
            z_dot_arr[i],
            theta_dot_arr[i],
            z_wf_dot_arr[i],
            z_wr_dot_arr[i]
        ], float)

        # Returning differentiated values (accelerations) from rhs_car (equations of motion) at each time step
        x_dot = rhs_car(ti, x_state, p, base)

        # From the x_dot, pick out all the relevant accelerations to then be used in graphing and later parts
        z_dot_dot[i] = x_dot[4]
        theta_dot_dot[i] = x_dot[5]
        z_wf_dot_dot[i] = x_dot[6]
        z_wr_dot_dot[i] = x_dot[7]

    return z_dot_dot, theta_dot_dot, z_wf_dot_dot, z_wr_dot_dot

# Vertical acceleration of passenger at position x from CG
def occupant_vertical_accel(z_dot_dot, theta_dot_dot, x_from_CG: float = 0.0):
    return z_dot_dot + x_from_CG * theta_dot_dot

# Root mean square speed of a variable
def rms(a):
    return float(np.sqrt(np.mean(a**2)))


# THIS NEXT SECTION OF THE CODE IS EXPLAINED IN THE "ROOT FINDING (BISECTION METHOD)" SHEET OF NOTES


# This function uses the root finding bisection method to find root of a function f(x) = 0 between a,b
def bisection(f, a, b, N): # N is the number of iterations

    fa = f(a)
    fb = f(b)

    # Very simple algorithm, explained in the maths notes
    if fa * fb >= 0:
        raise ValueError("Bisection requires f(a)*f(b) < 0")

    for _ in range(N):
        m = 0.5 * (a + b)
        fm = f(m)

        if fa * fm < 0:
            b = m
            fb = fm
        else:
            a = m
            fa = fm

    return 0.5 * (a + b)


# THIS NEXT SECTION OF THE CODE IS EXPLAINED IN THE "OPTIMISATION (FINDING OPTIMAL DAMPING)" SHEET OF NOTES


# This function returns the RMS acceleration for a passenger
def passenger_rms_for_damping(c_f, c_r, p, base, opts_local, seat_x):

    # Set parameters to our CarParams except for the damping which we are inputting ourselves
    pars = CarParams(
        body_M = p.body_M,
        body_inertia = p.body_inertia,
        body_a = p.body_a,
        body_b = p.body_b,
        FWS_k = p.FWS_k,
        FWD_c = c_f,
        RWS_k = p.RWS_k,
        RWD_c = c_r,
        m_wf = p.m_wf,
        m_wr = p.m_wr,
        k_tf = p.k_tf,
        k_tr = p.k_tr
    )

    # Run the solver to uniformly distribut positions and velocities later with sample_states
    sol = run_simulation(pars, base, opts_local)
    ts, z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot = sample_states(sol, opts_local.t_span, n = 300)

    # Find car body accelerations from the function to use in finding the acceleration of the passenger distance x from CG
    z_dot_dot, theta_dot_dot, _, _ = accelerations_from_rhs(ts, z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot, pars, base)
    a_pass = occupant_vertical_accel(z_dot_dot, theta_dot_dot, x_from_CG=seat_x)
    # Return the RMS of the passenger accel
    return rms(a_pass)

# The 2 functions below find the differentiated J(cf) and J(cr)

def dJ_dc_f(c_f, c_r, p, base, opts_local, seat_x):
    # Small step for differentiation (10% of damping coeff)
    h = 0.1 * c_f if c_f > 0 else 1.0

    J_plus = passenger_rms_for_damping(c_f + h, c_r, p, base, opts_local, seat_x)
    J_minus = passenger_rms_for_damping(c_f - h, c_r, p, base, opts_local, seat_x)

    # Returns the differentiated J(cf)
    return (J_plus - J_minus) / (2 * h)


def dJ_dc_r(c_f, c_r, p, base, opts_local, seat_x):
    # Small step for differentiation (10% of damping coeff)
    h = 0.1 * c_r if c_r > 0 else 1.0

    J_plus = passenger_rms_for_damping(c_f, c_r + h, p, base, opts_local, seat_x)
    J_minus = passenger_rms_for_damping(c_f, c_r - h, p, base, opts_local, seat_x)

    # Returns the differentiated J(cf)
    return (J_plus - J_minus) / (2 * h)


# This function finds an optimal c in the case that root finding fails 
# Root finding can fail in a road such as a motorway as there is not "too high" or "too low" c as it is flat
# It will just pick a value of derivative that is closest to 0 
def pick_optimal_c(f, a, b, N):

    fa = f(a)
    fb = f(b)

    # If there is a sign change, we can safely use bisection method (as before)
    if fa * fb < 0.0:
        return bisection(f, a, b, N)

    # If either endpoint already gives zero derivative, take it and use it
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b

    # Otherwise, derivative has the same sign on [a, b]:
    # If derivative > 0 everywhere, J increases with c -> minimum at a
    # If derivative < 0 everywhere, J decreases with c -> minimum at b
    if fa > 0.0 and fb > 0.0:
        # J increases with c, smallest J at lower bound
        return a
    if fa < 0.0 and fb < 0.0:
        # J decreases with c, smallest J at upper bound
        return b

    # Fallback (should not really be hit, but just in case of numerical quirks)
    return a if abs(fa) < abs(fb) else b


# This function finds optimal cf and cr using the bisection root finding methods on the derivates of J(cf, cr)
def optimise_damping(p, base, opts_local, seat_x,
                               c_f_range = (100, 3000),
                               c_r_range = (100, 3000),
                               N = 1):

    # Local convergence history for 1 optimisation
    c_f_opt_seq = []
    c_r_opt_seq = []

    # Take initial guesses damping coefficients
    c_f_opt = 1000.0
    c_r_opt = 1000.0

    for G in range(1, 9):

        # Optimise the front damping (fix rear)
        f_front = lambda c: dJ_dc_f(c, c_r_opt, p, base, opts_local, seat_x)
        a_f, b_f = c_f_range
        c_f_opt = pick_optimal_c(f_front, a_f, b_f, G)

        # Optimise the rear damping (fix front)
        f_rear = lambda c: dJ_dc_r(c_f_opt, c, p, base, opts_local, seat_x)
        a_r, b_r = c_r_range
        c_r_opt = pick_optimal_c(f_rear, a_r, b_r, G)

        c_f_opt_seq.append(c_f_opt)
        c_r_opt_seq.append(c_r_opt)

    # Return optimised damping coefficient values
    return c_f_opt, c_r_opt, np.array(c_f_opt_seq), np.array(c_r_opt_seq)


# THIS NEXT SECTION OF THE CODE IS EXPLAINED IN THE "OPTIMISATION (FINDING OPTIMAL DAMPING)" SHEET OF NOTES


# The Monte Carlo error approximation is a function I wrote as a personal project some time ago
# It is a method of estimating error when the true value is unknown
def damping_monte_carlo_error(p, base, opts_local, seat_x,
                                 n_samples = 4, rel_variation = 0.05):

    c_f_opt_nom, c_r_opt_nom, c_f_list, c_r_list = optimise_damping(p, base, opts_local, seat_x)

    # Start by taking N samples of the values you are trying to find
    c_f_samples = []
    c_r_samples = []

    for _ in range(n_samples):

        # Random +-5% changes of mass and spring stiffness
        scale_M = np.random.uniform(1 - rel_variation, 1 + rel_variation)
        scale_kf = np.random.uniform(1 - rel_variation, 1 + rel_variation)
        scale_kr = np.random.uniform(1 - rel_variation, 1 + rel_variation)

        pars = CarParams(
            body_M = p.body_M * scale_M,
            body_inertia = p.body_inertia,
            body_a = p.body_a,
            body_b = p.body_b,
            FWS_k = p.FWS_k * scale_kf,
            FWD_c = p.FWD_c,
            RWS_k = p.RWS_k * scale_kr,
            RWD_c = p.RWD_c,
            m_wf = p.m_wf,
            m_wr = p.m_wr,
            k_tf = p.k_tf,
            k_tr = p.k_tr
        )

        # Find the new optimum for this random setup
        c_f_ran, c_r_ran, _, _ = optimise_damping(pars, base, opts_local, seat_x)
        c_f_samples.append(c_f_ran)
        c_r_samples.append(c_r_ran)

    # Convert to arrays and compute Standard errors
    c_f_samples = np.array(c_f_samples)
    c_r_samples = np.array(c_r_samples)

    cf_mean = np.mean(c_f_samples)
    cr_mean = np.mean(c_r_samples)

    s2_cf = np.sum((c_f_samples - cf_mean)**2) / (n_samples - 1)
    s2_cr = np.sum((c_r_samples - cr_mean)**2) / (n_samples - 1)

    se_cf = np.sqrt(s2_cf / n_samples)
    se_cr = np.sqrt(s2_cr / n_samples)

    return c_f_opt_nom, c_r_opt_nom, se_cf, se_cr, c_f_list, c_r_list


# THIS NEXT SECTION OF THE CODE IS EXPLAINED IN THE "MAIN SCRIPT" SHEET OF NOTES


def main(p, road_base):

    # Initial conditions: body height, wheels on the road, zero velocity
    y_f0, y_r0, _, _ = road_base(0.0)
    y0 = [0.0, 0.0, y_f0, y_r0, 0.0, 0.0, 0.0, 0.0]

    # Solver inputs, set tolerances, time span (20s to match velocity of car)
    opts = SimulationOptions(t_span = (0.0, 20.0), 
                            y_0 = y0)
    opts_opt = SimulationOptions(t_span = (0.0, 8.0), 
                                y_0 = y0, 
                                r_tolerance = 1e-4, 
                                a_tolerance = 1e-6,
                                dense = True)

    # Passenger position (m) forward of CG
    seat_x = 0.8

    # Optimisation (root finding with bisection) + Monte Carlo error to give optimised damping coefficient and error
    c_f_opt, c_r_opt, se_cf, se_cr, c_f_opt_new, c_r_opt_new = damping_monte_carlo_error(
    p, road_base, opts_opt, seat_x,
    n_samples = 4, rel_variation = 0.05
    )

    # Update parameters to the optimal values
    p.FWD_c = c_f_opt
    p.RWD_c = c_r_opt

    # Find damping ratios and freqs
    freqs_hz, zetas = modal_damping_ratios(p)

    print(f"Optimal front damping  FWD_c = {c_f_opt:.1f} Ns/m")
    print(f"Optimal rear damping RWD_c = {c_r_opt:.1f} Ns/m")
    print("Body modal properties (with optimised damping):")
    print(f"Mode 1 (bounce-ish): f = {freqs_hz[0]:.2f} Hz, damping ratio = {zetas[0]:.3f}")
    print(f"Mode 2 (pitch-ish) : f = {freqs_hz[1]:.2f} Hz, damping ratio = {zetas[1]:.3f}")
    print(f"Std. error front = {se_cf:.2f} rear = {se_cr:.2f} Ns/m")
    print(f"95% CI front: {c_f_opt:.1f} +- {1.96*se_cf:.2f}")
    print(f"95% CI rear : {c_r_opt:.1f} +- {1.96*se_cr:.2f}\n")

    # Final simulation using optimal damping
    sol = run_simulation(p, road_base, opts)
    ts, z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot = sample_states(sol, opts.t_span, n = 2000)
    z_dot_dot, theta_dot_dot, z_wf_dot_dot, z_wr_dot_dot = accelerations_from_rhs(ts, z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot, p, road_base)
    a_pass = occupant_vertical_accel(z_dot_dot, theta_dot_dot, x_from_CG = seat_x)

    # Ride comfort metrics
    g = 9.81
    rms_z = rms(z_dot_dot)
    rms_pass = rms(a_pass)

    print(f"RMS body heave accel: {rms_z:.3f} m/s² ({rms_z/g:.3f} g)")
    print(f"RMS passenger accel : {rms_pass:.3f} m/s² ({rms_pass/g:.3f} g)")

    # PLOTS
    # Body heave
    def plot_heave(ts, z):
        plt.figure()
        plt.plot(ts, z, lw = 1.2)
        plt.xlabel("Time (s)")
        plt.ylabel("Heave z (m)")
        plt.title("Body heave (optimal damping)")
        plt.grid(True)

    # Body pitch
    def plot_pitch(ts, theta):
        plt.figure()
        plt.plot(ts, np.degrees(theta), lw = 1.2)
        plt.xlabel("Time (s)")
        plt.ylabel("Pitch θ (deg)")
        plt.title("Body pitch (optimal damping)")
        plt.grid(True)

    # Passenger vertical acceleration
    def plot_passenger_accel(ts, a_pass, seat_x):
        plt.figure()
        plt.plot(ts, a_pass, lw = 1.2)
        plt.xlabel("Time (s)")
        plt.ylabel("Vertical acceleration (m/s^2)")
        plt.title(f"Passenger vertical acceleration (x = {seat_x} m from CG)")
        plt.grid(True)

    def plot_pass_accel_conv(c_f_opt_new, c_r_opt_new):
        plt.figure(figsize = (10, 5))

        for i, (c_f, c_r) in enumerate(zip(c_f_opt_new, c_r_opt_new), start = 1):

            # Set bopth dampers in CarParams to the index i in the list
            p.FWD_c = c_f
            p.RWD_c = c_r

            # Run the simulation each time for each c
            sol = run_simulation(p, road_base, opts)
            ts, z, theta, z_wf, z_wr, z_dot, theta_dot, z_wf_dot, z_wr_dot = sample_states(sol, opts.t_span, n = 2000)
            z_dot_dot, theta_dot_dot, _, _ = accelerations_from_rhs(ts, z, theta, z_wf, z_wr,
                                                z_dot, theta_dot, z_wf_dot, z_wr_dot, p, road_base)
            a_pass = occupant_vertical_accel(z_dot_dot, theta_dot_dot, x_from_CG = seat_x)

            # Places the last iteration on top plus a slightly bolder alpha to make it clearer in graph
            plt.plot(ts, a_pass, lw = 1.1,
                    alpha = 1.0 if i == len(c_f_opt_new) else 0.8,
                    label = f"Step {i}: c_f={c_f:.0f}, c_r={c_r:.0f}  |  RMS={rms(a_pass):.3f} m/s^2")

        # Styling
        plt.axhline(0, ls = '--', lw = 1)
        plt.grid(True, alpha = 0.3)
        plt.xlabel("Time (s)")
        plt.ylabel("Passenger vertical acceleration (m/s^2)")
        plt.title("Passenger acceleration convergence")
        plt.legend(loc = "upper right", fontsize = 8, frameon = True)
        plt.tight_layout()

    plot_heave(ts, z)
    plot_pitch(ts, theta)
    plot_passenger_accel(ts, a_pass, seat_x,)
    plot_pass_accel_conv(c_f_opt_new, c_r_opt_new)