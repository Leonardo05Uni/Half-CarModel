# -*- coding: utf-8 -*-

# 2D car (4-DOF) suspension model:
#  - Sprung mass (body): heave z and pitch theta
#  - Two unsprung masses: front and rear wheels (z_wf, z_wr are the vertical positions)
#  - Linear springs/dampers in suspension, linear tyre stiffness
#  - Car drives along a measured road profile (from CSV)

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

# Import CSV
df = pd.read_csv("bumpy_road_cords.csv")
x = df["distance"].values # Distance labelled as x
y = df["height"].values # Road height labelled as y

# Smooth road profile y(x) using UnivariateSpline import 
spline = UnivariateSpline(x, y, s = 0.4) # s = 0.4 is used to define how smooth you want the curve to be, higher value more smooth
dsdx = spline.derivative() # dy/dx for road inputs later

# make_road_base function outputs the base(t) function which returns (yf, yr, yfdot, yrdot) as per the next bits of code
def make_road_base(p: CarParams,
                   spline: UnivariateSpline,
                   dsdx: UnivariateSpline,
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

# Finding natural frequencies of body
def undamped_naturals(p: CarParams):
    M, _, K = build_matrices_mck(p)
    lam, _ = np.linalg.eig(np.linalg.solve(M, K))  # Eigenvalues of M^(-1)K
    wn = np.sqrt(np.clip(lam, 0.0, None)) # Finding frequency (rad/s)
    return np.sort(wn / (2 * np.pi)) # Rad/s to Hz


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
