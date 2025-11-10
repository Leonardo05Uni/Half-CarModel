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
    body_M: float          # Mass (kg)
    body_inertia: float    # Inertia about CG (kgm^2)
    body_a: float          # Distance CG to front suspension (m)
    body_b: float          # Distance CG to back suspension (m)

    # Body Suspension
    FWS_k: float           # Front spring (N/m)
    FWD_c: float           # Front damper (Ns/m)
    RWS_k: float           # Rear spring (N/m)
    RWD_c: float           # Rear damper (Ns/m)

    # Wheel Masses
    m_wf: float            # Front wheel mass (kg)
    m_wr: float            # Rear wheel mass (kg)

    # Tyre stiffness
    k_tf: float            # Front tyre (N/m)
    k_tr: float            # Rear tyre (N/m)


@dataclass
class SimulationOptions:
    t_span: tuple          # t0 to t_end (s) start to end of simulation (length of road in our case)
    y_0: list              # Initial state, just the position of car at beginning - will be 0
    r_tolerance: float = 1e-7 # Tolerances required for solving later
    a_tolerance: float = 1e-9
    dense: bool = True     # Dense solution again required for solving later


# THIS NEXT SECTION OF THE CODE IS EXPLAINED IN THE "IMPORTING THE ROAD PROFILE" SHEET OF NOTES


# Road profile, the CSV has column "distance" (m) and "height" (m)
BaseInput = Callable[[float], Tuple[float, float, float, float]]

# Import CSV
df = pd.read_csv("bumpy_road_cords.csv")
x = df["distance"].values # Distance labelled as x
y = df["height"].values # Road height labelled as y

# Smooth road profile y(x) using UnivariateSpline import 
spline = UnivariateSpline(x, y, s = 0.4) # s = 0.4 is used to define how smooth you want the curve to be, higher value more smooth
dsdx = spline.derivative()   # dy/dx for road inputs later

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
                  [0.0, I]], dtype=float)

    # Damping matrix
    C = np.array([[c_f + c_r, a * c_f - b * c_r],
                  [a * c_f - b * c_r, a * a * c_f + b * b * c_r]], dtype=float)

    # Stiffness matrix
    K = np.array([[k_f + k_r, a * k_f - b * k_r],
                  [a * k_f - b * k_r, a * a * k_f + b * b * k_r]], dtype=float)

    return M, C, K

# Finding natural frequencies of body
def undamped_naturals(p: CarParams):
    M, _, K = build_matrices_mck(p)
    lam, _ = np.linalg.eig(np.linalg.solve(M, K))  # Eigenvalues of M^(-1)K
    wn = np.sqrt(np.clip(lam, 0.0, None))          # Finding frequency (rad/s)
    return np.sort(wn / (2 * np.pi))               # Rad/s to Hz


