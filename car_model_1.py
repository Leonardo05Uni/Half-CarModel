# -*- coding: utf-8 -*-

from dataclasses import dataclass
import numpy as np
from typing import Callable, Tuple

@dataclass
class CarParams:
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


def build_matrices_mck(p: CarParams): # p is the variable that references directly from Carparams class

    #Setting up variables to make easier typing withing function    
    m, I = p.body_M, p.body_inertia
    a, b = p.body_a, p.body_b
    k_f, k_r = p.FWS_k, p.RWS_k
    c_f, c_r = p.FWD_c, p.RWD_c
    
    #M Matrix
    M = np.array([[m, 0.0],
                  [0.0, I]], dtype = float)
    
    #C Matrix
    C = np.array([[c_f + c_r, a*c_f - b*c_r],
                  [a*c_f - b*c_r, a*a*c_f + b*b*c_r]], dtype = float)
    
    #K Matrix
    K = np.array([[k_f + k_r, a*k_f - b*k_r],
                  [a*k_f - b*k_r, a*a*k_f + b*b*k_r]], dtype = float)
    
    return M, C, K



#Ignore this, looked up a decent way to generate fake values for the road input. Flat road.
BaseInput = Callable[[float], Tuple[float, float, float, float]]

def zero_base(_: float) -> Tuple[float, float, float, float]:
    """Default: flat road, no motion (useful for free-decay tests)."""
    
    #The point of this function is that once the road profile has been made someone should be returning 4 y values.
    
    return 0.0, 0.0, 0.0, 0.0

#This function is literally just the maths in my sheets, lmk if its confusing
def rhs_car(t, x, p: CarParams, base: BaseInput):
    
    z, theta, z_dot, theta_dot = x
    y_f, y_r, y_f_dot, y_r_dot = base(t) #road input (flat for now)
    
    a, b = p.body_a, p.body_b
    k_f, k_r = p.FWS_k, p.RWS_k
    c_f, c_r = p.FWD_c, p.RWD_c
    
    #Relative deflections on wheels (change in lengths)
    dL_f = (z + a*theta) - y_f
    dL_r = (z - b*theta) - y_r
    dL_f_dot = (z_dot + a*theta_dot) - y_f_dot
    dL_r_dot = (z_dot - b*theta_dot) - y_r_dot
    
    #Forces on spring damper system (up positive)
    F_f = k_f*dL_f + c_f*dL_f_dot
    F_r = k_r*dL_r + c_r*dL_r_dot
    
    #Equations of motion
    z_dot_dot = -(F_f + F_r)/p.body_M
    theta_dot_dot = -(a*F_f - b*F_r)/p.body_inertia
    
    #Essentially returning velocity and acceleration
    return [z_dot, theta_dot, z_dot_dot, theta_dot_dot]
    
    
















