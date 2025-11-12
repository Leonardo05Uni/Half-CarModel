from car_model_5 import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#================Making the road profile from data==================
#in here you will need to call upon the road_profile module to generate the road profile and pass along parameters, then loop through all below code for each profile


p = CarParams(
    body_M = 1163 - (29 * 4), # body mass subtracting the wheel masses
    body_inertia = 3000.0,
    body_a = 0.996,
    body_b = 1.494,

    FWS_k = 30100, # Front wheel spring stiffness
    FWD_c = 2000.0, # Front wheel damping coefficient
    RWS_k = 32000, # Rear wheel spring stiffness
    RWD_c = 2000.0, # Rear wheel damping coefficient

    m_wf = 58, # Front wheel mass
    m_wr = 58, # Rear wheel mass
    k_tf = 200000, # Front tire stiffness
    k_tr = 200000, # Rear tire stiffness
)

#=================== run simulation starting here===================
# Car travels along the measured road at a constant velocity of 8 m/s
road_base = make_road_base(p, spline, dsdx, v = 8.0, x0 = 0.0)

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
c_f_opt, c_r_opt, se_cf, se_cr = damping_monte_carlo_error(
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


#=======================and ending the simulation here=====================
#this will need to be turned into one large function that can be looped back through for each road profile
# PLOTS

plot_heave(ts, z)
plot_pitch(ts, theta)
plot_passenger_accel(ts, a_pass, seat_x)
plt.show()
