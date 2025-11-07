# CMM-Group9

## Overview
The suspension system has a huge impact on how a car reacts to bumpy or uneven roads — it pretty much decides how stable and comfortable the ride feels.  

For this project, we made a simple model to test how different suspension setups, like spring stiffness and damping, change the car’s motion when it goes over rough ground.  
The simulation mainly looks at the car’s up-and-down (heave) and tilting (pitch) motion. It’s not a super-realistic car model, but it does a good job showing the main physics behind ride comfort.  

Instead of trying to copy every real-world detail, we just focused on what actually matters — how the suspension affects the way the car body moves.  

Our main goals are：  
- Study how the road surface, suspension, and vehicle body interact with each other.  
- Observe how stiffness and damping changes affect both vertical and angular accelerations.  
- Find a combination of parameters that gives the smoothest ride, using RMS acceleration as a measure of comfort.  

### System Component / Positions Names
(number reference those found in figure 1 on 'Model drawing - car v2.png')

The vehicle is represented as a 2-DOF rigid body model supported by front and rear spring–damper systems.  
Each wheel follows the input road profile independently, transmitting forces to the car body, which moves vertically and rotates about its center of gravity (CG).

| ID | Symbol | Description |
|----|---------|-------------|
| 1,5 | RWS / FWS | Rear / Front Wheel Spring |
| 2,6 | RWD / FWD | Rear / Front Wheel Damper |
| 3,7 | RWC / FWC | Rear / Front Wheel Centre |
| 4,8 | RWP / FWP | Rear / Front Wheel Connection (link) |
| 9 | Body | Vehicle Sprung Mass (rigid body) |

The **car body** is modeled as a rigid bar of mass \( M \) and rotational inertia \( I \).  
Front and rear suspensions act as spring–damper pairs located at distances \( a \) and \( b \) from the CG.

## Mathematical Modelling
### Road Input – Speed Bump Profile
**Definition:**  
The road bump is represented as a smooth half-sine profile:  

$$
h(x) =
\begin{cases}
\frac{H}{2}\,[1 - \cos\left(\frac{2\pi x}{B}\right)], & 0 \le x \le B \\
0, & \text{otherwise}
\end{cases}
$$

where  
- \(H\) is the bump height,  
- \(B\) is the bump base length,  
- \(x\) is the longitudinal position along the road profile.

### Vehicle Dynamic Equations
Let `z` = vertical displacement, `θ` = pitch angle,
`y_f`, `y_r` = front / rear wheel road input.  
**Relative Deflections**  

$$
\delta_f = (z + a\theta) - y_f
$$

$$
\delta_r = (z - b\theta) - y_r
$$

**Suspension Forces**  

$$
F_f = k_f \cdot \delta_f + c_f \cdot \dot{\delta_f}
$$

$$
F_r = k_r \cdot \delta_r + c_r \cdot \dot{\delta_r}
$$

**Equations of Motion**  

$$
M \ddot{z} = - (F_f + F_r)
$$

$$
I \ddot{\theta} = - (aF_f - bF_r)
$$ 


**The ODEs are solved using `scipy.integrate.solve_ivp` (Runge–Kutta RK45) with adaptive time-stepping and dense output.** 
### Numerical Simulation
- Programming Language: Python 3  
- Libraries: NumPy · SciPy · Matplotlib · Pandas  
- Integration method: RK45 (`solve_ivp`)  
- Sampling frequency: adaptive (dense output = True)  
- Post-processing: RMS acceleration calculation and plot generation


### Suffix Meaning

| Suffix        | Meaning |
|---------------|---------|
| `_P`          | Position vector |
| `_V`          | Velocity vector |
| `_M`          | Mass |
| `_theta`      | Orientation / Angle |
| `_phi`        | Damping coefficient |
| `_k`          | Stiffness coefficient |
| `_L`          | Length |
| `_a`          | Length to front wheel from CG |
| `_b`          | Length to back wheel from CG |
| `_dL`         | Change in length |
| `_inertia`    | Inertia of sprung mass |
| `_z`          | Height (z value of car movement) |
| `_dot`        | Differentiated (time derivative) |
| `_dot_dot`    | Second derivative (acceleration) |
| `_f`          | Front of car |
| `_r`          | Rear of car |
| `_c`          | Damping coefficient |
| `_base`       | Baseline input (e.g., `zero_base` flat road) |
| `_series`     | Time-dependent data series (e.g., `speed_bump_series`) |
| `_theta_dot`  | Angular velocity (pitch rate) |
| `_theta_dot_dot` | Angular acceleration (pitch angular rate) |
| `_z_dot`      | Vertical velocity |
| `_z_dot_dot`  | Vertical acceleration |
| `_pass`       | Passenger-related variable |
| `_x`          | Longitudinal position along car body |

## RMS Acceleration and Ride Comfort Analysis

The **Root Mean Square (RMS)** acceleration quantifies the average vibration intensity over time:

$$
a_{RMS} = \sqrt{\frac{1}{T} \int_0^T [a(t)]^2 \, dt}
$$

### ISO 2631-1 Ride Comfort Classification

| **RMS Acceleration (m/s²)** | **Comfort Level (ISO 2631-1)**      |
|-----------------------------:|:------------------------------------|
| < 0.315                     | Not uncomfortable                   |
| 0.315 – 0.63                | A little uncomfortable              |
| 0.5 – 1.0                   | Fairly uncomfortable                |
| 0.8 – 1.6                   | Uncomfortable                       |
| 1.25 – 2.5                  | Very uncomfortable                  |
| > 2.5                       | Extremely uncomfortable             |

### Example Results

| **Metric**                     | **Value (m/s²)** | **Comfort Rating**      |
|--------------------------------|:----------------:|:------------------------|
| RMS heave acceleration (CG)    | 0.045            | Not uncomfortable       |
| RMS passenger acceleration     | 0.060            | Not uncomfortable       |


## Project File Overview

### car_model_2.py  
Main script for the 2-DOF car body model.  
- This is the main script where the 2-DOF car model runs.  
- It sets the mass, damping and stiffness values, calls the motion equations from rhs_car()  
- Runs the simulation through `run_simulation()` using SciPy’s `solve_ivp`.

### road_profile.py  
Used for generating rough or random road surfaces.  
- It makes small bumps or potholes using simple math functions  
- Saves the road data into `bumpy_road_cords.csv` so it can be reused.

### speed_bump.py  
A quick test file for a single speed bump.  
- It calculates the bump height h(x) and time-based h_dot(x, speed).  
- to check how the suspension behaves at different car speeds.

### bumpy_road_cords.csv    
Just a CSV file with two columns — distance and height (metres).  
- It’s mainly for checking or plotting the road surface.

### Model drawing – car v2.png  
Diagram showing the simplified 2-DOF suspension layout.  
- Front and rear wheel points are marked, along with the basic reference positions.  








