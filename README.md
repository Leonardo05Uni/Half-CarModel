# CMM-Group9
## Group Members:  

Nigel Cheung, Boming Xiao, Floris Hijink, Leonardo Maffei Mercalli, Dimitri Rao  

## Overview
The suspension system has a huge impact on how a car reacts to bumpy or uneven roads — it pretty much decides how stable and comfortable the ride feels.  
This project simulates a car travelling over a rough road to investigate how variations in spring stiffness and damping coefficients affect the vertical (heave) and rotational (pitch) motion of the vehicle.  
The simulation provides a physically meaningful representation of vehicle dynamics, which is focused on rideing comfort and vibration controlling through suspension tuning.  

### Project Objectives  

Our main goals are：  
- Study how the road surface, tyres, suspension, and vehicle body interact with each other.  
- Observe how stiffness and damping changes affect both vertical and angular accelerations.  
- Find a combination of parameters that gives the smoothest ride, using RMS acceleration as a measure of comfort.  

### System Component / Positions Names
(number reference those found in figure 1 on 'Model drawing - car v2.png')  

The vehicle is modelled as a 4-DOF half-car system consisting of:
- **Sprung mass** which is the vehicle body  
- **Front and rear unsprung masses** about the wheels and axles   
- **Suspension springs/dampers** connecting the body to the wheels  
- **Tyre stiffness** connecting the wheels to the road surface  
 
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

### Assumptions and Limitations

- Small-angle motion assumed ($\sin\theta \approx \theta$)  
- Suspension follows a linear spring–damper model  
- Tyre flexibility represented through spring stiffness $k_{tf}$ and $k_{tr}$  
- Front and rear wheel masses included ($m_{wf}$, $m_{wr}$)  
- Car body treated as a rigid sprung mass with heave $z$ and pitch $\theta$  
- Road surface is defined by flat, bump, pothole, and random surfaces  
- Effects of air resistance and rolling friction are ignored  


## Mathematical Modelling
### Road Input – Speed Bump Profile
**Definition:**  
The road bump is represented as a smooth half-sine profile:  

$$
h(x) =
\begin{cases}
\frac{H}{2}\*[1 - \cos\left(\frac{2\pi x}{B}\right)], & 0 \le x \le B \\
0 & \text{otherwise}
\end{cases}
$$

where  
- \(H\) is the bump height,  
- \(B\) is the bump base length,  
- \(x\) is the longitudinal position along the road profile.

### Pothole profile: 
**Definition:**
Random small-scale perturbations are added to the baseline to emulate road imperfections and surface roughness.  

$$
P(x) =
\begin{cases}
0 & |x - p| \ge a \\
-d *\dfrac{a - |x - p|}{t} & a - t < |x - p| < a \\
-d & |x - p| \le a - t
\end{cases}
$$  

where  
- \(P(x)\) is the road vertical displacement (pothole depth function)   
- \(x\) is the longitudinal position along the road profile (m)  
- \(p\) is the center position of the pothole (m)   
- \(a\) is the half-width of the pothole (m), can also be explained as the distance from center to edge   
- \(t\) is the transition width (m) which defines the sloped region between road and pothole bottom   
- \(d\) is the maximum pothole depth (m), the positive value indicates downward deflection   

Interpretation:  
- When $|x - p| \ge a$: the surface is flat → no pothole $(P(x) = 0)$  
- When $a - t < |x - p| < a$: the surface slopes linearly downward  
- When $|x - p| \le a - t$: the surface stays at full depth $-d$  

### Vehicle Dynamic Equations
Let  
- `z` = vertical displacement of the car body’s centre of gravity 
- `θ` = pitch angle of the car body (positive in counterclockwise rotation)  
- `z_wf`, `z_wr` = vertical displacement of front and rear unsprung masses (wheel hubs)  
- `y_f`, `y_r` = front and rear road surface inputs  
- `a`, `b` = distances from the CG to the front and rear suspension connection points  
- `k_f`, `k_r` = suspension stiffness (front / rear)  
- `c_f`, `c_r` = suspension damping coefficients (front / rear)  
- `k_tf`, `k_tr` = tyre stiffness (front / rear)  
- `m_wf`, `m_wr` = unsprung masses (front / rear wheel assemblies)  
- `M` = sprung mass of the vehicle body  
- `I` = rotational inertia of the car body about the CG

**Kinematic Relationships**

The vertical displacement of the front and rear wheel hubs:

$$
z_f(t) = z(t) + a\theta(t)
$$

$$
z_r(t) = z(t) - b\theta(t)
$$

Differentiating gives the corresponding velocities:

$$
\dot{z_f}(t) = \dot{z}(t) + a\dot{\theta}(t)
$$

$$
\dot{z_r}(t) = \dot{z}(t) - b\dot{\theta}(t)
$$

**Suspension Deflections**  

$$
\delta_f = (z + a\theta) - z_{wf}
$$

$$
\delta_r = (z - b\theta) - z_{wr}
$$

**Suspension Forces**  

$$
F_f = k_f \cdot \delta_f + c_f \cdot \dot{\delta_f}
$$

$$
F_r = k_r \cdot \delta_r + c_r \cdot \dot{\delta_r}
$$

**Tyre deflections:**

$$
\eta_f = z_{wf} - y_f
$$

$$
\eta_r = z_{wr} - y_r
$$

**Tyre forces:**

$$
F_{t,f} = k_{tf}\eta_f
$$

$$
F_{t,r} = k_{tr}\eta_r
$$

**Equations of Motion**  

$$
M \ddot{z} = - (F_{s,f} + F_{s,r})
$$  

$$
I \ddot{\theta} = - (aF_{s,f} - bF_{s,r})
$$  

$$
m_{wf} \ddot{z}_{wf} = F_{s,f} - F_{t,f}
$$  

$$
m_{wr} \ddot{z}_{wr} = F_{s,r} - F_{t,r}
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
| `_wf`         | Front wheel (unsprung mass) displacement |
| `_wr`         | Rear wheel (unsprung mass) displacement |
| `_tf`         | Front tyre stiffness |
| `_tr`         | Rear tyre stiffness |
| `_eta_f`      | Tyre deflection at front (`z_wf - y_f`) |
| `_eta_r`      | Tyre deflection at rear (`z_wr - y_r`) |
| `_delta_f`    | Suspension deflection at front (`(z + aθ) - z_wf`) |
| `_delta_r`    | Suspension deflection at rear (`(z - bθ) - z_wr`) |
| `_F_s`        | Suspension force (spring + damper) |
| `_F_t`        | Tyre force |

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
- This is the main script where the 4-DOF car model runs.  
- It sets the mass, damping,tyre and stiffness values, calls the motion equations from rhs_car()  
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

## Outputs

- **Heave displacement:** `z(t)`  
- **Pitch angle:** `θ(t)`  
- **Wheel hub motion:** `z_wf(t)`, `z_wr(t)`  
- **Passenger acceleration:** `a_pass(t)`  
- **RMS comfort results and plots** 






