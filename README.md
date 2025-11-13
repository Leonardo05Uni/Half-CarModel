# CMM-Group9
## Group Members:  

Nigel Cheung, Boming Xiao, Floris Hijink, Leonardo Maffei Mercalli, Dimitri Rao  

## Overview
The suspension system has a huge impact on how a car reacts to bumpy or uneven roads — it pretty much decides how stable and comfortable the ride feels.  
This project simulates a car travelling over a rough road to investigate how variations in spring stiffness and damping coefficients affect the vertical (heave) and rotational (pitch) motion of the vehicle.  
The simulation provides a physically meaningful representation of vehicle dynamics, which is focused on rideing comfort and vibration controlling through suspension tuning.  

**The reallife example:**  
Car: 2018 Ford Fiesta ST Line 5dr, 1.0L.  
(https://media.ford.com/content/dam/fordmedia/Europe/documents/productReleases/Fiesta/FordFiesta2017_FiestaDrive_TechSpecs_EU.pdf)  


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

- Small-angle motion is assumed ($\sin\theta \approx \theta$)  
- Suspension follows a linear spring–damper model  
- Tyre flexibility represented through spring stiffness $k_{tf}$ and $k_{tr}$  
- Front and rear wheel masses included ($m_{wf}$, $m_{wr}$)  
- Car body has been treated as a rigid sprung mass with heave $z$ and pitch $\theta$  
- Road surface is defined by flat, bump, pothole, and random surfaces  
- Air resistance and rolling friction are ignored  

## Project File Overview

### MasterFile.py
Main control script
- imported parameters from "car_model_5.py"
- imported road layout from "road_profile.py" and "speed_bump.py"
- set up the initial conditions
- imported time integration through "run_simulation()"
- performed Monte Carlo–based damping optimization by (`damping_monte_carlo_error`)  
- computed RMS accelerations and ISO 2631-1 comfort levels
- generated plots including "Body Heave vs Time", "Body Pitch vs Time", "Passenger Vertical Accel (x={seat_x} m from CG)" and "Unsprung (wheel) accelerations".


### Model drawing – car v2.png  
Diagram showing the simplified 2-DOF suspension layout.  
- Front and rear wheel points are marked, along with the basic reference positions.

### car_model_2.py  
Main script for the 2-DOF car body model  
- This is the main script where the 4-DOF car model runs  
- It sets the mass, damping,tyre and stiffness values, calls the motion equations from rhs_car()  
- Runs the simulation through `run_simulation()` using SciPy’s `solve_ivp`  

### car_model_5.py  
More car model information added, includes:
- wheel masses (`m_wf`, `m_wr`), tyre stiffness (`k_tf`, `k_tr`)
- defined `CarParams` and  `SimulationOptions`
- Formulates the system differential equations in rhs_car() and defines the equations of motion
- Performs numerical integration using scipy.integrate.solve_ivp
- Includes modal analysis, RMS evaluation, and passenger acceleration output functions

### road_profile.py  
- generates road profiles such as bumps and potholes using mathematical functions
- supports composite and spline-based surface construction
- allows combination of potholes, surface roughness, and gradients for realistic road shapes
- exports height data to bumpy_road_cords.csv
- being imported in MasterFile.py to create realistic road inputs for the vehicle model

### speed_bump.py  
A quick test file for a single speed bump.  
- It calculates the bump height h(x) and time-based h_dot(x, speed)  
- to check how the suspension behaves at different car speeds

### bumpy_road_cords.csv    
Just a CSV file with two columns — distance and height (metres)  
- It’s mainly for checking or plotting the road surface
 
### Two images  
**image_2025-11-10_134556331.png / image_2025-11-10_134858877.png**  

Those two images displays the composition of unsprung mass and the tyre vertical stiffness range (150-300 kN/m)  

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
 
## How to Run the code? 
**Programming language: Python 3**  
**1. Open Masterfile.py and this model will automatically import all of the modules below:**
   - car_model_5.py
   - road_profile.py
   - speed_bump.py

**2. Python Libraries required:**  
   - `numpy`, `scipy`, `matplotlib`, `pandas`  
   **How to install?**  
   Via pip:  
   pip install numpy scipy matplotlib pandas

**3. What happens when the code is running?**  
   - First off, the script called `road_profile.py` would generate a road input —— `bumpy_road_cords.csv`
   - Secondly, it runs the dynamic model from `car_model_5.py` via `rhs_car()` by using `solve_ivp`
   - Then the RMS acceleration would be generated automatically
   - At last the plottings including `Body Heave vs Time`, `Body Pitch vs Time`, `Passenger Vertical Accel (x={seat_x} m from CG)` and `Unsprung (wheel) accelerations`  
     would be generated.
    

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


## Outputs

- **Heave displacement:** `z(t)`  
- **Pitch angle:** `θ(t)`  
- **Wheel hub motion:** `z_wf(t)`, `z_wr(t)`  
- **Passenger acceleration:** `a_pass(t)`  
- **RMS comfort results and plots** 






