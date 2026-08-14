# CMM-Group9
## Group Members:  

Nigel Cheung, Boming Xiao, Floris Hijink, Leonardo Maffei Mercalli, Dimitri Rao  

## Overview

This code takes several minutes to run because each Monte Carlo trial repeats
the damping optimisation. The sample count can be configured in
`MasterFile.py`; 20 samples are used by default.

The suspension system has a major influence on how a car reacts to bumpy or
uneven roads, affecting both vehicle stability and passenger comfort.

This project simulates a car travelling over different types of roads to investigate how variations in damping coefficients affect the dynamics of the vehicle.

The simulation provides a physically meaningful representation of vehicle dynamics, which is focused on riding comfort through suspension tuning.  

**Reference vehicle:**  
Car: 2018 Ford Fiesta ST Line 5dr, 1.0L.  
https://media.ford.com/content/dam/fordmedia/Europe/documents/productReleases/Fiesta/FordFiesta2017_FiestaDrive_TechSpecs_EU.pdf 


### Project Objectives

Our main goals are：  
- Study how the road surface, tyres, suspension, and vehicle body interact with each other.  
- Observe how changes in damping affect both vertical and angular accelerations.  
- Find a combination of parameters that gives the smoothest ride, using RMS acceleration as a measure of comfort.  

## Individual Contribution

Leonardo Maffei Mercalli developed the mathematical background and authored
`MasterFile.py` and `car_model_5.py`, including the coupled equations of motion,
ODE integration, damping optimisation, Monte Carlo parameter-sensitivity
analysis, modal calculations, passenger-acceleration evaluation and plotting.
The road-profile module and wider project documentation were completed
collaboratively within the five-person team.

### System Component / Positions Names
(number reference those found in figure 1 on `Model drawing - car v2.png`)  

The vehicle is modelled as a half-car system consisting of:
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
| 10，11 | KTR / KTF | Rear / Front Tyre Spring Stiffness |

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
- imported parameters from `car_model_5.py`
- imported road layout from "road_profile.py"
- set up the initial conditions
- imported time integration through `run_simulation()`
- performed bisection-based damping optimisation and Monte Carlo parameter-sensitivity analysis using `damping_monte_carlo_error()`  
- computed RMS accelerations for comparison with published ride-comfort ranges
- generated plots including `Body Heave vs Time`, `Body Pitch vs Time`, `Passenger Vertical Accel (x={seat_x} m from CG)` and `Unsprung (wheel) accelerations`.
- can be used for every road profile and car type by changing the parameters.

### Model drawing – car v2.png  
Diagram showing the simplified 4-DOF suspension layout.  
- Front and rear wheel points are marked, along with the basic reference positions.

### car_model_5.py  
More car model information added, includes:
- wheel masses (`m_wf`, `m_wr`), tyre stiffness (`k_tf`, `k_tr`)
- defined `CarParams` and  `SimulationOptions`
- Formulates the system differential equations in `rhs_car()` and defines the equations of motion
- Performs numerical integration using `scipy.integrate.solve_ivp`
- Includes modal analysis, RMS evaluation, and passenger acceleration output functions

### road_profile.py  
- generates road profiles such as bumps and potholes using mathematical functions
- supports composite and spline-based surface construction
- allows combination of potholes, surface roughness, and gradients for realistic road shapes
- exports height data to a corresponding CSV
- being imported in `MasterFile.py` to create realistic road inputs for the vehicle model  

### Road Profiles Folder
- These are screenshots of the pre-generated profile

### Plots Folder
- Screenshots of outputted graphs from each road
- Aso including the optimised response of the vehicle

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

**2. Python Libraries required:**  
   - `numpy`, `scipy`, `matplotlib`, `pandas`

**3. What if changing the road profile and car type?**
   - All parameters can be changed manually in `MasterFile.py`, allowing different vehicle specifications and road profiles to be simulated.  
   - For instance, by modifying the vehicle properties (mass, stiffness, damping, etc.) to match a specific car model, the program could automatically produce the corresponding
     simulated dynamic response for that vehicle.


## Numerical Methods Implementation  

1. **Root finding**  
   Location -- `car_model_5.py`, function: `J(c)`  

   J(c) defines the Root Mean Square Passenger Acceleration as a function of damping coefficient.
   The Root Finding Method is aimed at finding out the root of $\frac{dJ}{dC}=0$.
   Bisection Method is used by setting up $a$ and $b$ as two guessing roots when $f(a)f(b)<0$, then set $m=\frac{a+b}{2}$, when $f(a)f(m)<0$, $b=m$, $f(b)f(m)<0$, $a=m$ to minimize       the space between the intervals to find out the roots.
   
2. **Ordinary Differential Equations (ODEs)**   
   Location -- `car_model_5.py`, function: `rhs_car()`
   
   Defines the governing motion equations for the car model including heave, pitch and wheel vertical motion. etc.  
   Integrated over time by using `scipy.integrate.solve_ivp` with RK45 Method.  

3. **Interpolation**  
   Location -- `road_profile.py`, function: `generate_bumpy_road()`  

   It uses interpolation techniques (`scipy.interpolate.interp1d`) to make the height points of discrete road smoother.  
   It generates continuous road surfaces including bumps, potholes and random roughness.  

4. **Monte Carlo parameter-sensitivity analysis**
   Location -- `car_model_5.py`, function: `damping_monte_carlo_error()`

   Each trial independently perturbs the nominal sprung mass and front/rear
   spring stiffnesses by ±5%, then repeats the damping optimisation. A fixed
   random seed makes results reproducible. The output reports the mean,
   standard deviation and empirical 95% interval of the optimal front and rear
   damping coefficients. The number of trials is configurable in
   `MasterFile.py`.
   
   

## RMS Acceleration and Ride Comfort Analysis

The **Root Mean Square (RMS)** acceleration quantifies the average vibration intensity over time:

$$
a_{RMS} = \sqrt{\frac{1}{T} \int_0^T [a(t)]^2 \, dt}
$$

### Indicative ISO 2631-1 Ride Comfort Ranges

The ranges below provide useful context for the RMS acceleration results.
However, this model currently uses unweighted vertical acceleration; a formal
ISO 2631-1 assessment would additionally require the specified frequency
weighting and exposure-duration treatment.

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




