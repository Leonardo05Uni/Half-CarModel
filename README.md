# CMM-Group9

### Overview

The performance and comfort of a vehicle are strongly influenced by its suspension system, which governs how the car responds to uneven road surfaces, so
this project aims to simulate a car traveling over a bumpy road to investigate how different suspension parameters affect ride comfort and stability.  

The simulation provides a simplified but informative representation of vehicle dynamics which is focusing on vertical (heave) and rotational (pitch) motion caused by road irregularities.  

Through this, our project aims to:
- Understand the dynamic interaction between road surface, suspension, and vehicle body.  
- Quantify how spring stiffness and damping coefficients influence vertical and angular acceleration.  
- Identify parameter combinations that minimize the Root Mean Square acceleration.

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

### Project File Overview

**car_model_2.py**
  – Core 2-DOF body dynamics model.  
  - Defines mass, damping, and stiffness matrices  
  - Implements equations of motion via 'rhs_car()'  
  - Integrates system with 'run_simulation()' using 'scipy.solve_ivp'

**road_profile.py**
  – Random and parameterized road generation.  
  - Generates bumps and potholes via Gaussian or polynomial surfaces  
  - Outputs to CSV ('bumpy_road_cords.csv') for external visualization  

**speed_bump.py** 
  – Simple deterministic bump input generator.  
  - Calculates height function 'h(x)' and its time derivatives 'h_dot(x, speed)' 
  - Used to test suspension response at specific vehicle speeds  

**bumpy_road_cords.csv**
  – Exported dataset of (distance [m], height [m])  
  - Used for validation and 2D road surface plotting  

**Model drawing – car v2.png**
  – Schematic of the 2-DOF suspension layout  
  - Annotates front/rear wheel positions and reference points


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




