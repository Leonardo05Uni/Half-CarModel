# CMM-Group9
computer modeling and design group 9 project
This project models a 2-DOF vehicle body system responding to road surface inputs such as bumps and potholes.  
It simulates the heave (vertical) and pitch (rotational) motion of the vehicle using linear spring–damper dynamics,  
with equations of motion integrated numerically in Python.

Legend: 
Suffix meaning:
_P - Position vector
_V - Velocity vector
_M - Mass
_theta - orientation / Angle
_phi - Damping coefficient
_k - stiffness coefficient
_L - Length
_a - length to front wheel from CG
_b - length to back wheel from CG
_dL - Change in length
_inertia - Inertia of sprung mass
_z - height (z value of car movement)
_dot - differentiated
_f - front of car
_r - rear of car
_dot_dot – second derivative (acceleration)  
_c – damping coefficient 
_base – baseline input (e.g., zero_base flat road)  
_series – time-dependent data series (e.g., speed_bump_series)  
_theta_dot – angular velocity (pitch rate)  
_theta_dot_dot – angular acceleration (pitch angular rate)  
_z_dot – vertical velocity  
_z_dot_dot – vertical acceleration  
_pass – passenger-related variable  
_x – longitudinal position along car body

Component / Positions Names (number reference those found in figure 1 on 'Model drawing - car v2.png')
1 - RWS - rear wheel spring
2 - RWD - rear wheel damper
3 - RWC - rear wheel centre
4 - RWP - rear wheel connection (point)
5 - FWS - Front wheel spring
6 - FWD - Front Wheel Damper
7 - FWC - Front wheel centre
8 - FWP - Front wheel connection (point)
9 - DP - Drivers position
10 - body - Sprung mass
11 - F - Force

Project File Overview

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



