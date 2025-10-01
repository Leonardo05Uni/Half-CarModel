# CMM-Group9
computational modeling and design group 9 project

## Car Spring-Mass-Damper Modeling Framework

This framework provides a computational model for simulating car motion using spring-mass-damper systems to track wheel behavior over varying road profiles.

### Features

- **Spring-Mass-Damper Model**: Each wheel (front and rear) is modeled as a spring-mass-damper system
- **CSV Road Profile Input**: Load road profiles from CSV files with distance and height data
- **Two-Wheel Tracking**: Tracks front and rear wheels (assuming left/right wheels on same axle behave identically)
- **Time-Domain Simulation**: Uses numerical integration (scipy.odeint) for accurate dynamics
- **Visualization**: Generates plots showing wheel positions, velocities, and suspension displacement

### Installation

1. Install required dependencies:
```bash
pip install -r requirements.txt
```

### Usage

#### Basic Example

```python
from car_model import WheelModel, CarModel

# Create wheel models with mass, spring stiffness, and damping
front_wheel = WheelModel(mass=400, spring_stiffness=20000, damping_coefficient=2000)
rear_wheel = WheelModel(mass=500, spring_stiffness=25000, damping_coefficient=2500)

# Create car model with wheelbase and velocity
car = CarModel(front_wheel, rear_wheel, wheelbase=2.5, velocity=10.0)

# Load road profile from CSV
car.load_road_profile('example_road_profile.csv')

# Run simulation
results = car.simulate(duration=5.0, dt=0.01)
```

#### Run Example Simulation

Basic example:
```bash
python example_simulation.py
```

This will:
1. Load the example road profile from `example_road_profile.csv`
2. Simulate the car motion for 5 seconds
3. Generate plots showing wheel positions, velocities, and suspension displacement
4. Save results to `simulation_results.png`

Advanced examples with various scenarios:
```bash
python advanced_examples.py
```

This demonstrates:
1. **Sinusoidal road profile**: Smooth periodic bumps
2. **Speed bump scenario**: Response to a single large bump
3. **Random rough road**: Realistic road surface with multiple frequencies
4. **Damping comparison**: Effect of different damping coefficients on ride quality

### Road Profile CSV Format

The CSV file should contain two columns:
- **Column 1**: Distance along road (meters)
- **Column 2**: Road height (meters)

Example:
```csv
distance,height
0.0,0.0
5.0,0.0
10.0,0.05
15.0,0.1
20.0,0.05
```

### Model Description

The framework uses a quarter-car model for each wheel pair (front and rear). Each wheel follows the dynamics:

```
m * z̈ = -k * (z - z_road) - c * (ż - ż_road)
```

Where:
- `m`: Mass supported by the wheel (kg)
- `k`: Spring stiffness (N/m)
- `c`: Damping coefficient (N·s/m)
- `z`: Vertical position of the mass (m)
- `z_road`: Vertical position of the road surface (m)

### Customization

#### Adjusting Car Parameters

Modify the parameters in `example_simulation.py` or create your own script:

```python
# Adjust masses (kg)
front_mass = 400
rear_mass = 500

# Adjust spring stiffness (N/m)
front_spring = 20000
rear_spring = 25000

# Adjust damping (N·s/m)
front_damping = 2000
rear_damping = 2500

# Adjust car geometry and motion
wheelbase = 2.5  # meters
velocity = 10.0  # m/s
```

#### Creating Custom Road Profiles

You can either:
1. Create a CSV file with distance and height columns
2. Use the `set_road_profile()` method directly:

```python
import numpy as np

# Create sinusoidal road profile
distance = np.linspace(0, 100, 1000)
height = 0.1 * np.sin(2 * np.pi * distance / 10)

car.set_road_profile(distance, height)
```

### Files

- `car_model.py`: Core framework with WheelModel and CarModel classes
- `example_simulation.py`: Basic example script demonstrating the framework
- `advanced_examples.py`: Advanced examples with various road scenarios
- `test_framework.py`: Test suite for validating framework functionality
- `example_road_profile.csv`: Sample road profile data
- `requirements.txt`: Python package dependencies

### Requirements

- Python 3.6+
- numpy >= 1.21.0
- matplotlib >= 3.4.0
- scipy >= 1.7.0
