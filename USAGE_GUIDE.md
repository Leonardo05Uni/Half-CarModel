# Car Spring-Mass-Damper Framework - Usage Guide

## Overview

This framework provides a comprehensive solution for modeling car suspension dynamics using spring-mass-damper systems. It tracks the vertical motion of front and rear wheels as a car travels over a road profile.

## Quick Start

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Run Basic Example

```bash
python example_simulation.py
```

This will load a sample road profile and generate visualizations.

### 3. Run Advanced Examples

```bash
python advanced_examples.py
```

This demonstrates various scenarios including sinusoidal roads, speed bumps, rough terrain, and damping comparisons.

## Basic Usage

### Creating a Car Model

```python
from car_model import WheelModel, CarModel

# Define front wheel parameters
front_wheel = WheelModel(
    mass=400,                    # kg - mass supported by front wheels
    spring_stiffness=20000,      # N/m - spring constant
    damping_coefficient=2000     # N·s/m - damping coefficient
)

# Define rear wheel parameters
rear_wheel = WheelModel(
    mass=500,                    # kg
    spring_stiffness=25000,      # N/m
    damping_coefficient=2500     # N·s/m
)

# Create car model
car = CarModel(
    front_wheel=front_wheel,
    rear_wheel=rear_wheel,
    wheelbase=2.5,              # meters - distance between axles
    velocity=10.0               # m/s - forward velocity
)
```

### Loading Road Profile from CSV

```python
# Load from CSV file
car.load_road_profile('my_road_profile.csv')
```

CSV format:
```csv
distance,height
0.0,0.0
5.0,0.02
10.0,0.05
```

### Creating Road Profile Programmatically

```python
import numpy as np

# Example: Sinusoidal road
distance = np.linspace(0, 100, 500)
height = 0.1 * np.sin(2 * np.pi * distance / 10)
car.set_road_profile(distance, height)
```

### Running Simulation

```python
# Run simulation
results = car.simulate(
    duration=5.0,    # seconds
    dt=0.01          # time step (seconds)
)

# Access results
time = results['time']
front_position = results['front_position']
rear_position = results['rear_position']
front_velocity = results['front_velocity']
rear_velocity = results['rear_velocity']
front_road = results['front_road']
rear_road = results['rear_road']
distance = results['distance']
```

## Advanced Usage Examples

### Example 1: Sports Car Configuration

```python
# Stiff suspension for better handling
front_wheel = WheelModel(300, 30000, 3000)
rear_wheel = WheelModel(350, 35000, 3500)
sports_car = CarModel(front_wheel, rear_wheel, 2.3, 20.0)
```

### Example 2: Comfort-Oriented Setup

```python
# Softer suspension for comfort
front_wheel = WheelModel(500, 15000, 1500)
rear_wheel = WheelModel(600, 18000, 1800)
luxury_car = CarModel(front_wheel, rear_wheel, 2.8, 12.0)
```

### Example 3: Speed Bump Analysis

```python
import numpy as np

# Create speed bump
distance = np.linspace(0, 50, 500)
bump_center = 25.0
bump_width = 1.5
bump_height = 0.12
height = bump_height * np.exp(-((distance - bump_center) / bump_width) ** 2)

car.set_road_profile(distance, height)
results = car.simulate(duration=4.0)
```

### Example 4: Rough Road Simulation

```python
import numpy as np

# Random rough road
distance = np.linspace(0, 100, 1000)
np.random.seed(42)
height = (0.02 * np.sin(2 * np.pi * distance / 3) +
          0.01 * np.sin(2 * np.pi * distance / 1.5) +
          0.005 * np.random.randn(len(distance)))

car.set_road_profile(distance, height)
results = car.simulate(duration=8.0)
```

## Visualization

### Basic Plotting

```python
import matplotlib.pyplot as plt

# Plot wheel positions
plt.figure(figsize=(12, 6))
plt.plot(results['distance'], results['front_position'], label='Front Wheel')
plt.plot(results['distance'], results['rear_position'], label='Rear Wheel')
plt.plot(results['distance'], results['front_road'], label='Road', linestyle='--')
plt.xlabel('Distance (m)')
plt.ylabel('Height (m)')
plt.legend()
plt.grid(True)
plt.savefig('output.png')
```

### Suspension Displacement

```python
# Calculate suspension compression/extension
front_displacement = results['front_position'] - results['front_road']
rear_displacement = results['rear_position'] - results['rear_road']

plt.figure(figsize=(12, 6))
plt.plot(results['time'], front_displacement, label='Front')
plt.plot(results['time'], rear_displacement, label='Rear')
plt.xlabel('Time (s)')
plt.ylabel('Suspension Displacement (m)')
plt.legend()
plt.grid(True)
plt.savefig('suspension.png')
```

## Parameter Guidelines

### Mass (kg)
- Small car: 300-500 per wheel
- Mid-size car: 400-600 per wheel
- Large car/SUV: 500-800 per wheel

### Spring Stiffness (N/m)
- Soft (comfort): 15,000-20,000
- Medium: 20,000-25,000
- Stiff (sport): 25,000-35,000

### Damping Coefficient (N·s/m)
- Light damping: 500-1,500
- Medium damping: 1,500-2,500
- Heavy damping: 2,500-4,000

### Wheelbase (m)
- Compact car: 2.3-2.5
- Mid-size car: 2.5-2.8
- Large car: 2.8-3.2

### Velocity (m/s)
- Slow: 5-10 m/s (18-36 km/h)
- Medium: 10-20 m/s (36-72 km/h)
- Fast: 20-30 m/s (72-108 km/h)

## Physics Background

The framework uses a quarter-car model where each wheel follows the equation:

```
m * z̈ = -k * (z - z_road) - c * (ż - ż_road)
```

Where:
- `m`: Mass supported by the wheel
- `k`: Spring stiffness
- `c`: Damping coefficient
- `z`: Vertical position of the mass
- `z_road`: Road surface height

This is integrated numerically using `scipy.integrate.odeint` for accurate time-domain simulation.

## Testing

Run the test suite to verify the framework:

```bash
python test_framework.py
```

This tests:
- Model creation and initialization
- Road profile loading (CSV and arrays)
- Simulation execution
- Interpolation accuracy

## Troubleshooting

### Simulation is unstable
- Reduce time step (dt): try 0.005 instead of 0.01
- Check parameter values are reasonable
- Ensure road profile is smooth (no discontinuities)

### Results look wrong
- Verify units are correct (meters, seconds, kg)
- Check road profile scale (height should be < 1m typically)
- Ensure velocity matches road profile extent

### CSV file won't load
- Check format: two columns (distance, height)
- Remove any empty rows
- Ensure no invalid characters

## Further Customization

The framework can be extended to include:
- Tire dynamics
- Body pitch and roll
- Multiple road profiles (left/right wheels)
- Non-linear spring characteristics
- Active suspension systems

See the source code in `car_model.py` for implementation details.
