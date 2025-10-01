"""
Test script to verify the car modeling framework functionality.
"""

import numpy as np
from car_model import WheelModel, CarModel


def test_wheel_model():
    """Test basic wheel model creation."""
    print("Testing WheelModel creation...")
    wheel = WheelModel(mass=400, spring_stiffness=20000, damping_coefficient=2000)
    assert wheel.mass == 400
    assert wheel.spring_stiffness == 20000
    assert wheel.damping_coefficient == 2000
    assert wheel.position == 0.0
    assert wheel.velocity == 0.0
    print("  ✓ WheelModel creation successful")


def test_car_model():
    """Test car model creation."""
    print("Testing CarModel creation...")
    front_wheel = WheelModel(400, 20000, 2000)
    rear_wheel = WheelModel(500, 25000, 2500)
    car = CarModel(front_wheel, rear_wheel, wheelbase=2.5, velocity=10.0)
    
    assert car.wheelbase == 2.5
    assert car.velocity == 10.0
    assert car.front_wheel == front_wheel
    assert car.rear_wheel == rear_wheel
    print("  ✓ CarModel creation successful")


def test_road_profile_array():
    """Test setting road profile from arrays."""
    print("Testing road profile from arrays...")
    front_wheel = WheelModel(400, 20000, 2000)
    rear_wheel = WheelModel(500, 25000, 2500)
    car = CarModel(front_wheel, rear_wheel, wheelbase=2.5, velocity=10.0)
    
    distance = np.array([0, 10, 20, 30])
    height = np.array([0, 0.1, 0.05, 0])
    car.set_road_profile(distance, height)
    
    assert len(car.road_distance) == 4
    assert len(car.road_profile) == 4
    assert car.get_road_height(10) == 0.1
    print("  ✓ Road profile from arrays successful")


def test_road_profile_csv():
    """Test loading road profile from CSV."""
    print("Testing road profile from CSV...")
    front_wheel = WheelModel(400, 20000, 2000)
    rear_wheel = WheelModel(500, 25000, 2500)
    car = CarModel(front_wheel, rear_wheel, wheelbase=2.5, velocity=10.0)
    
    car.load_road_profile('example_road_profile.csv')
    
    assert car.road_distance is not None
    assert car.road_profile is not None
    assert len(car.road_distance) > 0
    print(f"  ✓ Road profile from CSV successful ({len(car.road_distance)} points)")


def test_simulation():
    """Test running a simulation."""
    print("Testing simulation...")
    front_wheel = WheelModel(400, 20000, 2000)
    rear_wheel = WheelModel(500, 25000, 2500)
    car = CarModel(front_wheel, rear_wheel, wheelbase=2.5, velocity=10.0)
    
    # Create simple sinusoidal road profile
    distance = np.linspace(0, 100, 200)
    height = 0.05 * np.sin(2 * np.pi * distance / 10)
    car.set_road_profile(distance, height)
    
    # Run short simulation
    results = car.simulate(duration=2.0, dt=0.01)
    
    assert 'time' in results
    assert 'front_position' in results
    assert 'rear_position' in results
    assert 'front_velocity' in results
    assert 'rear_velocity' in results
    assert 'front_road' in results
    assert 'rear_road' in results
    assert 'distance' in results
    
    assert len(results['time']) == len(results['front_position'])
    assert len(results['time']) == len(results['rear_position'])
    
    print(f"  ✓ Simulation successful ({len(results['time'])} time steps)")


def test_interpolation():
    """Test road height interpolation."""
    print("Testing road height interpolation...")
    front_wheel = WheelModel(400, 20000, 2000)
    rear_wheel = WheelModel(500, 25000, 2500)
    car = CarModel(front_wheel, rear_wheel, wheelbase=2.5, velocity=10.0)
    
    distance = np.array([0, 10, 20, 30])
    height = np.array([0, 0.1, 0.05, 0])
    car.set_road_profile(distance, height)
    
    # Test exact points
    assert abs(car.get_road_height(0) - 0) < 1e-6
    assert abs(car.get_road_height(10) - 0.1) < 1e-6
    assert abs(car.get_road_height(30) - 0) < 1e-6
    
    # Test interpolated point (midpoint between 10 and 20)
    mid_height = car.get_road_height(15)
    assert 0.05 < mid_height < 0.1  # Should be between the two values
    
    print("  ✓ Interpolation successful")


def main():
    """Run all tests."""
    print("=" * 60)
    print("Car Modeling Framework - Test Suite")
    print("=" * 60)
    print()
    
    try:
        test_wheel_model()
        test_car_model()
        test_road_profile_array()
        test_road_profile_csv()
        test_simulation()
        test_interpolation()
        
        print()
        print("=" * 60)
        print("All tests passed! ✓")
        print("=" * 60)
        return 0
    except Exception as e:
        print()
        print("=" * 60)
        print(f"Test failed: {e}")
        print("=" * 60)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    exit(main())
