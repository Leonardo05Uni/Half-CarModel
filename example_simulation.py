"""
Example Car Simulation
=======================

This script demonstrates how to use the car modeling framework to simulate
a car traveling over a road profile with spring-mass-damper systems for the wheels.
"""

import numpy as np
import matplotlib.pyplot as plt
from car_model import WheelModel, CarModel


def main():
    """Run example car simulation."""
    
    print("=" * 60)
    print("Car Spring-Mass-Damper Simulation")
    print("=" * 60)
    
    # Define car parameters
    front_mass = 400  # kg (mass supported by front wheels)
    rear_mass = 500   # kg (mass supported by rear wheels)
    
    # Spring stiffness (N/m)
    front_spring = 20000
    rear_spring = 25000
    
    # Damping coefficient (N·s/m)
    front_damping = 2000
    rear_damping = 2500
    
    # Create wheel models
    front_wheel = WheelModel(front_mass, front_spring, front_damping)
    rear_wheel = WheelModel(rear_mass, rear_spring, rear_damping)
    
    # Car parameters
    wheelbase = 2.5  # meters (distance between front and rear axles)
    velocity = 10.0  # m/s (car forward velocity)
    
    # Create car model
    car = CarModel(front_wheel, rear_wheel, wheelbase, velocity)
    
    print("\nCar Parameters:")
    print(f"  Front wheel mass: {front_mass} kg")
    print(f"  Rear wheel mass: {rear_mass} kg")
    print(f"  Front spring stiffness: {front_spring} N/m")
    print(f"  Rear spring stiffness: {rear_spring} N/m")
    print(f"  Front damping: {front_damping} N·s/m")
    print(f"  Rear damping: {rear_damping} N·s/m")
    print(f"  Wheelbase: {wheelbase} m")
    print(f"  Velocity: {velocity} m/s")
    
    # Load road profile from CSV
    print("\nLoading road profile from 'example_road_profile.csv'...")
    car.load_road_profile('example_road_profile.csv')
    print(f"  Road profile loaded: {len(car.road_distance)} points")
    
    # Run simulation
    simulation_duration = 5.0  # seconds
    print(f"\nRunning simulation for {simulation_duration} seconds...")
    results = car.simulate(simulation_duration, dt=0.01)
    
    print(f"  Simulation complete: {len(results['time'])} time steps")
    
    # Plot results
    print("\nGenerating plots...")
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # Plot 1: Wheel positions vs road profile
    axes[0].plot(results['distance'], results['front_position'], 
                 label='Front Wheel', linewidth=2)
    axes[0].plot(results['distance'], results['rear_position'], 
                 label='Rear Wheel', linewidth=2)
    axes[0].plot(results['distance'], results['front_road'], 
                 label='Road (Front)', linestyle='--', alpha=0.7)
    axes[0].plot(results['distance'], results['rear_road'], 
                 label='Road (Rear)', linestyle='--', alpha=0.7)
    axes[0].set_xlabel('Distance (m)')
    axes[0].set_ylabel('Height (m)')
    axes[0].set_title('Wheel Positions vs Road Profile')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: Wheel velocities
    axes[1].plot(results['time'], results['front_velocity'], 
                 label='Front Wheel', linewidth=2)
    axes[1].plot(results['time'], results['rear_velocity'], 
                 label='Rear Wheel', linewidth=2)
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylabel('Velocity (m/s)')
    axes[1].set_title('Wheel Vertical Velocities')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # Plot 3: Displacement from road
    front_displacement = results['front_position'] - results['front_road']
    rear_displacement = results['rear_position'] - results['rear_road']
    
    axes[2].plot(results['time'], front_displacement, 
                 label='Front Wheel', linewidth=2)
    axes[2].plot(results['time'], rear_displacement, 
                 label='Rear Wheel', linewidth=2)
    axes[2].set_xlabel('Time (s)')
    axes[2].set_ylabel('Displacement (m)')
    axes[2].set_title('Suspension Displacement (Wheel Position - Road Height)')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_file = 'simulation_results.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"  Results saved to '{output_file}'")
    
    # Show plot (optional, may not work in all environments)
    # plt.show()
    
    print("\n" + "=" * 60)
    print("Simulation complete!")
    print("=" * 60)


if __name__ == '__main__':
    main()
