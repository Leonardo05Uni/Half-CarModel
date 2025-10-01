"""
Advanced Examples for Car Modeling Framework
============================================

This script demonstrates various use cases and customizations of the framework:
1. Sinusoidal road profile (smooth bumps)
2. Speed bump scenario
3. Random rough road
4. Comparison of different damping coefficients
"""

import numpy as np
import matplotlib.pyplot as plt
from car_model import WheelModel, CarModel


def example_sinusoidal_road():
    """Example with smooth sinusoidal road profile."""
    print("\n" + "=" * 60)
    print("Example 1: Sinusoidal Road Profile")
    print("=" * 60)
    
    # Create car model
    front_wheel = WheelModel(400, 20000, 2000)
    rear_wheel = WheelModel(500, 25000, 2500)
    car = CarModel(front_wheel, rear_wheel, wheelbase=2.5, velocity=15.0)
    
    # Create sinusoidal road profile
    distance = np.linspace(0, 100, 500)
    height = 0.08 * np.sin(2 * np.pi * distance / 8)  # 8m wavelength, 0.08m amplitude
    car.set_road_profile(distance, height)
    
    # Simulate
    print("Running simulation...")
    results = car.simulate(duration=6.0, dt=0.01)
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(results['distance'], results['front_position'], 
            label='Front Wheel', linewidth=2)
    ax.plot(results['distance'], results['rear_position'], 
            label='Rear Wheel', linewidth=2)
    ax.plot(results['distance'], results['front_road'], 
            label='Road Profile', linestyle='--', linewidth=2, alpha=0.7)
    ax.set_xlabel('Distance (m)', fontsize=12)
    ax.set_ylabel('Height (m)', fontsize=12)
    ax.set_title('Car Response to Sinusoidal Road Profile', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('example_sinusoidal.png', dpi=150)
    print("Saved plot: example_sinusoidal.png")


def example_speed_bump():
    """Example with speed bump."""
    print("\n" + "=" * 60)
    print("Example 2: Speed Bump")
    print("=" * 60)
    
    # Create car model
    front_wheel = WheelModel(400, 20000, 2000)
    rear_wheel = WheelModel(500, 25000, 2500)
    car = CarModel(front_wheel, rear_wheel, wheelbase=2.5, velocity=8.0)
    
    # Create speed bump profile
    distance = np.linspace(0, 50, 500)
    height = np.zeros_like(distance)
    
    # Speed bump at 20m (Gaussian-shaped)
    bump_center = 20.0
    bump_width = 1.5
    bump_height = 0.12
    height = bump_height * np.exp(-((distance - bump_center) / bump_width) ** 2)
    
    car.set_road_profile(distance, height)
    
    # Simulate
    print("Running simulation...")
    results = car.simulate(duration=6.0, dt=0.01)
    
    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    # Position plot
    axes[0].plot(results['distance'], results['front_position'], 
                 label='Front Wheel', linewidth=2)
    axes[0].plot(results['distance'], results['rear_position'], 
                 label='Rear Wheel', linewidth=2)
    axes[0].plot(results['distance'], results['front_road'], 
                 label='Speed Bump', linestyle='--', linewidth=2, alpha=0.7, color='red')
    axes[0].set_ylabel('Height (m)', fontsize=12)
    axes[0].set_title('Speed Bump Response - Position', fontsize=14, fontweight='bold')
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3)
    
    # Velocity plot
    axes[1].plot(results['time'], results['front_velocity'], 
                 label='Front Wheel', linewidth=2)
    axes[1].plot(results['time'], results['rear_velocity'], 
                 label='Rear Wheel', linewidth=2)
    axes[1].set_xlabel('Time (s)', fontsize=12)
    axes[1].set_ylabel('Vertical Velocity (m/s)', fontsize=12)
    axes[1].set_title('Speed Bump Response - Velocity', fontsize=14, fontweight='bold')
    axes[1].legend(fontsize=11)
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('example_speed_bump.png', dpi=150)
    print("Saved plot: example_speed_bump.png")


def example_rough_road():
    """Example with random rough road."""
    print("\n" + "=" * 60)
    print("Example 3: Random Rough Road")
    print("=" * 60)
    
    # Create car model
    front_wheel = WheelModel(400, 20000, 2000)
    rear_wheel = WheelModel(500, 25000, 2500)
    car = CarModel(front_wheel, rear_wheel, wheelbase=2.5, velocity=12.0)
    
    # Create rough road profile with random noise
    distance = np.linspace(0, 80, 800)
    np.random.seed(42)  # For reproducibility
    
    # Combine multiple frequency components for realistic roughness
    height = (0.02 * np.sin(2 * np.pi * distance / 3) +
              0.01 * np.sin(2 * np.pi * distance / 1.5) +
              0.005 * np.random.randn(len(distance)))
    
    car.set_road_profile(distance, height)
    
    # Simulate
    print("Running simulation...")
    results = car.simulate(duration=6.5, dt=0.01)
    
    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    # Position plot
    axes[0].plot(results['distance'], results['front_position'], 
                 label='Front Wheel', linewidth=1.5, alpha=0.8)
    axes[0].plot(results['distance'], results['rear_position'], 
                 label='Rear Wheel', linewidth=1.5, alpha=0.8)
    axes[0].plot(results['distance'], results['front_road'], 
                 label='Road Profile', linestyle='--', linewidth=1, alpha=0.5, color='gray')
    axes[0].set_ylabel('Height (m)', fontsize=12)
    axes[0].set_title('Rough Road Response - Position', fontsize=14, fontweight='bold')
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3)
    
    # Acceleration (approximate from velocity)
    dt = results['time'][1] - results['time'][0]
    front_accel = np.gradient(results['front_velocity'], dt)
    rear_accel = np.gradient(results['rear_velocity'], dt)
    
    axes[1].plot(results['time'], front_accel, 
                 label='Front Wheel', linewidth=1.5, alpha=0.8)
    axes[1].plot(results['time'], rear_accel, 
                 label='Rear Wheel', linewidth=1.5, alpha=0.8)
    axes[1].set_xlabel('Time (s)', fontsize=12)
    axes[1].set_ylabel('Vertical Acceleration (m/s²)', fontsize=12)
    axes[1].set_title('Rough Road Response - Acceleration', fontsize=14, fontweight='bold')
    axes[1].legend(fontsize=11)
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('example_rough_road.png', dpi=150)
    print("Saved plot: example_rough_road.png")


def example_damping_comparison():
    """Compare different damping coefficients."""
    print("\n" + "=" * 60)
    print("Example 4: Damping Coefficient Comparison")
    print("=" * 60)
    
    # Create road profile (single bump)
    distance = np.linspace(0, 60, 500)
    height = 0.1 * np.exp(-((distance - 30) / 2) ** 2)
    
    # Test different damping values
    damping_values = [500, 2000, 5000]
    colors = ['red', 'blue', 'green']
    labels = ['Low Damping (500)', 'Medium Damping (2000)', 'High Damping (5000)']
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    print("Testing different damping coefficients...")
    for damping, color, label in zip(damping_values, colors, labels):
        # Create car model with specific damping
        front_wheel = WheelModel(400, 20000, damping)
        rear_wheel = WheelModel(500, 25000, damping)
        car = CarModel(front_wheel, rear_wheel, wheelbase=2.5, velocity=10.0)
        car.set_road_profile(distance, height)
        
        # Simulate
        results = car.simulate(duration=6.0, dt=0.01)
        
        # Plot front wheel response
        axes[0].plot(results['time'], results['front_position'], 
                     label=f'{label}', linewidth=2, color=color, alpha=0.8)
        
        # Plot suspension displacement
        displacement = results['front_position'] - results['front_road']
        axes[1].plot(results['time'], displacement, 
                     label=f'{label}', linewidth=2, color=color, alpha=0.8)
    
    # Plot road profile reference
    axes[0].plot(results['distance'] / 10, results['front_road'], 
                 label='Road Profile', linestyle='--', linewidth=2, 
                 color='black', alpha=0.5)
    
    axes[0].set_ylabel('Height (m)', fontsize=12)
    axes[0].set_title('Effect of Damping on Wheel Position', fontsize=14, fontweight='bold')
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_xlabel('Time (s)', fontsize=12)
    axes[1].set_ylabel('Suspension Displacement (m)', fontsize=12)
    axes[1].set_title('Effect of Damping on Suspension Displacement', fontsize=14, fontweight='bold')
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.3)
    axes[1].axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('example_damping_comparison.png', dpi=150)
    print("Saved plot: example_damping_comparison.png")


def main():
    """Run all examples."""
    print("=" * 60)
    print("Advanced Examples - Car Modeling Framework")
    print("=" * 60)
    
    # Run examples
    example_sinusoidal_road()
    example_speed_bump()
    example_rough_road()
    example_damping_comparison()
    
    print("\n" + "=" * 60)
    print("All examples completed!")
    print("=" * 60)
    print("\nGenerated files:")
    print("  - example_sinusoidal.png")
    print("  - example_speed_bump.png")
    print("  - example_rough_road.png")
    print("  - example_damping_comparison.png")


if __name__ == '__main__':
    main()
