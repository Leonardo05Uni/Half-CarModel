"""
Car Model with Spring-Mass-Damper System
=========================================

This module implements a simplified car model using spring-mass-damper systems
for the wheels. The model tracks two wheel positions (front and rear), assuming
left and right wheels on the same axle behave identically.

The quarter-car model uses the following differential equations:
    m * z_ddot = -k * (z - z_road) - c * (z_dot - z_road_dot)

where:
    - m: mass supported by the wheel
    - k: spring stiffness
    - c: damping coefficient
    - z: vertical position of the mass
    - z_road: vertical position of the road surface
"""

import numpy as np
from scipy.integrate import odeint
import csv


class WheelModel:
    """
    Represents a single wheel with spring-mass-damper system.
    
    Attributes:
        mass (float): Mass supported by this wheel (kg)
        spring_stiffness (float): Spring constant (N/m)
        damping_coefficient (float): Damping coefficient (N·s/m)
        position (float): Current vertical position (m)
        velocity (float): Current vertical velocity (m/s)
    """
    
    def __init__(self, mass, spring_stiffness, damping_coefficient):
        """
        Initialize wheel model parameters.
        
        Args:
            mass (float): Mass supported by the wheel (kg)
            spring_stiffness (float): Spring constant (N/m)
            damping_coefficient (float): Damping coefficient (N·s/m)
        """
        self.mass = mass
        self.spring_stiffness = spring_stiffness
        self.damping_coefficient = damping_coefficient
        self.position = 0.0
        self.velocity = 0.0
    
    def dynamics(self, state, t, road_height, road_velocity):
        """
        Compute the derivative of the state for the spring-mass-damper system.
        
        Args:
            state (array): [position, velocity]
            t (float): Time
            road_height (float): Road surface height at this time
            road_velocity (float): Road surface velocity at this time
        
        Returns:
            array: [velocity, acceleration]
        """
        position, velocity = state
        
        # Spring-mass-damper equation
        acceleration = (
            -self.spring_stiffness * (position - road_height) 
            - self.damping_coefficient * (velocity - road_velocity)
        ) / self.mass
        
        return [velocity, acceleration]


class CarModel:
    """
    Car model with front and rear wheels using spring-mass-damper systems.
    
    This model assumes:
    - Front left and right wheels behave identically
    - Rear left and right wheels behave identically
    - The car moves forward at constant velocity
    - Road profile is provided as a function of distance or time
    
    Attributes:
        front_wheel (WheelModel): Front wheel model
        rear_wheel (WheelModel): Rear wheel model
        wheelbase (float): Distance between front and rear axles (m)
        velocity (float): Forward velocity of the car (m/s)
    """
    
    def __init__(self, front_wheel, rear_wheel, wheelbase, velocity):
        """
        Initialize car model.
        
        Args:
            front_wheel (WheelModel): Front wheel model
            rear_wheel (WheelModel): Rear wheel model
            wheelbase (float): Distance between front and rear axles (m)
            velocity (float): Forward velocity of the car (m/s)
        """
        self.front_wheel = front_wheel
        self.rear_wheel = rear_wheel
        self.wheelbase = wheelbase
        self.velocity = velocity
        self.road_profile = None
        self.road_distance = None
    
    def load_road_profile(self, csv_file):
        """
        Load road profile from a CSV file.
        
        The CSV file should have two columns:
        - Column 1: Distance along road (m)
        - Column 2: Road height (m)
        
        Args:
            csv_file (str): Path to CSV file containing road profile
        """
        distance = []
        height = []
        
        with open(csv_file, 'r') as f:
            reader = csv.reader(f)
            next(reader, None)  # Skip header if present
            for row in reader:
                if len(row) >= 2:
                    try:
                        distance.append(float(row[0]))
                        height.append(float(row[1]))
                    except ValueError:
                        continue  # Skip invalid rows
        
        self.road_distance = np.array(distance)
        self.road_profile = np.array(height)
    
    def set_road_profile(self, distance, height):
        """
        Set road profile directly from arrays.
        
        Args:
            distance (array): Distance along road (m)
            height (array): Road height at each distance (m)
        """
        self.road_distance = np.array(distance)
        self.road_profile = np.array(height)
    
    def get_road_height(self, distance):
        """
        Get road height at a specific distance using interpolation.
        
        Args:
            distance (float): Distance along road (m)
        
        Returns:
            float: Road height at the given distance (m)
        """
        if self.road_profile is None or self.road_distance is None:
            return 0.0
        
        return np.interp(distance, self.road_distance, self.road_profile)
    
    def simulate(self, duration, dt=0.01):
        """
        Simulate the car motion over the road profile.
        
        Args:
            duration (float): Simulation duration (seconds)
            dt (float): Time step for simulation (seconds)
        
        Returns:
            dict: Dictionary containing simulation results with keys:
                - 'time': Time array
                - 'front_position': Front wheel vertical position
                - 'front_velocity': Front wheel vertical velocity
                - 'rear_position': Rear wheel vertical position
                - 'rear_velocity': Rear wheel vertical velocity
                - 'front_road': Road height at front wheel
                - 'rear_road': Road height at rear wheel
                - 'distance': Distance traveled
        """
        time = np.arange(0, duration, dt)
        n_steps = len(time)
        
        # Initialize state arrays
        front_position = np.zeros(n_steps)
        front_velocity = np.zeros(n_steps)
        rear_position = np.zeros(n_steps)
        rear_velocity = np.zeros(n_steps)
        front_road = np.zeros(n_steps)
        rear_road = np.zeros(n_steps)
        distance = np.zeros(n_steps)
        
        # Initial conditions
        front_state = [self.front_wheel.position, self.front_wheel.velocity]
        rear_state = [self.rear_wheel.position, self.rear_wheel.velocity]
        
        for i in range(n_steps):
            t = time[i]
            dist = self.velocity * t
            distance[i] = dist
            
            # Get road heights for front and rear wheels
            front_road_height = self.get_road_height(dist)
            rear_road_height = self.get_road_height(dist - self.wheelbase)
            
            front_road[i] = front_road_height
            rear_road[i] = rear_road_height
            
            # Store current state
            front_position[i] = front_state[0]
            front_velocity[i] = front_state[1]
            rear_position[i] = rear_state[0]
            rear_velocity[i] = rear_state[1]
            
            if i < n_steps - 1:
                # Estimate road velocity (finite difference)
                next_dist = self.velocity * time[i + 1]
                front_road_vel = (self.get_road_height(next_dist) - front_road_height) / dt
                rear_road_vel = (self.get_road_height(next_dist - self.wheelbase) - rear_road_height) / dt
                
                # Integrate one time step for front wheel
                t_span = [t, t + dt]
                front_sol = odeint(
                    self.front_wheel.dynamics,
                    front_state,
                    t_span,
                    args=(front_road_height, front_road_vel)
                )
                front_state = front_sol[-1]
                
                # Integrate one time step for rear wheel
                rear_sol = odeint(
                    self.rear_wheel.dynamics,
                    rear_state,
                    t_span,
                    args=(rear_road_height, rear_road_vel)
                )
                rear_state = rear_sol[-1]
        
        return {
            'time': time,
            'front_position': front_position,
            'front_velocity': front_velocity,
            'rear_position': rear_position,
            'rear_velocity': rear_velocity,
            'front_road': front_road,
            'rear_road': rear_road,
            'distance': distance
        }
