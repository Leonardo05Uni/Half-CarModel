import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random

def generate_bumpy_road(length=100, resolution=0.01, incline=0,
                        bump_height=0.1, bump_width=0.2,
                        pothole_depth=0.1, pothole_width=2,
                        num_bumps=1, num_potholes=4, imperfection=0.01):
    """

    Parameters:
        length (float): Total length of the road (m)
        resolution (float): Distance step (m)
        incline (float): Slope of the road (rise/run)
        bump_height (float): Height of speed bumps (m)
        bump_width (float): Width of speed bumps (m)
        pothole_depth (float): Depth of potholes (m)
        pothole_width (float): Width of potholes (m)
        num_bumps (int): Number of speed bumps
        num_potholes (int): Number of potholes
        imperfection (float): Random road grain error (m)
    """

    # Generate base road profile
    x = np.arange(0, length, resolution)
    y = incline * x 

    # Randomly place bumps and potholes
    np.random.seed(42)  # for repeatability
    bump_positions = np.random.uniform(10, length - 10, num_bumps)
    pothole_positions = np.random.uniform(10, length - 10, num_potholes)

    # Add speed bumps (modeled as Gaussian curves)
    for pos in bump_positions:
        y += bump_height * np.exp(-0.5 * ((x - pos) / (bump_width / 2)) ** 2)

    # Flat potholes
    for pos in pothole_positions:
        half_w = pothole_width / 2
        transition = half_w * 0.10  # 10% of width used to model smoothish dip

        # Flat section at the bottom of pothole
        inside = (x > pos - half_w + transition) & (x < pos + half_w - transition)
        y[inside] -= pothole_depth

        # Smooth transition edges as linear lines
        left_edge = (x >= pos - half_w) & (x < pos - half_w + transition)
        right_edge = (x > pos + half_w - transition) & (x <= pos + half_w)

        # Linear slope from 0 → pothole_depth
        y[left_edge] -= pothole_depth * ((x[left_edge] - (pos - half_w)) / transition)
        y[right_edge] -= pothole_depth * (1 - (x[right_edge] - (pos + half_w - transition)) / transition)

    # Add random imperfections
    bumpy_y = []
    for pos in range(len(y)):
        bumpy_y.append(random.uniform(y[pos] - imperfection, y[pos] + imperfection))

    # Create DataFrame with proper structure
    df = pd.DataFrame({'distance': x, 'height': bumpy_y})
    df.to_csv('bumpy_road_cords.csv', index=False)

    return df, x, y, bump_positions, pothole_positions


def plot_road(df, x_base, y_base, bump_positions, pothole_positions):
    plt.figure(figsize=(12, 5))
    

    
    # Plot bumpy road from dataframe
    plt.plot(df['distance'], df['height'], color='black', linewidth=2, label='Road surface')
    plt.title("2D Bumpy Road Profile", fontsize=14)
    plt.xlabel("Distance (m)")
    plt.ylabel("Height (m)")
    plt.grid(True, linestyle='--', alpha=0.5)#

    # Plot baseline (without imperfections)
    plt.plot(x_base, y_base, color='red', linewidth=1.5, alpha=0.7, label='Baseline', linestyle='--')

    plt.legend()
    plt.show()


# Generate road and plot
df, x, y, bumps, potholes = generate_bumpy_road()
plot_road(df, x, y, bumps, potholes)