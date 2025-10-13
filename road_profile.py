import numpy as np
import matplotlib.pyplot as plt

def generate_bumpy_road(length=100, resolution=0.01, incline=0,
                        bump_height=0.1, bump_width=0.2,
                        pothole_depth=0.1, pothole_width=2,
                        num_bumps=1, num_potholes=4):
    """
    Generate a 2D bumpy road profile.

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
    """

    # Generate base road profile
    x = np.arange(0, length, resolution)
    y = incline * x  # Base incline line

    # Randomly place bumps and potholes
    np.random.seed(42)  # for repeatability
    bump_positions = np.random.uniform(10, length - 10, num_bumps)
    pothole_positions = np.random.uniform(10, length - 10, num_potholes)

    # Add speed bumps (modeled as Gaussian curves)
    for pos in bump_positions:
        y += bump_height * np.exp(-0.5 * ((x - pos) / (bump_width / 2)) ** 2)

    # flat potholes
    for pos in pothole_positions:
        half_w = pothole_width / 2
        transition = half_w * 0.25  # 25% of width used for smooth edges

        # Flat section (bottom of pothole)
        inside = (x > pos - half_w + transition) & (x < pos + half_w - transition)
        y[inside] -= pothole_depth

        # Smooth transition edges (linear ramps)
        left_edge = (x >= pos - half_w) & (x < pos - half_w + transition)
        right_edge = (x > pos + half_w - transition) & (x <= pos + half_w)

        # Linear slope from 0 → pothole_depth
        y[left_edge] -= pothole_depth * ((x[left_edge] - (pos - half_w)) / transition)
        y[right_edge] -= pothole_depth * (1 - (x[right_edge] - (pos + half_w - transition)) / transition)

    return x, y, bump_positions, pothole_positions


def plot_road(x, y, bump_positions, pothole_positions):
    """Plot the road profile."""
    plt.figure(figsize=(12, 5))
    plt.plot(x, y, color='black', linewidth=2)
    plt.title("2D Bumpy Road Profile", fontsize=14)
    plt.xlabel("Distance (m)")
    plt.ylabel("Height (m)")
    plt.grid(True, linestyle='--', alpha=0.5)

    # Mark bump and pothole positions
    plt.scatter(bump_positions, [y[np.argmin(abs(x - pos))] for pos in bump_positions],
                color='orange', label='Speed bumps', zorder=5)
    plt.scatter(pothole_positions, [y[np.argmin(abs(x - pos))] for pos in pothole_positions],
                color='blue', label='Potholes', zorder=5)
    
    plt.legend()
    plt.show()


x, y, bumps, potholes = generate_bumpy_road()

plot_road(x, y, bumps, potholes)
