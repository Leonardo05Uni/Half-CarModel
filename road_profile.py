import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def generate_bumpy_road(csv_name='road.csv',length=100, resolution=0.01, incline=0,
                        bump_height=0.1, bump_width=1.2,
                        pothole_depth=0.1, pothole_width=2,
                        num_bumps=2, num_potholes=2, imperfection=0.01):
    """

    Parameters:
        length (float): Total length of the road (m)
        resolution (float): Distance step (m)
        incline (float): Slope of the road (rise/run)
        bump_height (float): Height of speed bumps (m) 
            - Max speed bump height is 0.1m (https://streetsolutionsuk.co.uk/blogs/news/the-ultimate-guide-to-speed-bump-regulations-uk)
        bump_width (float): Width of speed bumps (m)
            - Minimum length of speed bumps is 0.9cm (https://streetsolutionsuk.co.uk/blogs/news/the-ultimate-guide-to-speed-bump-regulations-uk)
        pothole_depth (float): Depth of potholes (m) 
            - Larger than 40mm is considered a pothole (https://www.darlington.gov.uk/transport-roads-and-parking/highways/potholes)
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
        
        x_local = (x - pos) / (bump_width / 2)
    
        #Step height is 0.025m (https://streetsolutionsuk.co.uk/blogs/news/the-ultimate-guide-to-speed-bump-regulations-uk)
        step_height = 0.025

        n = 2  #Modelled as an x^2 graph
        core_bump = (1 - x_local**n) * (bump_height - step_height)
        core_bump[np.abs(x_local) > 1] = 0  

        y_bump = np.copy(core_bump)


        left_edge = pos - bump_width / 2
        right_edge = pos + bump_width / 2

        
        y_bump[x >= left_edge] += step_height

        
        y_bump[x >= right_edge] -= step_height

        y += y_bump

    # Flat potholes
    for pos in pothole_positions:
        half_w = pothole_width / 2
        transition = half_w * 0.10  # 10% of width used to model smooth dip

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
    bumpy_y = np.copy(y)

    # Define a mask for flat regions (regions not affected by bumps or potholes)
    flat_region = np.ones_like(y, dtype=bool)

    # Keeping bumps perfectly smooth
    for pos in bump_positions:
        left_edge = pos - bump_width / 2
        right_edge = pos + bump_width / 2
        flat_region[(x >= left_edge) & (x <= right_edge)] = False

    # Keeping potholes perfectly smooth
    for pos in pothole_positions:
        left_edge = pos - pothole_width / 2
        right_edge = pos + pothole_width / 2
        flat_region[(x >= left_edge) & (x <= right_edge)] = False

    # Apply imperfections only to flat areas
    bumpy_y[flat_region] += np.random.uniform(-imperfection, imperfection, size=np.sum(flat_region))

    # Create DataFrame with proper structure
    df = pd.DataFrame({'distance': x, 'height': bumpy_y})
    df.to_csv(csv_name, index=False)

    return df, x, y, bump_positions, pothole_positions


def plot_road(df, x_base, y_base):
    plt.figure(figsize=(12, 5))
    

    
    # Plot bumpy road from dataframe
    plt.plot(df['distance'], df['height'], color='black', linewidth=2, label='Road surface')
    plt.title("2D Bumpy Road Profile", fontsize=14)
    plt.xlabel("Distance (m)")
    plt.ylabel("Height (m)")
    plt.grid(True, linestyle='--', alpha=0.5)#

    # Plot baseline (without imperfections)
    plt.plot(x_base, y_base, color='red', linewidth=1.5, alpha=0.7, label='Baseline', linestyle='--')

    plt.ylim(-0.5, 0.5)

    plt.legend()
    plt.show()


# Pre generating roads
df, x, y, bumps, potholes = generate_bumpy_road(csv_name='bumpy_road_cords.csv',length=100, resolution=0.01, incline=0,
                        bump_height=0.1, bump_width=1.2,
                        pothole_depth=0.1, pothole_width=2,
                        num_bumps=2, num_potholes=2, imperfection=0.01)
plot_road(df, x, y)

df, x, y, bumps, potholes = generate_bumpy_road(csv_name='speedbump.csv',length=10, resolution=0.1, incline=0,
                        bump_height=0.1, bump_width=1.2,
                        pothole_depth=0.1, pothole_width=2,
                        num_bumps=1, num_potholes=0, imperfection=0)
plot_road(df, x, y)

df, x, y, bumps, potholes = generate_bumpy_road(csv_name='motorway.csv',length=500, resolution=0.01, incline=0,
                        bump_height=0.1, bump_width=1.2,
                        pothole_depth=0.1, pothole_width=2,
                        num_bumps=0, num_potholes=0, imperfection=0.01)
                        
plot_road(df, x, y)

