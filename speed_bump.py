import numpy as np # Importing numpy
import matplotlib.pyplot as plt # Importing matplotlib for plotting graphs

bump_H = 0.5 # Maximum height of speed bump
bump_L = 2 # Maximum length of speed bump

margin=bump_L*0.75 # the margin on each side of the plot for good looking
x_center=bump_L/2 # centre of the speed bump
x=np.linspace(-margin,bump_L+margin,400) # range of x axis, 400 displays the amount of dots
h=np.where((x>=0)&(x<=bump_L), 0.5*bump_H*(1-np.cos(2*np.pi*x/bump_L)),0) # function of speed bump
plt.figure(figsize=(bump_H+5,bump_L+5)) # figure size that displayed
plt.plot(x,h,color='black',linewidth=1) # color and linewidth of the graph
plt.xlim(x_center-margin, x_center+margin) # centre the plot, range of x axis
plt.ylim(0,2.5) # range of y axis
plt.tight_layout()
plt.axis("off") # hide the axis, make it clearer
plt.show()

# Height function h(x)
def h(x): # Height of the bump at position x
    x = np.array(x)
    h_val = np.where((x >= 0) & (x <= bump_L),0.5 * bump_H * (1 - np.cos(2 * np.pi * x / bump_L)),0)
    return h_val
# Time derivative h_dot(x)
"""
Notes:
    dh/dx = (pi*H/L)*sin(2*pi*x/L)
    to get dh/dt, multiply by dx/dt (dh/dx * dx/dt).
    dx/dt is just the speed.
    so y_f_dot and y_r_dot will just be their value of dh/dx * speed.
"""
def h_dot(x, speed): #Vertical velocity (dh/dt) at position x for given horizontal speed.
    x = np.array(x)
    dh_dx = np.where((x >= 0) & (x <= bump_L),(np.pi * bump_H / bump_L) * np.sin(2 * np.pi * x / bump_L),0)
    return dh_dx * speed

# Main function: series output
def speed_bump_series(L_a, L_b, speed_mph, dt=0.01):
    """
    Return continuous series of (y_f, y_r, y_f_dot, y_r_dot) as arrays.
    - L_a, L_b: distance from CG to front/rear axle (m)
    - speed_mph: vehicle speed in mph
    - dt: timestep (s)
    """
    speed = speed_mph / 2.237  # mph → m/s
    car_L = L_a + L_b

    # simulate from front entering to rear leaving bump
    total_time = (bump_L + car_L) / speed
    t = np.arange(0, total_time, dt)

    # wheel positions over time
    x_f = speed * t
    x_r = speed * t - car_L

    # heights and derivatives
    y_f = h(x_f)
    y_r = h(x_r)
    y_f_dot = h_dot(x_f, speed)
    y_r_dot = h_dot(x_r, speed)

    return t, y_f, y_r, y_f_dot, y_r_dot

np.set_printoptions(threshold=np.inf)

t, y_f, y_r, y_f_dot, y_r_dot = speed_bump_series(1.2, 1.3, 15)

print("Time array:\n", t)
print("Front wheel height:\n", y_f)
print("Rear wheel height:\n", y_r)
print("Front wheel velocity:\n", y_f_dot)
print("Rear wheel velocity:\n", y_r_dot)
