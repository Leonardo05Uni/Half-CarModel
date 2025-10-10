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

def speed_bump_height(x):
    if (x>=0) and (x<=bump_L): # Checks if position is within length of speed bump
        height = (bump_H/2)*(1 - np.cos(2*np.pi*x/bump_L)) # Formula for speed bump
    else:
        height = 0 # 0 if out of range 

    return height # output of the height of the speed bump at that point (x)
