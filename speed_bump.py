import numpy as np # Importing numpy
import matplotlib.pyplot as plt # Importing matplotlib for plotting graphs

bump_H = 0.5 # Maximum height of speed bump
bump_L = 2 # Maximum length of speed bump

margin=bump_L*0.75
x_center=bump_L/2
x=np.linspace(-margin,bump_L+margin,400)
h=np.where((x>=0)&(x<=bump_L), 0.5*bump_H*(1-np.cos(2*np.pi*x/bump_L)),0)
plt.figure(figsize=(bump_H+5,bump_L+5))
plt.plot(x,h,color='black',linewidth=1)
plt.xlim(x_center-margin, x_center+margin)
plt.ylim(0,5)
plt.tight_layout()
plt.axis("off")
plt.show()

def speed_bump_height(x):
    if (x>=0) and (x<=bump_L): # Checks if position is within length of speed bump
        height = (bump_H/2)*(1 - np.cos(2*np.pi*x/bump_L)) # Formula for speed bump
    else:
        height = 0 # 0 if out of range 

    return height # output of the height of the speed bump at that point (x)
