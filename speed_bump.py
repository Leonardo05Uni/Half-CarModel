import numpy as np

_H = 0.5 # Maximum height of speed bump
_B = 2 # Maximum length of speed bump

def speed_bump_height(x):
    if (x>=0) and (x<=_B): # Checks if position is within length of speed bump
        height = (_H/2)*(1 - np.cos(2*np.pi*x/_B)) # Formula for speed bump
    else:
        height = 0 # 0 if out of range 

    return height # output of the height of the speed bump at that point (x)
