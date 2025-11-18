from car_model_5 import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#================Making the road profile from data==================
#in here you will need to call upon the road_profile module to generate the road profile and pass along parameters, then loop through all below code for each profile

p = CarParams(
    body_M = 1163 - (29 * 4), # body mass subtracting the wheel masses
    body_inertia = 3000.0,
    body_a = 0.996, # Distance from CG to front axle
    body_b = 1.494, # Distance from CG to rear axle

    FWS_k = 30100, # Front wheel spring stiffness
    FWD_c = 2000.0, # Placeholder Values
    RWS_k = 32000, # Rear wheel spring stiffness
    RWD_c = 2000.0, # Placeholder Values

    m_wf = 58, # Front wheel mass
    m_wr = 58, # Rear wheel mass
    k_tf = 200000, # Front tire stiffness
    k_tr = 200000, # Rear tire stiffness
)

#=================== run simulation starting here===================

main(p)
plt.show()