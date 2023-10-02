import sys 
sys.path.append('/Users/annikacarlson/Documents/ControlTheory/vtolSim')
import matplotlib.pyplot as plt
import numpy as np
import parameters.vtolParam as P
from tools.signalGenerator import signalGenerator
from viewer.vtolAnimation import vtolAnimation
from viewer.dataPlotter import dataPlotter
from dynamics.vtolDynamics import vtolDynamics

# instantiate reference input classes
VTOL = vtolDynamics(alpha=0.0)

# Run animation
animation = vtolAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        if t < 2:
            fr = 1.5*9.81/2
            fl = 1.5*9.81/2
        elif t < 5:
            fr = (1.5*9.81/2)*1.01
            fl = (1.5*9.81/2)*1.01
        elif t < 5.25:
            fr = (1.5*9.81/2)*0.98
            fl = (1.5*9.81/2)*1.02
        elif t < 5.50:
            fr = (1.5*9.81/2)*1.02
            fl = (1.5*9.81/2)*0.98     
        elif t < 10:
            fr = (1.5*9.81/2)*1.02
            fl = (1.5*9.81/2)*1.02
        else: 
            fl = (1.5*9.81/2)*1.00
            fr = (1.5*9.81/2)*1.00

        y = VTOL.update(fr, fl)
        t += P.Ts

    animation.update(VTOL.state)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
