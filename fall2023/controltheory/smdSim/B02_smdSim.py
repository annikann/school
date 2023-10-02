import sys 
sys.path.append('/Users/annikacarlson/Documents/ControlTheory/smdSim')
import matplotlib.pyplot as plt
import numpy as np
import parameters.smdParam as P
from tools.signalGenerator import signalGenerator
from viewer.smdAnimation import smdAnimation
from viewer.dataPlotter import dataPlotter
from dynamics.smdDynamics import smdDynamics
import keyboard as key

# instantiate reference input classes
SMD = smdDynamics(alpha=0.0)
reference = signalGenerator(amplitude=0.5, frequency=0.1)
force = P.f0    # Force input

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = smdAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propogate dynamics
    t_next_plot = t + P.t_plot
    # set variables
    while t < t_next_plot:
        r = reference.square(t)
        u = force           # Constant force input 
        y = SMD.update(u)  
        t = t + P.t_plot     # advance time by t_plot 
    # update animations
    # state = np.array([[z], [0.0], [0.0]])
    animation.update(SMD.state)
    dataPlot.update(t, r, SMD.state, u)

    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
