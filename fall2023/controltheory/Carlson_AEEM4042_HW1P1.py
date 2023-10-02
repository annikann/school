import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
import matplotlib.pyplot as plt
import library.smdParam as P
from library.signalGenerator import signalGenerator
from library.smdAnimation import smdAnimation
from library.dataPlotter import dataPlotter
from library.smdDynamics import smdDynamics

# instantiate pendulum, controller, and reference classes
smd = smdDynamics(alpha=0.0)
reference = signalGenerator(amplitude=0.5, frequency=0.02)
force = signalGenerator(amplitude=1, frequency=1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = smdAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics at rate Ts
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        u = force.sin(t)
        y = smd.update(u)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts
    # update animation and data plots at rate t_plot
    animation.update(smd.state)
    dataPlot.update(t, r, smd.state, u)
    plt.pause(0.0001)  # allows time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
