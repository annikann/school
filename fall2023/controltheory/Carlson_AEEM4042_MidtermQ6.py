import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
import matplotlib.pyplot as plt
import numpy as np
import library.bbParam as P
from library.bbForceAnimation import bbForceAnimation
from library.bbDynamics import bbDynamics
from library.bbController import bbController
from library.signalGenerator import signalGenerator
import keyboard

# instantiate VTOL, controller, and animation
ballbeam = bbDynamics(alpha=0.0)
control = bbController()
reference = signalGenerator(amplitude=0.15, frequency=0.01)
disturbance = signalGenerator(amplitude=1.0)
animation = bbForceAnimation(limits=0.7, multfigs=True)

# add subplots
z_plot = animation.fig.add_subplot(2, 2, 2)
f_plot = animation.fig.add_subplot(2, 2, 4)

# empty lists for plotting
sim_times = []
zs = []
fs = []
zrs = []

# set initial values
# z = 0
# u = 0
# zr = 0

print('Press Q to end simulation')
t = P.t_start  # time starts at t_start
# y = ballbeam.h()
while t < P.t_end:
    t_next_plot = t + P.t_plot

    while t < t_next_plot:

        zr = 0.25 + reference.square(t)
        d = disturbance.step(t)
        F = control.update(zr, ballbeam.state)
        y = ballbeam.update(F + d)  # Propagate the dynamics

        sim_times.append(t)
        zs.append(y[0][0])
        fs.append(F)
        zrs.append(zr)

        t += P.Ts

    animation.update(ballbeam.state)


    z_plot.clear(); f_plot.clear()
    # z_plot.plot(sim_times, zs, label="state", color='c')
    z_plot.plot(sim_times, zrs, label="target", color='m')
    z_plot.legend(loc="upper left")
    z_plot.set_ylabel("z (m)")
    z_plot.grid()

    f_plot.plot(sim_times, fs, color = 'c')
    f_plot.set_ylabel("Total Force (N)")
    f_plot.grid()

    plt.pause(0.0001)

    if keyboard.is_pressed("q"): break
