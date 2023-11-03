import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
import matplotlib.pyplot as plt
import library.vtolParam3 as P
import numpy as np
from library.vtolAnimation import vtolAnimation
from library.vtolDynamics2 import vtolDynamics
from library.vtolController3 import control
import keyboard

# instantiate VTOL, controller, and animation
VTOL = vtolDynamics(alpha=0.0)
ctrl = control()
animation = vtolAnimation(limits=10, multfigs=True)

# Set parameters to tune controller
kPth = 0.3721
kDth = 0.19134
kPz = -0.00771
kDz = -0.03285
kPh = 0.1134
kDh = 0.5833
kIh = 0.
kIth = 0.5
kIz = 0.

# add subplots
z_plot = animation.fig.add_subplot(4, 2, 2)
h_plot = animation.fig.add_subplot(4, 2, 4)
theta_plot = animation.fig.add_subplot(4, 2, 6)
f_plot = animation.fig.add_subplot(4, 2, 8)

# create lists for plotting
zs = [P.z0]
thetas = [0.]
hs = [P.h0]
frs = [0]
fls = [0]
h_target = 5.
z_target = 5.
h_targets = [h_target]
z_targets = [z_target]

t = P.t_start  # time starts at t_start
sim_times = [t]
while t < P.t_end:

    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        if t > 1.0:
            h_target = 5.
            z_target = 5.
    
        fr, fl = ctrl.update(h_target, z_target, VTOL.h())
        y = VTOL.update(fr, fl)
        t += P.Ts

    sim_times.append(t)
    zs.append(y[0][0])
    hs.append(y[1][0])
    thetas.append(np.rad2deg(y[2][0]))
    frs.append(fr)
    fls.append(fl)
    h_targets.append(h_target)
    z_targets.append(z_target)

    animation.update(VTOL.state)

    z_plot.clear(); h_plot.clear(); theta_plot.clear(); f_plot.clear()

    z_plot.plot(sim_times, zs, label="state", color='c')
    z_plot.plot(sim_times, z_targets, label="target", color='m')
    z_plot.legend(loc="upper left")
    z_plot.set_ylabel("Z (m)")
    z_plot.grid()

    h_plot.plot(sim_times, hs, label="state", color='c')
    h_plot.plot(sim_times, h_targets, label="target", color='m')
    h_plot.legend(loc="upper left")
    h_plot.set_ylabel("Height (m)")
    h_plot.grid()

    theta_plot.plot(sim_times, thetas, color='c')
    theta_plot.set_ylabel("Theta (deg)")
    theta_plot.grid()

    f_plot.plot(sim_times, frs, label="right", color = 'c')
    f_plot.plot(sim_times, fls, label="left", color = 'm')
    f_plot.legend(loc="upper left")
    f_plot.set_ylabel("Forces (N)")
    f_plot.grid()

    plt.pause(0.01)

    if keyboard.is_pressed("q"): break
