import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
import matplotlib.pyplot as plt
import library.vtolParam2 as P
import numpy as np
from library.vtolAnimation import vtolAnimation
from library.vtolDynamics2 import vtolDynamics
from library.vtolController2 import vtolController2
from library.signalGenerator import signalGenerator
import keyboard

# instantiate VTOL, controller, and animation
VTOL = vtolDynamics(alpha=0.0)
control = vtolController2()
animation = vtolAnimation(limits=10, multfigs=True)

# signal input
z_val = signalGenerator(amplitude=2.5, frequency=0.08)

# Set parameters to tune controller
# h
tr_h = 7.5
wn_h = 2.2/tr_h
damprat = 0.707
b = 2*damprat*wn_h
c = wn_h**2
control.kDh = b*(P.mc + 2*P.mr)
control.kPh = c*(P.mc + 2*P.mr)
# theta
tr_th = 0.8
wn_th = 2.2/tr_th
b = 2*damprat*wn_th
c = wn_th**2
control.kDth = b*(P.Jc + 2*P.mr*P.d**2)
control.kPth = c*(P.Jc + 2*P.mr*P.d**2)
# z
tr_z = 4.5
wn_z = 2.2/tr_z
b = 2*damprat*wn_z
c = wn_z**2
control.kDz = (b - P.mu/(P.mc + 2*P.mr))/-P.g
control.kPz = c/-P.g

# add subplots
z_plot = animation.fig.add_subplot(4, 2, 2)
h_plot = animation.fig.add_subplot(4, 2, 4)
theta_plot = animation.fig.add_subplot(4, 2, 6)
f_plot = animation.fig.add_subplot(4, 2, 8)

# create lists for plotting
zs = [0.]
thetas = [0.]
hs = [0.1]
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
            z_target = 3 + z_val.square(t)
    
        fr, fl = control.update(z_target, h_target, VTOL.state)
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
    # z_plot.set_ylim(-1,7)
    z_plot.legend(loc="upper left")
    z_plot.set_ylabel("Z (m)")
    z_plot.grid()

    h_plot.plot(sim_times, hs, label="state", color='c')
    h_plot.plot(sim_times, h_targets, label="target", color='m')
    # h_plot.set_ylim(0,10)
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
