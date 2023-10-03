import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
import matplotlib.pyplot as plt
import library.vtolParam as P
from library.vtolAnimation import vtolAnimation
from library.vtolDynamics import vtolDynamics
from library.vtolController import vtolController
import keyboard

# instantiate VTOL, controller, and animation
VTOL = vtolDynamics(alpha=0.0)
control = vtolController()
animation = vtolAnimation(limits=10, multfigs=True)

# add subplots
z_plot = animation.fig.add_subplot(2, 2, 2)
f_plot = animation.fig.add_subplot(2, 2, 4)
z_plot.set_ylabel("z (m)")
f_plot.set_ylabel("Force (N)")

# empty lists for plotting
sim_times = []
zs = []
fs = []
z_targets = []

# set initial values
z = 0
u = 0
z_target = 0.25

t = P.t_start  # time starts at t_start
while t < P.t_end:
    if t <= 2:
        z_target = 0.25
    elif t <= 30:
        z_target = 5
    elif t <= 50:
        z_target = 8
    elif t <= 65:
        z_target = 0.5
    else:
        z_target = 0.25

    F = control.update(z_target, VTOL.state)
    y = VTOL.update(F/2., F/2.)  # Propagate the dynamics

    sim_times.append(t)
    zs.append(y[1][0])
    fs.append(F)
    z_targets.append(z_target)

    animation.update(VTOL.state)

    z_plot.clear(); f_plot.clear()
    z_plot.plot(sim_times, zs, label="state", color='c')
    z_plot.plot(sim_times, z_targets, label="target", color='m')
    z_plot.legend(loc="upper left")
    z_plot.set_ylabel("z (m)")
    z_plot.grid()

    f_plot.plot(sim_times, fs, color = 'c')
    f_plot.set_ylabel("Total Force (N)")
    f_plot.grid()

    plt.pause(0.01)

    t += P.Ts
    if keyboard.is_pressed("q"): break
