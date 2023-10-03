import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
import matplotlib.pyplot as plt
import library.smdParam as P
from library.signalGenerator import signalGenerator
from library.smdAnimation import smdAnimation
from library.smdDynamics import smdDynamics
from library.smdController import smdController
import keyboard

# instantiate SMD, controller, and animation
smd = smdDynamics(alpha=0.0)
control = smdController()
animation = smdAnimation()

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
z_target = 0

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    if t <= 2:
        z_target = 0
    elif t <= 20:
        z_target = 1
    elif t <= 30:
        z_target = 1.5
    elif t <= 40:
        z_target = 0.5
    else:
        z_target = 0

    u = control.update(z_target, smd.state)
    z = smd.update(u)

    sim_times.append(t)
    zs.append(z)
    fs.append(u)
    z_targets.append(z_target)

    animation.update(smd.state)

    z_plot.clear(); f_plot.clear()
    z_plot.plot(sim_times, zs, label="state", color='c')
    z_plot.plot(sim_times, z_targets, label="target", color='pink')
    z_plot.legend(loc="lower right")
    z_plot.set_ylabel("z (m)")
    z_plot.grid()

    f_plot.plot(sim_times, fs, color = 'c')
    f_plot.set_ylabel("Force (N)")
    f_plot.grid()

    plt.pause(0.001)
    t += P.Ts

    if keyboard.is_pressed('q'): break

# # Keeps the program from closing until the user presses a button.
# print('Press key to close')
# plt.waitforbuttonpress()
# plt.close()
