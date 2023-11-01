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
control = smdController(Fmax=6, sigma=0.05, flag=False)
animation = smdAnimation(limits=2, multfigs=True)

# set values for controller
tr = 2
wn = 2.2/tr
damprat = 0.7
a = 2*damprat*wn
c = wn**2
control.kD = P.m*(a - P.b/P.m)
control.kP = P.m*(c - P.k/P.m)
control.kI = 0.5

# add subplots
z_plot = animation.fig.add_subplot(3, 2, 2)
zdot_plot = animation.fig.add_subplot(3, 2, 4)
f_plot = animation.fig.add_subplot(3, 2, 6)

# empty lists for plotting
sim_times = []
zs = []
zdots = []
fs = []
z_targets = []

# set initial values
z = 0
u = 0
z_target = 1

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop

    u = control.updatePID(z_target, smd.h())
    z = smd.update(u)

    sim_times.append(t)
    zs.append(z[0])
    zdots.append(smd.state[1][0])
    fs.append(u[0])
    z_targets.append(z_target)

    animation.update(smd.state)

    z_plot.clear(); zdot_plot.clear(), f_plot.clear()

    z_plot.plot(sim_times, zs, label="state", color='c')
    z_plot.plot(sim_times, z_targets, label="target", color='m')
    z_plot.legend(loc="upper left")
    z_plot.set_ylabel("z (m)")
    z_plot.grid()
    
    zdot_plot.plot(sim_times, zdots, color='c')
    zdot_plot.set_ylabel("zdot (m/s)")
    zdot_plot.grid()

    f_plot.plot(sim_times, fs, color = 'c')
    f_plot.set_ylabel("Force (N)")
    f_plot.grid()

    plt.pause(0.01)
    t += P.Ts

    if keyboard.is_pressed('q'): break
