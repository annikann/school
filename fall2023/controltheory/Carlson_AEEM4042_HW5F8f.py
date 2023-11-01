import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
import matplotlib.pyplot as plt
import library.vtolParam as P
from library.vtolAnimation import vtolAnimation
from library.vtolDynamics import vtolDynamics
from library.vtolController2 import vtolController2
import keyboard

# instantiate VTOL, controller, and animation
VTOL = vtolDynamics(alpha=0.0)
control = vtolController2()
animation = vtolAnimation(limits=10, multfigs=True)

# Set parameters to tune controller
# h
tr_h = 2.61
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
tr_z = 2.9
wn_z = 2.2/tr_z
b = 2*damprat*wn_z
c = wn_z**2
control.kDz = (b - P.u/(P.mc + 2*P.mr))/-P.g
control.kPz = c/-P.g

# add subplots
z_plot = animation.fig.add_subplot(2, 2, 2)
f_plot = animation.fig.add_subplot(2, 2, 4)

# empty lists for plotting
sim_times = []
zs = []
fs = []
z_targets = []

# set initial values
z = 3
u = 0
z_target = 0.0

t = P.t_start  # time starts at t_start
while t < P.t_end:
    if t <= 2:
        z_target = 0.0
    elif t <= 25:
        z_target = 5
    elif t <= 50:
        z_target = 7
    elif t <= 65:
        z_target = 2
    elif t <= 80:
        z_target = 5
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
