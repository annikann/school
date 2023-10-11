import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
import matplotlib.pyplot as plt
import numpy as np
import library.bbParam as P
from library.bbForceAnimation import bbAnimation
from library.bbDynamics import bbDynamics
import keyboard

# instantiate VTOL, controller, and animation
ballbeam = bbDynamics(alpha=0.0)
animation = bbAnimation(limits=0.7, multfigs=True)

# set initial values
state = np.array([[0.5],   # z
                  [0.],    # zdot
                  [0.],    # theta
                  [0.]])   # thetadot
z0 = [0][0]
zdot0 = [1][0]
theta0 = [2][0]
thetadot0 = [3][0]

# add subplots
z_plot = animation.fig.add_subplot(4, 2, 2)
zdot_plot = animation.fig.add_subplot(4, 2, 4)
theta_plot = animation.fig.add_subplot(4, 2, 6)
thetadot_plot = animation.fig.add_subplot(4, 2, 8)
F_plot = animation.fig.add_subplot(4, 2, 7)
# axpos = F_plot.get_position(); axpos.x0 += 0.0; axpos.y0 -= 0.02; F_plot.set_position(axpos)

# empty lists for plotting
sim_times = []
zs = []; zdots = []
thetas = []; thetadots = []
Fs = []

t = P.t_start  # time starts at t_start
while t < P.t_end:
    if t <= 3:
        F = P.g
    elif t <= 5:
        F = 11
    elif t <= 7:
        F = 10

    y = ballbeam.update(F)  # Propagate the dynamics

    sim_times.append(t)
    Fs.append(F)
    zs.append(y[0][0])
    zdots.append(y[1][0])
    thetas.append(y[2][0])
    thetadots.append(y[3][0])

    animation.update(ballbeam.state)

    z_plot.clear(); zdot_plot.clear(); theta_plot.clear(); thetadot_plot.clear(); F_plot.clear()
    z_plot.plot(sim_times, zs, color='c', label='z'); z_plot.set_ylabel("Position (m)"); z_plot.legend(loc='upper right')
    zdot_plot.plot(sim_times, zdots, color='c', label='$\dot{z}$'); zdot_plot.set_ylabel('Velocity (m/s)'); zdot_plot.legend(loc='upper right')
    theta_plot.plot(sim_times, thetas, color='c', label='$\\theta$'); theta_plot.set_ylabel('Angle (rads)'); theta_plot.legend(loc='upper right')
    thetadot_plot.plot(sim_times, thetadots, color='c', label='$\dot{\\theta}$'); thetadot_plot.set_ylabel('Angular Velocity (rads/sec)')
    thetadot_plot.legend(loc='upper right')
    F_plot.plot(sim_times, Fs, color='c'); F_plot.set_ylabel('Force (N)')
    z_plot.grid(); zdot_plot.grid(); theta_plot.grid(); thetadot_plot.grid(); F_plot.grid()

    plt.pause(0.01)

    t += P.Ts
    if keyboard.is_pressed("q"): break
