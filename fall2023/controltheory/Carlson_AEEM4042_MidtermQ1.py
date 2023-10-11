import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
import matplotlib.pyplot as plt
import numpy as np
import library.bbParam as P
from library.signalGenerator import signalGenerator
from library.bbAnimation import bbAnimation
import keyboard

# Initialize state array
state = np.array([[P.z0],
                  [P.zdot0],
                  [P.theta0],
                  [P.thetadot0]])

# Call animation
animation = bbAnimation(0.7, True)

# Signal inputs
z_val = signalGenerator(amplitude=0.5, frequency=0.1) 
theta_val = signalGenerator(amplitude=np.deg2rad(45), frequency=0.01)

# Create lists for plots, add subplots
sim_times = []
zs = []
thetas = []
z_plot = animation.fig.add_subplot(2, 2, 2)
theta_plot = animation.fig.add_subplot(2, 2, 4)

print('Press Q to end simulation')
t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    t_next_plot = t + P.t_plot
    
    while t < t_next_plot:
        state[0][0] = z_val.sin(t)
        if state[0][0] < 0: state[0][0] = -z_val.sin(t)
        state[2][0] = theta_val.sin(t)

        zs.append(state[0][0])
        thetas.append(state[2][0])
        sim_times.append(t)

        animation.update(state)

        z_plot.clear(); theta_plot.clear()
        z_plot.plot(sim_times, zs, color='c')
        z_plot.set_title('Position')
        z_plot.set_ylabel("z (m)")
        z_plot.grid()
        theta_plot.plot(sim_times, thetas, color = 'c')
        theta_plot.set_title('Beam Angle')
        theta_plot.set_ylabel("theta (rads)")
        theta_plot.grid()

        t = t + P.Ts  # advance time by Ts

    plt.pause(0.001)  # allows time for animation to draw

    if keyboard.is_pressed('q'): break

