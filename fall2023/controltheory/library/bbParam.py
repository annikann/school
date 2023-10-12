# Inverted Pendulum Parameter File
import numpy as np

# Physical parameters 
m1 = 0.35    # Mass of ball, kg
m2 = 2.0     # Mass of beam, kg
l = 0.5      # Length of beam, m
g = 9.8      # Gravity, m/s**2
r = 0.015    # Radius of the ball, m

# Initial Conditions
z0 = 0.5                # ,m
zdot0 = 0.0             # ,m/s
theta0 = 0.0            # ,rads
thetadot0 = 0.0         # ,rads/sec

# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 500.0  # End time of simulation
Ts = 0.1       # sample time for simulation
t_plot = 0.1   # the plotting and animation is updated at this rate

Fmax = 15    # N

