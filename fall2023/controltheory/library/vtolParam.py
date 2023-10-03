# Inverted Pendulum Parameter File
import numpy as np

# Physical parameters 
mc = 1       # Mass of center body, kg
mr = 0.25    # Mass of rotor, kg
g = 9.8      # Gravity, m/s**2
Jc = 0.0042  # Moment of inertia, kg*m^2
d = 0.3      # Arm distance, m
mu = 0.1     # Something! kg/s

# parameters for animation
w = 0.5       # Width of the body, m
h = 0.5       # Height of the body, m

# Initial Conditions
z0 = 0.0                # ,m
zdot0 = 0.0             # ,m/s
h0 = 0.0                # ,m
hdot0 = 0.0             # ,m/s    
theta0 = 0.0            # ,rads
thetadot0 = 0.0         # ,rads/sec

# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 500.0  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.1  # the plotting and animation is updated at this rate

# saturation limits
F_max = 40.0                # Max Force, N

