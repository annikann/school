# Inverted Pendulum Parameter File
import numpy as np

# Physical parameters 
mr = 0.25    # Mass of rotor, kg
g = 9.81     # Gravity, m/s**2

# Variable parameters
mc = np.random.uniform(0.8*1, 1.2*1)          # Mass of center body, kg
Jc = np.random.uniform(0.8*0.0042, 1.2*0.0042)     # Moment of inertia, kg*m^2
d = np.random.uniform(0.8*0.3, 1.2*0.3)         # Arm distance, m
mu = np.random.uniform(0.8*0.1, 1.2*0.1)         # Something! kg/s

kIth = 0.5
kIz = 0.5
kIh = 0.0
kDh = 0.5833
kPh = 0.1134
kDz = -0.03286
kPz = -0.007709
kDth = 0.191
kPth = 0.372

# parameters for animation
w = 0.5       # Width of the body, m
h = 0.5       # Height of the body, m

# Initial Conditions
z0 = 0.0                # ,m
zdot0 = 0.0             # ,m/s
h0 = 0.1                # ,m
hdot0 = 0.0             # ,m/s    
theta0 = 0.0            # ,rads
thetadot0 = 0.0         # ,rads/sec

# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 500.0  # End time of simulation
Ts = 0.1  # sample time for simulation
t_plot = 0.1  # the plotting and animation is updated at this rate

# saturation limits
F_max = 10.0                # Max Force, N
theta_max = 0.5*np.pi

sigma = 0.05
beta = (2.0*sigma - Ts)/(2.0*sigma + Ts)

