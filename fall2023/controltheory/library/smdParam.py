# Inverted Pendulum Parameter File
import numpy as np

# Physical parameters
g = 9.8    # Gravity, m/s**2
# m = 5      # Mass, kg
# k = 3      # Spring coefficient
# b = 0.5    # Damping coefficient, Ns

# Randomize parameters
m = np.random.uniform(0.8*5, 1.2*5)         # mass of mass (kg)
k = np.random.uniform(0.8*3, 1.2*3)         # spring constant (N/m)
b = np.random.uniform(0.8*0.5, 1.2*0.5)     # damping coefficient (Ns/m)

# parameters for animation
w = 0.5          # Width of the cart, m
h = 0.5          # Height of the cart, m
gap = 0.005      # Gap between the cart and x-axis

# Initial Conditions
f0 = 3                  # ,N
z0 = 0.0                # ,m
zdot0 = 0.0             # ,m/s

# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 1000.0  # End time of simulation
Ts = 0.1  # sample time for simulation
t_plot = 0.1  # the plotting and animation is updated at this rate

# saturation limits
Fmax = 6.0                # Max Force, N

