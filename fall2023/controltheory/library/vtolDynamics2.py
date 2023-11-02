import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
import numpy as np 
import library.vtolParam2 as P

class vtolDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.z0],         # z initial position
            [P.zdot0],      # zdot initial velocity
            [P.h0],         # h initial position
            [P.hdot0],      # hdot initial velocity
            [P.theta0],     # theta initial position
            [P.thetadot0]   # theta initial velocity
        ])
        # simulation time step
        self.Ts = P.Ts
        self.mc = P.mc
        self.mr = P.mr
        self.Jc = P.Jc
        self.d = P.d
        self.mu = P.mu
        self.g = P.g
        self.force_limit = P.F_max       

    def update(self, fr, fl):
        # saturate the input force
        self.rk4_step(fr, fl)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y
    
    def f(self, state, fr, fl):
        # Set initial states
        z = state[0][0]
        zdot = state[1][0]
        h = state[2][0]
        hdot = state[3][0]
        theta = state[4][0]
        thetadot = state[5][0]
        # Create equations of motion.
        zddot = (-1*(fr + fl)*np.sin(theta) - self.mu*zdot)/(self.mc + 2*self.mr)
        # hddot = ((fr + fl)*np.cos(theta) - (self.mc + 2*self.mr)*self.g)/(self.mc + 2*self.mr)
        hddot = ((fr + fl)*np.cos(theta))/(self.mc + 2*self.mr) - self.g
        thetaddot = (self.d*(fr - fl))/(self.Jc + 2*self.mr*self.d**2)
        # Build and return xdot
        xdot = np.array([[zdot], [zddot], [hdot], [hddot], [thetadot], [thetaddot]])
        return xdot

    def h(self):
        z = self.state[0][0]
        h = self.state[2][0]
        theta = self.state[4][0]
        y = np.array([[z], [h], [theta]])
        return y

    def rk4_step(self, fr, fl):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, fr, fl)
        F2 = self.f(self.state + self.Ts / 2 * F1, fr, fl)
        F3 = self.f(self.state + self.Ts / 2 * F2, fr, fl)
        F4 = self.f(self.state + self.Ts * F3, fr, fl)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        