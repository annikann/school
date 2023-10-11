import sys
sys.path.append('/Users/annikacarlson/Documents/school/controltheory/library')
import numpy as np 
import library.bbParam as P

class bbDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([[P.z0],          # z initial position
                               [P.zdot0],       # zdot initial velocity
                               [P.theta0],      # theta initial position
                               [P.thetadot0]])  # theta initial velocity
        
        # Set initial conditions
        self.Ts = P.Ts
        self.m1 = P.m1
        self.m2 = P.m2
        self.l = P.l
        self.r = P.r
        self.g = P.g
        self.force_limit = P.F_max       

    def update(self, F):
        # saturate the input force
        self.rk4_step(F)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y
    
    def f(self, state, F):
        # Set initial states
        z = state[0][0]
        zdot = state[1][0]
        theta = state[2][0]
        thetadot = state[3][0]
        # Create equations of motion.
        zddot = z*(thetadot**2) - self.g*np.sin(theta)
        thetaddot = (1/(((self.m2*(self.l**2))/3) + self.m1*(z**2)))*(F*self.l*np.cos(theta) - 2*self.m1*z*zdot*thetadot - \
                                                                    self.m1*self.g*z*np.cos(theta) - (((self.m2*self.g*self.l)/2)*np.cos(theta)))
        # Build and return xdot
        xdot = np.array([[zdot], [zddot], [thetadot], [thetaddot]])
        return xdot

    def h(self):
        z = self.state[0][0]
        zdot = self.state[1][0]
        theta = self.state[2][0]
        thetadot = self.state[3][0]
        y = np.array([[z], [zdot], [theta], [thetadot]])
        return y

    def rk4_step(self, F):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, F)
        F2 = self.f(self.state + self.Ts / 2 * F1, F)
        F3 = self.f(self.state + self.Ts / 2 * F2, F)
        F4 = self.f(self.state + self.Ts * F3, F)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        