import sys
sys.path.append('/Users/annikacarlson/Documents/FlightMech/flight_sim')
import numpy as np 
from numpy import sin as s
from numpy import cos as c
from numpy import tan as t
import library.aerosonde_parameters as P

class mavDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = P.states0
        # Setting parameters
        self.Ts = P.ts_simulation
        self.g = P.gravity
        self.mass = P.mass
        self.jx = P.jx
        self.jy = P.jy
        self.jz = P.jz
        self.jxz = P.jxz
        self.gamma = P.gamma
        self.gamma1 = P.gamma1
        self.gamma2 = P.gamma2
        self.gamma3 = P.gamma3
        self.gamma4 = P.gamma4
        self.gamma5 = P.gamma5
        self.gamma6 = P.gamma6
        self.gamma7 = P.gamma7
        self.gamma8 = P.gamma8   

    def update(self, fx, fy, fz, l, m, n):
        # saturate the input force
        self.rk4_step(fx, fy, fz, l, m, n)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y
    
    def f(self, state, fx, fy, fz, l, m, n):
        # Set initial states
        pn = state[0][0]
        pe = state[1][0]
        pd = state[2][0]
        u = state[3][0]
        v = state[4][0]
        w = state[5][0]
        phi = state[6][0]
        theta = state[7][0]
        psi = state[8][0]
        p = state[9][0]
        q = state[10][0]
        r = state[11][0]

        eom1 = np.array([[c(theta)*c(psi), s(phi)*s(theta)*c(psi) - c(phi)*s(psi), c(phi)*s(theta)*c(psi) + s(phi)*s(psi)],
                         [c(theta)*s(psi), s(phi)*s(theta)*s(psi) + c(phi)*c(psi), c(phi)*s(theta)*s(psi) - s(phi)*c(psi)],
                         [-s(theta),       s(phi)*c(theta),                        c(phi)*c(theta)]])
        eom1x = np.array([[u], [v], [w]])
        eom1s = np.dot(eom1, eom1x)

        pndot = [[eom1s[0][0]]]
        pedot = [[eom1s[1][0]]]
        pddot = [[eom1s[2][0]]]
        
        udot = np.array([r*v - q*w]) + np.array([fx/self.mass])
        vdot = np.array([p*w - r*u]) + np.array([fy/self.mass])
        wdot = np.array([q*u - p*w]) + np.array([fz/self.mass])

        eom3 = np.array([[1, s(phi)*t(theta), c(phi)*t(theta)],
                         [0, c(phi), -s(phi)],
                         [0, s(phi)/c(theta), c(phi)/c(theta)]])
        eom3x = np.array([[p], [q], [r]])
        eom3s = np.dot(eom3, eom3x)

        phidot = [[eom3s[0][0]]]
        thetadot = [[eom3s[1][0]]]
        psidot = [[eom3s[2][0]]]

        pdot = np.array([self.gamma1*p*q - self.gamma2*q*r]) + np.array([self.gamma3*l + self.gamma4*n])
        qdot = np.array([self.gamma5*p*r - self.gamma6*(p**2 - r**2)]) + np.array([(m/self.jy)])
        rdot = np.array([self.gamma7*p*q - self.gamma1*q*r]) + np.array([self.gamma4*l + self.gamma8*n])

        # Build and return xdot
        xdot = np.array([[pndot[0][0]], 
                         [pedot[0][0]], 
                         [pddot[0][0]], 
                         [udot[0]], 
                         [vdot[0]], 
                         [wdot[0]], 
                         [phidot[0][0]], 
                         [thetadot[0][0]], 
                         [psidot[0][0]], 
                         [pdot[0]], 
                         [qdot[0]], 
                         [rdot[0]]], dtype=float)
        return xdot   

    def h(self):
        pn = self.state[0][0]
        pe = self.state[1][0]
        pd = self.state[2][0]
        u = self.state[3][0]
        v = self.state[4][0]
        w = self.state[5][0]
        phi = self.state[6][0]
        theta = self.state[7][0]
        psi = self.state[8][0]
        p = self.state[9][0]
        q = self.state[10][0]
        r = self.state[11][0]
        y = np.array([[pn], [pe], [pd], [u], [v], [w], [phi], [theta], [psi], [p], [q], [r]], dtype=float)
        return y

    def rk4_step(self, fx, fy, fz, l, m, n):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, fx, fy, fz, l, m, n)
        F2 = self.f(self.state + self.Ts / 2 * F1, fx, fy, fz, l, m, n)
        F3 = self.f(self.state + self.Ts / 2 * F2, fx, fy, fz, l, m, n)
        F4 = self.f(self.state + self.Ts * F3, fx, fy, fz, l, m, n)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        