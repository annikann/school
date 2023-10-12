import library.bbParam as P
import numpy as np

class bbController:
    def __init__(self):
        self.m1 = P.m1
        self.m2 = P.m2
        self.l = P.l
        self.Fmax = P.Fmax
        self.thetamax = 30*np.pi/180

        self.kPth = 1.825
        self.kDth = 1.173
        self.kPz = -0.00494
        self.kDz = -0.0317
        self.kDC = 1.

    def update(self, zr, state):
        z = state[0][0]
        zdot = state[1][0]
        theta = state[2][0]
        thetadot = state[3][0]

        # Equilibrium force
        feq = (P.m1*P.g*z + P.m2*P.g*(0.5*P.l))/P.l
        thetar = self.kPz*(zr - z) - self.kDz*zdot
        fc = self.kPth*(thetar - theta) - self.kDth*thetadot

        f = feq + fc
        F = self.saturate(f, self.Fmax)
        return F
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u
    