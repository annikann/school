import library.bbParam as P
import numpy as np

class bbController:
    def __init__(self):
        self.m1 = P.m1
        self.m2 = P.m2
        self.l = P.l
        self.Fmax = P.Fmax

        self.kPth = 1.825
        self.kDth = 1.173
        self.kPz = -0.00494
        self.kDz = -0.0317
        self.kDC = 1.
        self.t_rz = 10.
        self.t_rth = 1.

    # SPLIT BALL AND BEAM??
    def update(self, zr, state):
        z = state[0][0]
        zdot = state[1][0]
        theta = state[2][0]
        thetadot = state[3][0]

        # Figure out how these work
        # feq = P.k*z
        # fc = self.kPz*(zr - z)  - self.kD*zdot
        # f = feq + fc

        # Think this is all wrong, see pg 119 in the book
        # fball = (P.m1*P.g*z)/P.l
        # fc = self.kPz*(zr - z) - self.kDz*zdot
        # f = fball + fc

        thetar = self.kPz*(zr - z) - self.kDz*zdot
        F = self.kPth*(thetar - theta) - self.kDth*thetadot
        F = self.saturate(F, self.Fmax)
        return F
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u