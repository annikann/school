import library.bbParam as P
import numpy as np

class bbController:
    def __init__(self):
        self.m1 = P.m1
        self.m2 = P.m2
        self.l = P.l
        self.Fmax = P.Fmax

        ## ORIGINAL GAINS
        # self.kPth = 1.825
        # self.kDth = 1.173
        # self.kPz = -0.004949
        # self.kDz = -0.031743

        self.kPth = 1.825
        self.kDth = 1.173
        self.kPz = -0.0064 #last -0.007
        self.kDz = -0.044  #last -0.05

    def update(self, zr, state):
        z = state[0][0]
        zdot = state[1][0]
        theta = state[2][0]
        thetadot = state[3][0]

        # Equilibrium force
        feq = (P.m1*P.g*z / self.l) + (P.m2*P.g / 2)
        thetar = self.kPz*(zr - z) - self.kDz*zdot
        fc = self.kPth*(thetar - theta) - self.kDth*thetadot
        f = feq + fc
        F = self.saturate(f, self.Fmax)
        return F
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u
    