import library.smdParam as P
import numpy as np

class smdController:
    def __init__(self):
        self.m = P.m
        self.k = P.k
        self.b = P.b
        self.Fmax = P.Fmax

        ## ORIGINAL GAINS
        # self.kP = 3.05
        # self.kD = 7.2

        ## ITERATED
        self.kP = 5.9
        self.kD = 8.6

    def update(self, zc, state):
        z = state[0]
        zdot = state[1]
        feq = P.k*z
        fc = self.kP*(zc - z)  - self.kD*zdot
        f = feq + fc
        f = self.saturate(f, self.Fmax)
        return f[0]
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u