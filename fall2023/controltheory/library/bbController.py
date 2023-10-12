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
        # self.filter = zeroCancelingFilter()

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

        # tmp = self.kPz*(zr - z) - self.kDz*zdot
        # thetar = self.filter.update(tmp)
        # F = self.kPth*(thetar - theta) - self.kDth*thetadot
        # F = self.saturate(F, self.Fmax)

        thetar = self.kPz*(zr - z) - self.kDz*zdot + zr
        F = self.kPth*(thetar - theta) - self.kDth*thetadot
        F = self.saturate(F, self.Fmax)
        return F

    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u
    
# class zeroCancelingFilter:
#     def __init__(self):
#         self.a = 2.652
#         self.b = 2.652
#         self.state = 0.0

#     def update(self, input):
#         # integrate using RK1
#         self.state = self.state \
#                     + P.Ts*(-self.b*self.state + self.a*input)
#         return self.state