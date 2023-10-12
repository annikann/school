import library.bbParam as P

class bbController:
    def __init__(self):
        self.m1 = P.m1
        self.m2 = P.m2
        self.Fmax = P.Fmax

        self.kPtheta = 1.825
        self.kDtheta = 1.173
        self.KPz = -0.00494
        self.kDz = -0.0317
        self.kDC = 1.
        self.t_rz = 10.
        self.t_rtheta = 1.

    def update(self, zr, thetar, state):
        z = state[0][0]
        zdot = state[1][0]
        theta = state[2][0]
        thetadot = state[3][0]
        feq = P.k*z
        fc = self.kP*(zc - z)  - self.kD*zdot
        f = feq + fc
        return f[0]