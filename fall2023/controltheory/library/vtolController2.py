import library.vtolParam as P
import numpy as np

class vtolController2:
    def __init__(self):
        self.mr = P.mr
        self.mc = P.mc
        self.d = P.d
        self.g = P.g
        self.mu = P.mu
        self.jc = P.Jc
        self.d = P.d
        self.kPh = 0.1134
        self.kDh = 0.5833
        self.kPth = 0.3721
        self.kDth = 0.1913
        self.kPz = -0.00771
        self.kDz = -0.03285
        self.Fmax = P.F_max

    def update(self, zc, hc, state):
        z, zdot, h, hdot, theta, thetadot = state.flatten()
        
        # Theta
        taueq = 0.
        zs = self.kPz*(zc - z) - self.kDz*zdot
        tauc = self.kPth*(zs - theta)  - self.kDth*thetadot
        tau = taueq + tauc
        
        # Height
        feq = self.g*(self.mc + 2*self.mr)
        fc = self.kPh*(hc - h)  - self.kDh*hdot
        f = feq + fc
        
        # Forces
        fr = (tau + self.d*f)/(2*self.d)
        fl = f - fr

        fr, fl = self.saturate(fr, fl)

        return fr, fl
    
    def saturate(self, u1, u2):
        if abs(u1) > self.Fmax:
            u1 = self.Fmax*np.sign(u1)
        if abs(u2) > self.Fmax:
            u2 = self.Fmax*np.sign(u2)
        return u1, u2