import library.vtolParam3 as P
import numpy as np
from library.vtolControllerPID2 import PIDControl

class control:
    def __init__(self):
        self.mc = P.mc
        self.mr = P.mr
        self.kpt = P.kPth
        self.kdt = P.kDth
        self.kpz = P.kPz
        self.kdz = P.kDz
        self.kph = P.kPh
        self.kdh = P.kDh
        self.d = P.d
        self.g = P.g
        self.force_limit = P.F_max

        self.torqueCtrl = PIDControl(P.kPth, P.kIth, P.kDth, P.F_max, P.beta, P.Ts)
        self.thetaCtrl = PIDControl(P.kPz, P.kIz, P.kDz, P.theta_max, P.beta, P.Ts)
        self.forceCtrl = PIDControl(P.kPh, P.kIh, P.kDh, P.F_max, P.beta, P.Ts)

        self.limit = P.F_max

    def update(self, hr, zr, q):
        z = q[0][0]
        h = q[1][0]
        theta = q[2][0]

        tau_eq = 0.
        theta_r = self.thetaCtrl.PD(zr, z)
        tau_s = self.torqueCtrl.PID(theta_r, theta)
        tau = tau_eq + tau_s

        F_eq = (self.mc + (2.0*self.mr)) * self.g
        F_s = self.forceCtrl.PID(hr, h)
        F = F_eq + F_s

        fr = (tau + self.d*F)/(2*self.d)
        fl = F - fr

        fr = self.saturate(fr)
        fl = self.saturate(fl)

        return fr, fl
    
    def saturate(self, u):
        if abs(u) > self.force_limit:
            u = self.force_limit*np.sign(u)
        return u