# Wind Generation

import numpy as np
from library.rotations import Rvb
import control.matlab as mat
import library.aerosonde_parameters as P

class wind():

    def __init__(self, Vs):
        self.Vs = Vs

    # Dryden Gust Function
    def drydenGust(self, Va, t):
        # Set values from Dryden gust model parameters
        Lu = Lv = 200
        Lw = 50
        sigmau = sigmav = 1.06
        sigmaw = 0.7

        # Transfer functions
        au=sigmau*np.sqrt(2*Va/Lu)
        av=sigmav*np.sqrt(3*Va/Lv)
        aw=sigmaw*np.sqrt(3*Va/Lw)
        Hu = mat.tf([0, au],[1, Va/Lu])
        Hv = mat.tf([av, av*Va/(np.sqrt(3)*Lv)],[1, 2*Va/Lv, (Va/Lv)**2])
        Hw = mat.tf([aw, aw*Va/(np.sqrt(3)*Lw)],[1, 2*Va/Lw, (Va/Lw)**2])

        # Add white noise and solve
        T=[0, t]
        white_noise_u = np.random.normal(0,1,1)
        white_noise_v = np.random.normal(0,1,1)
        white_noise_w = np.random.normal(0,1,1)
        y_u, _, _ = mat.lsim(Hu, white_noise_u[0], T)
        y_v, _, _ = mat.lsim(Hv, white_noise_v[0], T)
        y_w, _, _ = mat.lsim(Hw, white_noise_w[0], T)

        # gust components
        wg_u = y_u[1]
        wg_v = y_v[1]
        wg_w = y_w[1]

        return np.array([[wg_u],[wg_v],[wg_w]])

    def windout(self, states, Va, sim_time):
        pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = states.flatten()
        gusts = self.drydenGust(Va, sim_time)
        windtot = Rvb(phi, theta, psi)*self.Vs + gusts
        Var = np.array([[u - windtot[0][0]],
                        [v - windtot[1][0]],
                        [w - windtot[2][0]]])
        ur, vr, wr = Var.flatten()
        Va = np.sqrt(ur**2 + vr**2 + wr**2)
        alpha = np.arctan(wr/ur)
        beta = np.arcsin(vr/Va)

        print(ur, vr, wr)

        return Va, alpha, beta