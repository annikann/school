import numpy as np
from control.matlab import tf
import library.aerosonde_parameters as  P

class ComputeGains:
    def __init__(self):
        self.pn = P.states0[0][0]
        self.pe = P.states0[1][0]
        self.pd = P.states0[2][0]
        self.u = P.states0[3][0]
        self.v = P.states0[4][0]
        self.w = P.states0[5][0]
        self.phi = P.states0[6][0]
        self.theta = P.states0[7][0]
        self.psi = P.states0[8][0]
        self.p = P.states0[9][0]
        self.q = P.states0[10][0]
        self.r = P.states0[11][0]

        self.jx = P.jx
        self.jy = P.jy
        self.jz = P.jz
        self.jxz = P.jxz
        self.M = P.mass
        self.gravity = P.gravity
        self.Ts = P.ts_simulation

        self.S_wing = P.S_wing
        self.b = P.b
        self.c = P.c
        self.S_prop = P.S_prop
        self.rho = P.rho
        self.e = P.e
        self.AR = P.AR
        self.C_L_0 = P.C_L_0
        self.C_D_0 = P.C_D_0
        self.C_m_0 = P.C_m_0
        self.C_L_alpha = P.C_L_alpha
        self.C_D_alpha = P.C_D_alpha
        self.C_m_alpha = P.C_m_alpha
        self.C_L_q = P.C_L_q
        self.C_D_q = P.C_D_q
        self.C_m_q = P.C_m_q
        self.C_L_delta_e = P.C_L_delta_e
        self.C_D_delta_e = P.C_D_delta_e
        self.C_m_delta_e = P.C_m_delta_e
        self.m = P.M
        self.alpha0 = P.alpha0
        self.epsilon = P.epsilon
        self.C_D_p = P.C_D_p
        self.C_Y_0 = P.C_Y_0
        self.C_ell_0 = P.C_ell_0
        self.C_n_0 = P.C_n_0
        self.C_Y_beta = P.C_Y_beta
        self.C_ell_beta = P.C_ell_beta
        self.C_n_beta = P.C_n_beta
        self.C_Y_p = P.C_Y_p
        self.C_ell_p = P.C_ell_p
        self.C_n_p = P.C_n_p
        self.C_Y_r = P.C_Y_r
        self.C_ell_r = P.C_ell_r
        self.C_n_r = P.C_n_r
        self.C_Y_delta_a = P.C_Y_delta_a
        self.C_ell_delta_a = P.C_ell_delta_a
        self.C_n_delta_a = P.C_n_delta_a
        self.C_Y_delta_r = P.C_Y_delta_r
        self.C_ell_delta_r = P.C_ell_delta_r
        self.C_n_delta_r = P.C_n_delta_r
        self.C_prop = P.C_prop
        self.k_motor = P.k_motor
        self.k_T_p = P.k_T_p
        self.k_omega = P.k_omega

    def compute_tfs(self, x_trim, u_trim):

        Va_trim = np.sqrt(x_trim[3]**2 + x_trim[4]**2 + x_trim[5]**2)
        alpha_trim = np.arctan(x_trim[5]/x_trim[3])
        beta_trim = np.arctan(x_trim[4]/Va_trim)
        theta_trim = x_trim[7]
        gamma = self.jx * self.jz - self.jxz ** 2
        gamma3 = self.jz / gamma
        gamma4 = self.jxz / gamma
        Cpp = gamma3 * self.C_ell_p + gamma4 * self.C_n_p
        Cpdelta_a = gamma3 * self.C_ell_delta_a + gamma4 * self.C_n_delta_a

        # Calcule constants
        a_phi1 = -.5 * self.rho * Va_trim**2 * self.S_wing * self.b * Cpp * (self.b / (2 * Va_trim))
        a_phi2 = .5 * self.rho * Va_trim**2 * self.S_wing * self.b * Cpdelta_a
        a_theta1 = -(self.rho * Va_trim**2 * self.c * self.S_wing)/(2 * self.jy) * self.C_m_q * (self.c/(2 * Va_trim))
        a_theta2 = -(self.rho * Va_trim**2 * self.c * self.S_wing)/(2 * self.jy) * self.C_m_alpha
        a_theta3 = (self.rho * Va_trim**2 * self.c * self.S_wing)/(2 * self.jy) * self.C_m_delta_e

        a_V1 = ((self.rho * Va_trim * self.S_wing)/self.M) * (self.C_D_0 + (self.C_D_alpha * alpha_trim) \
                + (self.C_D_delta_e * u_trim[0])) + (self.rho * self.S_prop)/(self.M) * self.C_prop * Va_trim
        a_V2 = (self.rho * self.S_prop)/(self.M) * self.C_prop * self.k_motor**2 * u_trim[3]
        a_V3 = self.gravity * np.cos(theta_trim - alpha_trim)

        a_beta1 = -(self.rho * Va_trim * self.S_wing)/(2 * self.M) * self.C_Y_beta
        a_beta2 = (self.rho * Va_trim * self.S_wing)/(2 * self.M) * self.C_Y_delta_r

        # Calculate transfer functions
        T_phi_delta_a = tf([a_phi2], [1, a_phi1, 0])
        T_chi_phi = tf([self.gravity/Va_trim], [1, 0])
        T_theta_delta_e = tf([a_theta3], [1, a_theta1, a_theta2])
        T_h_theta = tf([Va_trim], [1, 0])
        T_h_Va = tf([theta_trim], [1, 0])
        T_Va_delta_t = tf([a_V2], [1, a_V1])
        T_Va_theta = tf([-a_V3], [1, a_V1])
        T_beta_delta_r = tf([a_beta2], [1, a_beta1])

        print('~~~~~~~ Open Loop Transfer Functions ~~~~~~~')
        print('T_phi_delta_a = ', T_phi_delta_a)
        print('T_chi_phi = ', T_chi_phi)
        print('T_theta_delta_e = ', T_theta_delta_e)
        print('T_h_theta = ', T_h_theta)
        print('T_beta_delta_r = ', T_beta_delta_r)
        print('T_Va_delta_t = ', T_Va_delta_t)
        print('T_Va_theta = ', T_Va_theta)
        print('T_h_Va = ', T_h_Va)

        return (T_phi_delta_a, T_chi_phi, T_theta_delta_e, T_h_theta, T_h_Va, T_Va_delta_t, T_Va_theta, T_Va_theta, T_beta_delta_r)


    def statespace(self, x_trim, u_trim):

        pn = x_trim.item(0)
        pe = x_trim.item(1)
        pd = x_trim.item(2)
        u = x_trim.item(3)
        v = x_trim.item(4)
        w = x_trim.item(5)
        phi = x_trim.item(6)
        theta = x_trim.item(7)
        psi = x_trim.item(8)
        p = x_trim.item(9)
        q = x_trim.item(10)
        r = x_trim.item(11)

        d_e = u_trim[0]
        d_a = u_trim[1]
        d_r = u_trim[2]
        d_t = u_trim[3]
        Va = np.sqrt(u**2 + v**2 + w**2)
        aoa = np.arctan(w/u)
        Beta = np.arctan(v/Va)

        gamma = self.jx * self.jz - self.jxz ** 2
        gamma1 = (self.jxz * (self.jx - self.jy + self.jz)) / gamma
        gamma2 = (self.jz * (self.jz - self.jy) + self.jxz ** 2) / gamma
        gamma3 = self.jz / gamma
        gamma4 = self.jxz / gamma
        gamma5 = (self.jz - self.jx) / self.jy
        gamma6 = self.jxz / self.jy
        gamma7 = ((self.jx - self.jy) * self.jx + self.jxz ** 2) / gamma
        gamma8 = self.jx / gamma

        C_p_0 = gamma3*self.C_ell_0 + gamma4*self.C_n_0
        C_p_Beta = gamma3*self.C_ell_beta + gamma4*self.C_n_beta
        C_p_p = gamma3*self.C_ell_p + gamma4*self.C_n_p
        C_p_r = gamma3*self.C_ell_r + gamma4*self.C_n_r
        C_p_delta_a = gamma3*self.C_ell_delta_a + gamma4*self.C_n_delta_a
        C_p_delta_r = gamma3*self.C_ell_delta_r + gamma4*self.C_n_delta_r
        C_r_0 = gamma4*self.C_ell_0 + gamma8*self.C_n_0
        C_r_Beta = gamma4*self.C_ell_beta + gamma8*self.C_n_beta
        C_r_p = gamma4*self.C_ell_p + gamma8*self.C_n_p
        C_r_r = gamma4*self.C_ell_r + gamma8*self.C_n_r
        C_r_delta_a = gamma4*self.C_ell_delta_a + gamma8*self.C_n_delta_a
        C_r_delta_r = gamma4*self.C_ell_delta_r + gamma8*self.C_n_delta_r

        Y_v = ((self.rho*self.S_wing*v)/(4*self.M*Va))*(self.C_Y_p*p + self.C_Y_r*r) \
              + ((self.rho*self.S_wing*v)/self.M)*(self.C_Y_0 + self.C_Y_beta*Beta + self.C_Y_delta_a*d_a + self.C_Y_delta_r*d_r) \
                + ((self.rho*self.S_wing*self.C_Y_beta)/(2*self.M))*np.sqrt(u**2 + w**2)
        Y_p = w + ((self.rho*Va*self.S_wing*self.b)/(4*self.M))*self.C_Y_p
        Y_r = -u + ((self.rho*Va*self.S_wing*self.b)/(4*self.M))*self.C_Y_r
        Y_delta_a = ((self.rho*Va**2*self.S_wing)/(2*self.M)) * self.C_Y_delta_a
        Y_delta_r = ((self.rho*Va**2*self.S_wing)/(2*self.M)) * self.C_Y_delta_r
        L_v = ((self.rho*self.S_wing*self.b**2*v)/(4*Va))*(C_p_p*p + C_p_r*r) \
              + (self.rho*self.S_wing*self.b*v)*(C_p_0 + C_p_Beta*Beta + C_p_delta_a*d_a + C_p_delta_r*d_r) \
                + (self.rho*self.S_wing*self.b*C_p_Beta/2)*np.sqrt(u**2 + w**2)
        L_p = gamma1*q + (self.rho*Va*self.S_wing*self.b**2/4)*C_p_p
        L_r = -gamma2*q + (self.rho*Va*self.S_wing*self.b**2/4)*C_p_r
        L_delta_a = (self.rho*Va**2*self.S_wing*self.b/2)*C_p_delta_a
        L_delta_r = (self.rho*Va**2*self.S_wing*self.b/2)*C_p_delta_r
        N_v = ((self.rho*self.S_wing*self.b**2*v)/(4*Va))*(C_r_p*p + C_r_r*r) \
              + (self.rho*self.S_wing*self.b*v)*(C_r_0 + C_r_Beta*Beta + C_r_delta_a*d_a + C_r_delta_r*d_r) \
                + (self.rho*self.S_wing*self.b*C_r_Beta/2)*np.sqrt(u**2 + w**2)
        N_p = gamma7*q + (self.rho*Va*self.S_wing*self.b**2 / 4)*C_r_p
        N_r = -gamma1*q + (self.rho*Va*self.S_wing*self.b**2 / 4)*C_r_r
        N_delta_a = (self.rho*Va**2*self.S_wing*self.b / 2)*C_r_delta_a
        N_delta_r = (self.rho*Va**2*self.S_wing*self.b / 2)*C_r_delta_r

        C_D = P.C_D_0 + (P.C_D_alpha * aoa)
        C_L = P.C_L_0 + (P.C_L_alpha * aoa)
        C_x_a = -P.C_D_alpha * np.cos(aoa) + P.C_L_alpha * np.sin(aoa)
        C_x_0 = -P.C_D_0 * np.cos(aoa) + P.C_L_0 * np.sin(aoa)
        C_x_d_e = -P.C_D_delta_e * np.cos(aoa) + P.C_L_delta_e * np.sin(aoa)
        C_x_q = -P.C_D_q * np.cos(aoa) + P.C_L_q * np.sin(aoa)
        C_Z = -C_D * np.sin(aoa) - C_L * np.cos(aoa)
        C_Z_q = -P.C_D_q * np.sin(aoa) - P.C_L_q * np.cos(aoa)
        C_Z_delta_e = -P.C_D_delta_e * np.sin(aoa) - P.C_L_delta_e * np.cos(aoa)
        C_Z_0 = -P.C_D_0 * np.sin(aoa) - P.C_L_0 * np.cos(aoa)
        C_Z_alpha = - P.C_D_alpha * np.sin(aoa) - P.C_L_alpha * np.cos(aoa)

        X_u = ((u * P.rho * P.S_wing) / P.M) * (C_x_0 + (C_x_a * d_a) + (C_x_d_e * d_e)) - ((P.rho * P.S_wing * w * C_x_a) / (2 * P.M)) \
            + ((P.rho * P.S_wing * P.c * C_x_q * u * q) / (4 * P.M * Va)) - ((P.rho * P.S_prop * P.C_prop *u) / P.M)
        X_w = -q + ((w * P.rho * P.S_wing) / P.M) * (C_x_0 + (C_x_a * d_a) + (C_x_d_e * d_e)) \
            + ((P.rho * P.S_wing * P.c * C_x_q * w * q) / (4 * P.M * Va)) + ((P.rho * P.S_wing * u * C_x_a) / (2 * P.M)) \
              - ((P.rho * P.S_prop * P.C_prop * w) / P.M)
        X_q = -w + ((P.rho * Va * P.S_wing * C_x_q * P.c) / (4 * P.M))
        X_delta_e = (P.rho * (Va ** 2) * P.S_wing * C_x_d_e) / (2 * P.M)
        X_delta_t = (P.rho * P.S_prop * P.C_prop * (P.k_motor ** 2) * d_t) / P.M
        Z_u = q + ((u * P.rho * P.S_wing) / (P.M)) * (C_Z_0 + (C_Z_alpha * aoa) + (C_Z_delta_e * d_e)) \
              - ((P.rho * P.S_wing * C_Z_alpha *w) / (2 * P.M)) + ((u * P.rho * P.S_wing * C_Z_q * P.c * q) / (4 * P.M * Va))
        Z_w = ((w * P.rho * P.S_wing) / (P.M)) * (C_Z_0 + (C_Z_alpha * aoa) + (C_Z_delta_e * d_e)) \
            + ((P.rho * P.S_wing * C_Z_alpha * u) / (2 * P.M)) + ((w * P.rho * P.S_wing * C_Z_q * P.c * q) / (4 * P.M * Va))
        Z_q = u + (P.rho * Va * P.S_wing * C_Z_q * P.c) / (4 * P.M)
        Z_delta_e = (P.rho * (Va ** 2) * P.S_wing * C_Z_delta_e) / (2 * P.M)
        M_u = ((u * P.rho * P.S_wing * P.c) / P.jy) * (P.C_m_0 + (P.C_m_alpha * aoa) + (P.C_m_delta_e * d_e)) \
            - ((P.rho * P.S_wing * P.c * P.C_m_alpha * w) / (2 * P.jy)) + ((P.rho * P.S_wing * (P.c ** 2) * P.C_m_q * q * u) / (4 * P.jy * Va))
        M_w = ((w * P.rho * P.S_wing * P.c) / P.jy) * (P.C_m_0 + P.C_m_alpha * aoa + P.C_m_delta_e * d_e) \
            + ((P.rho * P.S_wing * P.c * P.C_m_alpha * u) / (2 * P.jy)) + ((P.rho * P.S_wing * P.c ** 2 * P.C_m_q * q * w) / (4 * P.jy * Va))
        M_q = (P.rho * Va * P.c ** 2 * P.S_wing * P.C_m_q) / (4 * P.jy)
        M_delta_e = (P.rho * (Va ** 2) * P.S_wing * P.c * P.C_m_delta_e) / (2 * P.jy)

        Alat = np.array([[Y_v, Y_p, Y_r, self.gravity*np.cos(theta)*np.cos(phi), 0],
                      [L_v, L_p, L_r, 0, 0], 
                      [N_v, N_p, N_r, 0, 0],
                      [0, 1, np.cos(phi)*np.tan(theta), q*np.cos(phi)*np.tan(theta)-r*np.sin(phi)*np.tan(theta), 0],
                      [0, 0, np.cos(phi)*(1/np.cos(theta)), p*np.cos(phi)*(1/np.cos(theta)) - r*np.sin(phi)*(1/np.cos(theta)), 0]])
        Blat = np.array([[Y_delta_a, Y_delta_r], 
                         [L_delta_a, L_delta_r],
                         [N_delta_a, N_delta_r], 
                         [0, 0], 
                         [0, 0]])

        Along = np.array([[X_u, X_w, X_q, -self.gravity*np.cos(theta), 0],
                          [Z_u, Z_w, Z_q, -self.gravity*np.sin(theta), 0],
                          [M_u, M_w, M_q, 0, 0],
                          [0, 0, 1, 0, 0],
                          [np.sin(theta), -np.cos(theta), 0, u*np.cos(theta) + w*np.sin(theta), 0]])

        Blong = np.array([[X_delta_e, X_delta_t], 
                          [Z_delta_e, 0], 
                          [M_delta_e, 0], 
                          [0, 0], 
                          [0, 0]])

        print('Lateral A matrix: ', Alat)
        print('Lateral B matrix: ', Blat)
        print('Longitudinal A matrix: ', Along)
        print('Longitudinal B matrix: ', Blong)

        elatvalue, elatvect = np.linalg.eig(Alat)
        elongvalue, elongvect = np.linalg.eig(Along)
        print('Lateral eigenvalues: ', elatvalue)
        print('Longitudinal eigenvalues: ', elongvalue)

        def categorize_eigenvalues(eigenvalues):
            modes = []

            for eigenvalue in eigenvalues:
                real_part = np.real(eigenvalue)
                imag_part = np.imag(eigenvalue)

                # Determine mode based eigenvalues
                if np.isclose(real_part, 0.0, atol=1e-3):
                    if imag_part < 0:
                        modes.append("Phugoid Mode")
                    else:
                        modes.append("Spiral Mode")
                elif np.isclose(imag_part, 0.0, atol=1e-3):
                    if real_part < 0:
                        modes.append("Short Period Mode")
                    else:
                        modes.append("Roll Mode")
                else:
                    if real_part < 0:
                        modes.append("Dutch Roll Mode")
                    else:
                        modes.append("Other")

            return modes

        # Categorize eigenvalues of Alon
        alon_modes = categorize_eigenvalues(elongvalue)

        # Categorize eigenvalues of Alat
        alat_modes = categorize_eigenvalues(elatvalue)

        # Print the modes associated with Alon and Alat
        print("Modes associated with Alon eigenvalues:")
        for i, mode in enumerate(alon_modes):
            print(f"Eigenvalue {i + 1}: {mode}")

        print("\nModes associated with Alat eigenvalues:")
        for i, mode in enumerate(alat_modes):
            print(f"Eigenvalue {i + 1}: {mode}")

        return Alat, Blat, elatvalue, elongvalue