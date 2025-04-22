# Inlet Design - Cycle Calc Integration
# Annika Carlson
# carlsoai@mail.uc.edu

import numpy as np
import sys, os
sys.path.append(os.getcwd() + r"/..")
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from scipy.optimize import fsolve
from utils.tbm import shockRelation
from utils.normshock import normshock
from utils.compflow import compflow
import utils.stdatmos as std
std = std.stdAtmos()
import utils.units as uu

def ramp(M0:float, M2:float, theta1:float, theta2:float, A1:float, A2:float,  pi_d:float):
    """
    Function to calculate various inlet conditions including capture area,
    additive drag, and pressure recovery for both the subsonic and supersonic cases.
    
    Parameters
    ----------
    M0 : float
        Freestream Mach number.
    M2 : float
        Fan face Mach number.
    theta1 : float
        First oblique shock deflection angle.
    theta2 : float
        Second oblique shock deflection angle.
    A1 : float
        Inlet throat area.
    A2 : float
        Fan face area.
    pi_d : float
        Diffuser total pressure ratio.
        
    Returns
    -------
    A0 : float
        Capture area.
    M1 : float
        Mach number at the throat.
    mdot1 : float
        Mass flow rate through A1.
    D_add : float
        Additive inlet drag.
    eta_r : float
        Total pressure loss across all shocks.
    """

    # set some values
    y = 1.4                 # specific heat ratio (cold)
    h = 40000.              # altitude (ft)
    R = 53.34               # gas constant (ft*lbf/lbm*R)
    T0 = std.T(h)           # freestream static temp (degF) (SA at 40,000 ft)
    P0 = std.P(h)           # freestream static pressure (psf) (SA at 40,000 ft)
    rho = std.rho(h)        # freestream density (slugs/ft^3)

    # freestream total conditions
    Tt0 = uu.degF2r(T0)/compflow(M0, y)[0]
    Pt0 = P0/compflow(M0, y)[1]*uu.psf2psi

    # calculate M1 to achieve design M2 and A1/A2 
    def equation(M1):
        return pi_d*(M2/M1)*((1 + ((y - 1)/2)*(M1**2))/(1 + ((y - 1)/2)*(M2**2)))**((y + 1)/(2 * (y - 1))) - (A1/A2)
    M1_guess = 0.1
    M1_design = fsolve(equation, M1_guess)[0]

    # supersonic cases
    if M0 > 1.0:
        # calculate the max turning angle of the flow before the first shock
        betas = np.linspace(1, 90, 50)
        thetas1 = np.array([shockRelation(M=M0, B=beta) for beta in betas])
        theta1_max = np.max(thetas1)

        # if 1st turn angle is too high, need a normal shock instead
        if theta1 > theta1_max:
            M01, Pty1_Ptx1, Py1_Px1, _, _ = normshock(M0, y)

            # set M1 from M2 and A1/A2 
            M1 = M1_design

            # calculate A0
            A0 = A1*(M1/M01)*((1 + ((y - 1)/2)*M01**2)/(1 + ((y - 1)/2)*M1**2))**((y + 1) / (2 * (y - 1)))

            # set total pressure loss across the shock
            eta_r = Pty1_Ptx1
            Pt1 = Pt0*eta_r

            # calculate inlet additive drag for the streamtube (psf)
            D_add = (P0*(Py1_Px1)*uu.psf2psi)*A1*(y*M1*
                    ((1 + 0.5 * (y - 1) * M01**2) / (1 + 0.5 * (y - 1) * M1**2))**(y / (2 * (y - 1))) *
                    (M1 * ((1 + 0.5 * (y - 1) * M01**2) / (1 + 0.5 * (y - 1) * M1**2))**0.5 - M01) +
                    ((1 + 0.5 * (y - 1) * M01**2) / (1 + 0.5 * (y - 1) * M1**2))**(y / (y - 1)) - 1)

        # if 1st turn angle is achievable, solve oblique shock
        else:
            # determine shock wave angle of the first oblique shock
            B1 = shockRelation(M=M0, theta=theta1, y=y)

            # get normal component of freestream mach number + solve like normal shock
            M1n_s1 = M0*np.sin(np.deg2rad(B1))
            M2n_s1, Pty1_Ptx1, _, _, _ = normshock(M1n_s1, y)

            # calculate actual mach number after the first shock, M01
            M01 = M2n_s1/np.sin(np.deg2rad(B1) - np.deg2rad(theta1))

            # calculate the max turning angle of the flow before the second shock
            thetas2 = np.array([shockRelation(M=M01, B=beta) for beta in betas])
            theta2_max = np.max(thetas2)

            # if 2nd turn angle is too high, need a normal shock instead
            if theta2 > theta2_max:
                M02, Pty2_Ptx2, Py2_Px2, _, _ = normshock(M01, y)

                # set M1 from M2 and A1/A2 
                M1 = M1_design

                # calculate A0
                A0 = A1*(M1/M02)*((1 + ((y - 1)/2)*M02**2)/(1 + ((y - 1)/2)*M1**2))**((y + 1) / (2 * (y - 1)))

                # calculate total pressure loss across all of the shocks
                eta_r = Pty1_Ptx1*Pty2_Ptx2
                Pt1 = Pt0*eta_r

                # calculate inlet additive drag for the streamtube (psf)
                D_add = (P0*(Py2_Px2)*uu.psf2psi)*A1*(y*M1*
                        ((1 + 0.5 * (y - 1) * M02**2) / (1 + 0.5 * (y - 1) * M1**2))**(y / (2 * (y - 1))) *
                        (M1 * ((1 + 0.5 * (y - 1) * M02**2) / (1 + 0.5 * (y - 1) * M1**2))**0.5 - M02) +
                        ((1 + 0.5 * (y - 1) * M02**2) / (1 + 0.5 * (y - 1) * M1**2))**(y / (y - 1)) - 1)

            # if 2nd turn angle is achievable, solve oblique shock
            else:
                # solve across the second oblique shock the same way as last
                B2 = shockRelation(M=M01, theta=theta2, y=y)
                M1n_s2 = M01*np.sin(np.deg2rad(B2))
                M2n_s2, Pty2_Ptx2, Py2_Px2, _, _ = normshock(M1n_s2, y)
                M02 = M2n_s2/np.sin(np.deg2rad(B2) - np.deg2rad(theta2))
    
                # check to see if flow becomes subsonic, and no normal shock occurs
                if M02 < 1:
                    # set M1 from M2 and A1/A2 
                    M1 = M1_design

                    # calculate A0
                    A0 = A1*(M1/M02)*((1 + ((y - 1)/2)*M02**2)/(1 + ((y - 1)/2)*M1**2))**((y + 1) / (2 * (y - 1)))

                    # calculate total pressure loss across all of the shocks
                    eta_r = Pty1_Ptx1*Pty2_Ptx2
                    Pt1 = Pt0*eta_r

                    # calculate inlet additive drag for the streamtube (psf)
                    D_add = (P0*(Py2_Px2)*uu.psf2psi)*A1*(y*M1*
                            ((1 + 0.5 * (y - 1) * M02**2) / (1 + 0.5 * (y - 1) * M1**2))**(y / (2 * (y - 1))) *
                            (M1 * ((1 + 0.5 * (y - 1) * M02**2) / (1 + 0.5 * (y - 1) * M1**2))**0.5 - M02) +
                            ((1 + 0.5 * (y - 1) * M02**2) / (1 + 0.5 * (y - 1) * M1**2))**(y / (y - 1)) - 1)

                # solve normal shock if flow is still supersonic
                else: 
                    # solve across the terminal normal shock
                    M1x = M02
                    M1y, Pty3_Ptx3, _, _, _ = normshock(M1x, y)

                    # final mach number after the terminal shock, going into the diffuser
                    M1 = M1y

                    # set A0
                    A0 = A1

                    # calculate total pressure loss across all of the shocks
                    eta_r = Pty1_Ptx1*Pty2_Ptx2*Pty3_Ptx3
                    Pt1 = Pt0*eta_r

                    # no additive drag
                    D_add = 0.0

        print(M01)
        print(M02)
        print(M1)

        # calculate mass flow rate at A1
        MFP1 = compflow(M1, y)[4]/1.28758
        mdot1 = (A1*(Pt1*MFP1))/np.sqrt(Tt0)  

    # subsonic cases
    else:
        # set M1 from M2 and A1/A2 
        M1 = M1_design

        # calculate mass flow rate at A1
        MFP1 = compflow(M1, y)[4]/1.28758
        mdot1 = (A1*(Pt0*MFP1))/np.sqrt(Tt0)
        
        # inlet additive drag (psf)
        D_add = (P0*uu.psf2psi)*A1*(y*M1*
                ((1 + 0.5 * (y - 1) * M0**2) / (1 + 0.5 * (y - 1) * M1**2))**(y / (2 * (y - 1))) *
                (M1 * ((1 + 0.5 * (y - 1) * M0**2) / (1 + 0.5 * (y - 1) * M1**2))**0.5 - M0) +
                ((1 + 0.5 * (y - 1) * M0**2) / (1 + 0.5 * (y - 1) * M1**2))**(y / (y - 1)) - 1)

        # calculate capture area, A0
        A0 = A1*(M1/M0)*((1 + ((y - 1)/2)*M0**2)/(1 + ((y - 1)/2)*M1**2))**((y + 1) / (2 * (y - 1)))

        eta_r = 1.0

    results = {
            "A0": A0,
            "M1": M1,
            "mdot1": mdot1,
            "D_add": D_add,
            "eta_r": eta_r,
            }

    return results

M0 = 2.0
M2 = 0.65
A2 = 1749.209
A1 = 1523.499
theta1 = 14.362
theta2 = 11.065
pi_d = 0.98

results = ramp(M0, M2, theta1, theta2, A1, A2, pi_d)
# print(results)
for key, value in results.items():
    print(f"{key}: {value}")


# y = 1.4
# # Define the equation
# def equation(M1):
#     return pi_d * (M2 / M1) * ((1 + ((y - 1)/2) * M1**2) / (1 + ((y - 1)/2) * M2**2))**((y + 1) / (2 * (y - 1))) - (A1 / A2)

# def Ax_Ay(Mx, My, Ptx_Pty):
#     Arat = (1/Ptx_Pty)*(My/Mx)*((1 + ((y-1)/2)*Mx**2)/(1 + ((y-1)/2)*My**2))**((y + 1)/(2*(y - 1)))
#     return Arat

# # Create a range of M1 values (avoid dividing by zero!)
# M1_vals = np.linspace(0.1, 3, 300)

# # Evaluate the function
# y_vals = Ax_Ay(M1_vals, M2, pi_d)

# # Plot
# import matplotlib.pyplot as plt
# plt.plot(M1_vals, y_vals, label=r'$f(M_1)$')
# plt.axhline(0, color='gray', linestyle='--')
# plt.xlabel('M1')
# plt.ylabel('f(M1)')
# plt.title('Plot of equation(M1)')
# plt.grid(True)
# plt.legend()
# plt.show()
