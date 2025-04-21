# Inlet Design - Cycle Calc Integration
# Annika Carlson
# carlsoai@mail.uc.edu

import numpy as np
import sys, os
sys.path.append(os.getcwd() + r"/..")
from scipy.optimize import fsolve
from utils.tbm import shockRelation
from utils.normshock import normshock
from utils.compflow import compflow
import utils.stdatmos as std
std = std.stdAtmos()
import utils.units as uu

def ramp(M0:float, M2:float, theta1:float, theta2:float, A1:float, A2:float):
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
    Pt0 = P0/compflow(0.9, y)[1]*uu.psf2psi

    # supersonic cases
    if M0 > 1.0:
        # determine shock wave angle of the first oblique shock
        B1 = shockRelation(M=M0, theta=theta1, y=y)

        # get normal component of freestream mach number + solve like normal shock
        M1n_s1 = M0*np.sin(np.deg2rad(B1))
        M2n_s1, Pty1_Ptx1, _, _, _ = normshock(M1n_s1, y)

        # calculate actual mach number after the first shock, M01
        M01 = M2n_s1/np.sin(np.deg2rad(B1) - np.deg2rad(theta1))

        # solve across the second oblique shock the same way
        B2 = shockRelation(M=M01, theta=theta2, y=y)
        M1n_s2 = M0*np.sin(np.deg2rad(B2))
        M2n_s2, Pty2_Ptx2, _, _, _ = normshock(M1n_s2, y)
        M02 = M2n_s2/np.sin(np.deg2rad(B2) - np.deg2rad(theta2))

        # solve across the terminal normal shock
        M1x = M02
        M1y, Pty3_Ptx3, _, _, _ = normshock(M1x, y)

        # final mach number after the terminal shock, going into the diffuser
        M1 = M1y

        # calculate total pressure loss across all of the shocks
        eta_r = Pty1_Ptx1*Pty2_Ptx2*Pty3_Ptx3

        # calculate mass flow rate at A1
        MFP1 = compflow(M1, y)[4]/1.28758
        mdot1 = (A1*(Pt0*MFP1))/np.sqrt(Tt0)

        D_add = 0.0
        A0 = A1

    # subsonic cases
    else:
        # calculate M1 from M2 and A1/A2 
        def equation(M1):
            return (M2/M1)*((1 + ((y - 1)/2)*M1**2)/(1 + ((y - 1)/2)*M2**2))**((y + 1) / (2 * (y - 1))) - (A1/A2)
        M1_guess = 0.1
        M1 = fsolve(equation, M1_guess)[0]

        # calculate mass flow rate at A1
        MFP1 = compflow(M1, y)[4]/1.28758
        mdot1 = (A1*(Pt0*MFP1))/np.sqrt(Tt0)
        
        # nondimensionalized inlet additive drag
        D_add_p0A1 = ( y * M1 *
                        ((1 + 0.5 * (y - 1) * M0**2) / (1 + 0.5 * (y - 1) * M1**2))**(y / (2 * (y - 1))) *
                        (M1 * ((1 + 0.5 * (y - 1) * M0**2) / (1 + 0.5 * (y - 1) * M1**2))**0.5 - M0) +
                        ((1 + 0.5 * (y - 1) * M0**2) / (1 + 0.5 * (y - 1) * M1**2))**(y / (y - 1)) - 1)

        # dimensionalize that guy
        D_add = D_add_p0A1*(P0*uu.psf2psi)*A1   # psf

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