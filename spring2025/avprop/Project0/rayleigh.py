# Rayleigh FLow Functions
# Annika Carlson
# carlsoai@mail.uc.edu

import numpy as np
from scipy.optimize import fsolve

def Mach_Rayleigh(M1:float, y:float, To_Tostar:float):
    """
    Function to calculate mach number after heat addition.
    """
    def equation(M, To_Tostar, y):
        return ((((y + 1) * (M**2)) / ((1 + y * (M**2))**2)) * (2 + (y - 1) * (M**2))) - To_Tostar

    def solve_M(To_Tostar, y):
        M_subsonic_guess = 0.1
        M_supersonic_guess = 1.01
        M_subsonic = fsolve(equation, M_subsonic_guess, args=(To_Tostar, y))[0]
        M_supersonic = fsolve(equation, M_supersonic_guess, args=(To_Tostar, y))[0]
        return np.array([M_subsonic, M_supersonic])
    
    return solve_M(To_Tostar, y)[0] if M1 < 1.0 else solve_M(To_Tostar, y)[1]

def rayleigh(M1, y, To_Tostar_2):
    """
    Function to calculate change in flow parameters with heat addition (Rayleigh Flow).
    
    Parameters
    ----------
    M1 : float
        Mach number before heat addition.
    y : float
        Specific heat ratio.
    To_Tostar_2 : float
        Sonic reference condition total temperature ratio, after heat addition.

    Returns
    -------
    M : float
        Mach number.
    phi_Msqd : float
        Entropy change parameter.
    P_Pstar : float
        Sonic reference condition static pressure ratio.
    T_Tstar : float
        Sonic reference condition static temperature ratio.
    rho_rhostar : float
        Sonic reference condition density ratio.
    Po_Postar : float
        Sonic reference condition total pressure ratio.
    To_Tostar: float
        Sonic reference condition total pressure ratio.
    V_Vstar : float
        Sonic reference condition velocity ratio.
    """
    # Compute Mach number or To_Tostar ratio based on inputs
    if To_Tostar_2 is None:
        M = M1
        To_Tostar = ((M**2)*(1 + y)**2 / (1 + y*(M**2))**2) * (1 + (M**2)*(y - 1)/2)/(1 + (y - 1)/2)
    else:
        M = Mach_Rayleigh(M1, y, To_Tostar_2)
        To_Tostar = To_Tostar_2

    # Compute all ratios
    phiMsqd = ((M**2)*(1 + ((y - 1)/2)*(M**2)))/((1 + y*(M**2))**2)
    P_Pstar = (1 + y)/(1 + y*M**2)
    T_Tstar = ((M**2)*(1 + y)**2) / ((1 + y*(M**2))**2)
    rho_rhostar = (1/(M**2))*((1 + y*(M**2))/(1 + y))
    Po_Postar = P_Pstar*((1+(M**2)*(y - 1)/2)**(y/(y - 1)))/(1 + (y - 1)/2)**(y/(y - 1))
    V_Vstar = T_Tstar/P_Pstar 

    return M, np.round(phiMsqd, 6), np.round(To_Tostar, 6), np.round(T_Tstar, 6), np.round(Po_Postar, 6), \
            np.round(P_Pstar, 6), np.round(rho_rhostar, 6), np.round(V_Vstar, 6)



