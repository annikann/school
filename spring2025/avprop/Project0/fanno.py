# Fanno FLow Functions
# Annika Carlson
# carlsoai@mail.uc.edu

import numpy as np
from scipy.optimize import fsolve

def fanno(M1:float, y:float, cf:float, L:float, D:float):
    """
    Function to calculate change in flow parameters with friction (Fanno Flow).
    
    Parameters
    ----------
    M1 : float
        Mach number before heat addition.
    y : float
        Ratio of specific heats.
    cf : float
        Coefficient of friction.
    L : float
        Flow length.
    D : float
        Flow diameter.

    Returns
    -------
    Ms[i] : float
        Mach number.
    fLmax_D : float
        Friction factor at max length L to achieve sonic flow from M.
    I_Istar : float
        Sonic reference condition impulse function ratio.
    T_Tstar : float
        Sonic reference condition static temperature ratio.
    Po_Postar : float
        Sonic reference condition total pressure ratio.
    P_Pstar : float
        Sonic reference condition static pressure ratio.
    rho_rhostar : float
        Sonic reference condition density ratio.
    """
    # calculate friction factor for max L to become sonic
    if cf == None and L == None and D == None:
        fLmax_D_1 = ((1 - M1**2)/(y*(M1**2))) + ((y + 1)/(2*y))*np.log(((y + 1)*M1**2)/(2 + (y - 1)*M1**2))
        fLmax_D = fLmax_D_1
        Ms = [M1]

    # use cf, L and D values if given
    else:
        fL_D_act = (4*cf*L)/D
        fLmax_D_2 = fLmax_D_1 - fL_D_act
        fLmax_D = fLmax_D_2

        # solve for M2 based on calculated fLmax_D_2
        def equation(M2, fLmax_D_2, y):
            return ((1 - M2**2)/(y*(M2**2))) + ((y + 1)/(2*y))*np.log(((y + 1)*M2**2)/(2 + (y - 1)*M2**2)) - fLmax_D_2

        def solve_M(fLmax_D_2, y):
            M_subsonic_guess = 0.1
            M_supersonic_guess = 1.01
            M_subsonic = fsolve(equation, M_subsonic_guess, args=(fLmax_D_2, y))[0]
            M_supersonic = fsolve(equation, M_supersonic_guess, args=(fLmax_D_2, y))[0]
            return np.array([M_subsonic, M_supersonic])
        
        # return sub or supersonic based on M1
        M2 = solve_M(fLmax_D_2, y)[0] if M1 < 1.0 else solve_M(fLmax_D_2, y)[1]
        Ms = [M1, M2]

    results = []
    for i in range(len(Ms)):
        # calculate sonic impulse ratio
        I_Istar = (1 + y*(Ms[i]**2))/(Ms[i]*np.sqrt(2*(y + 1)*(1 + ((y - 1)/2)*(Ms[i]**2))))

        # calculate sonic pressure ratio
        P_Pstar = (1/Ms[i])*np.sqrt((y + 1)/(2 + (y - 1)*(Ms[i]**2)))

        # calculate sonic temperature ratio
        T_Tstar = (y + 1)/(2 + (y - 1)*(Ms[i]**2))

        # calculate sonic total pressure ratio
        Po_Postar = (1/Ms[i])*((2 + (y - 1)*(Ms[i]**2))/(y + 1))**((y + 1)/(2*(y - 1)))

        # calculate sonic density ratio
        rho_rhostar = (1/Ms[i])*np.sqrt((2 + (y - 1)*(Ms[i]**2))/(y + 1))

        results.append([Ms[i], fLmax_D, I_Istar, T_Tstar, Po_Postar, P_Pstar, rho_rhostar])

    return  results