# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#         Fanno Flow        
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
from scipy.optimize import fsolve
from tabulate import tabulate

def FannoFlow(M1:float, cf:float, L:float, D:float, y:float):
    """
    Function to calculate change in flow parameters with friction (Fanno Flow).
    
    Parameters
    ----------
    M1 : float
        Mach number before heat addition.
    cf : float
        Coefficient of friction.
    L : float
        Flow length.
    D : float
        Flow diameter.
    y : float
        Ratio of specific heats.

    Returns
    -------
    Ms[i] : float
        Mach number.
    P_Pstar : float
        Sonic reference condition static pressure ratio.
    T_Tstar : float
        Sonic reference condition static temperature ratio.
    Po_Postar : float
        Sonic reference condition total pressure ratio.
    rho_rhostar : float
        Sonic reference condition density ratio.
    """
    fL_D_act = (4*cf*L)/D
    fLmax_D_1 = ((1 - M1**2)/(y*(M1**2))) + ((y + 1)/(2*y))*np.log(((y + 1)*M1**2)/(2 + (y - 1)*M1**2))
    fLmax_D_2 = fLmax_D_1 - fL_D_act

    def equation(M2, fLmax_D_2, y):
        return ((1 - M2**2)/(y*(M2**2))) + ((y + 1)/(2*y))*np.log(((y + 1)*M2**2)/(2 + (y - 1)*M2**2)) - fLmax_D_2

    def solve_M(fLmax_D_2, y):
        M_subsonic_guess = 0.1
        M_supersonic_guess = 1.01
        M_subsonic = fsolve(equation, M_subsonic_guess, args=(fLmax_D_2, y))[0]
        M_supersonic = fsolve(equation, M_supersonic_guess, args=(fLmax_D_2, y))[0]
        return np.array([M_subsonic, M_supersonic])
    
    M2 = solve_M(fLmax_D_2, y)[0] if M1 < 1.0 else solve_M(fLmax_D_2, y)[1]
    
    Ms = [M1, M2]
    results = []
    for i in range(len(Ms)):
        I_Istar = (1 + y*(Ms[i]**2))/(Ms[i]*np.sqrt(2*(y + 1)*(1 + ((y - 1)/2)*(Ms[i]**2))))
        P_Pstar = (1/Ms[i])*np.sqrt((y + 1)/(2 + (y - 1)*(Ms[i]**2)))
        T_Tstar = (y + 1)/(2 + (y - 1)*(Ms[i]**2))
        Po_Postar = (1/Ms[i])*((2 + (y - 1)*(Ms[i]**2))/(y + 1))**((y + 1)/(2*(y - 1)))
        rho_rhostar = (1/Ms[i])*np.sqrt((2 + (y - 1)*(Ms[i]**2))/(y + 1))

        if i == 0: fLmax_D = fLmax_D_1
        else: fLmax_D = fLmax_D_2

        headers = ["Parameter", "Value"]
        table = [
            ["Mach Number (M)", f"{Ms[i]:.4f}"],
            ["4cfL* / D", f"{fLmax_D}"],
            ["I / I*", f"{I_Istar:.4f}"],
            ["P / P*", f"{P_Pstar:.4f}"],
            ["Po / Po*", f"{Po_Postar:.4f}"],
            ["T / T*", f"{T_Tstar:.4f}"],
            ["ρ / ρ*", f"{rho_rhostar:.4f}"],
        ]

        print("\nFanno Flow Results")
        print('Station ' + str(i + 1) + ':')
        print(tabulate(table, headers, tablefmt="grid", floatfmt=".4f"))

        results.append([Ms[i], P_Pstar, T_Tstar, Po_Postar, rho_rhostar])

    return  results