# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#         Fanno Flow        
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
from scipy.optimize import fsolve
from tabulate import tabulate

def T_Tstar_Fanno(M, gamma):
    T_Tstar = (gamma + 1)/(2 + (gamma - 1)*(M**2))
    return T_Tstar

def P_Pstar_Fanno(M, gamma):
    P_Pstar = (1/M)*np.sqrt((gamma + 1)/(2 + (gamma - 1)*(M**2)))
    return P_Pstar

def rho_rhostar_Fanno(M, gamma):
    rho_rhostar = (1/M)*np.sqrt((2 + (gamma - 1)*(M**2))/(gamma + 1))
    return rho_rhostar

def Po_Postar_Fanno(M, gamma):
    Po_Postar = (1/M)*((2 + (gamma - 1)*(M**2))/(gamma + 1))**((gamma + 1)/(2*(gamma - 1)))
    return Po_Postar

def FannoFlow(M1, cf, L, D, gamma):

    fL_D_act = (4*cf*L)/D
    fLmax_D_1 = ((1 - M1**2)/(gamma*(M1**2))) + ((gamma + 1)/(2*gamma))*np.log(((gamma + 1)*M1**2)/(2 + (gamma - 1)*M1**2))
    fLmax_D_2 = fLmax_D_1 - fL_D_act

    def equation(M2, fLmax_D_2, gamma):
        return ((1 - M2**2)/(gamma*(M2**2))) + ((gamma + 1)/(2*gamma))*np.log(((gamma + 1)*M2**2)/(2 + (gamma - 1)*M2**2)) - fLmax_D_2

    def solve_M(fLmax_D_2, gamma):
        M_subsonic_guess = 0.1
        M_supersonic_guess = 1.01
        M_subsonic = fsolve(equation, M_subsonic_guess, args=(fLmax_D_2, gamma))[0]
        M_supersonic = fsolve(equation, M_supersonic_guess, args=(fLmax_D_2, gamma))[0]
        return np.array([M_subsonic, M_supersonic])
    
    M2 = solve_M(fLmax_D_2, gamma)[0] if M1 < 1.0 else solve_M(fLmax_D_2, gamma)[1]
    
    Ms = [M1, M2]
    results = []
    for i in range(len(Ms)):
        P_Pstar = P_Pstar_Fanno(Ms[i], gamma)
        T_Tstar = T_Tstar_Fanno(Ms[i], gamma)
        Po_Postar = Po_Postar_Fanno(Ms[i], gamma)
        rho_rhostar = rho_rhostar_Fanno(Ms[i], gamma)

        headers = ["Parameter", "Value"]
        table = [
            ["Mach Number (M)", f"{Ms[i]:.4f}"],
            ["P / P*", f"{P_Pstar:.4f}"],
            ["T / T*", f"{T_Tstar:.4f}"],
            ["Po / Po*", f"{Po_Postar:.4f}"],
            ["ρ / ρ*", f"{rho_rhostar:.4f}"],
        ]

        print("\nFanno Flow Results")
        print('Station ' + str(i + 1) + ':')
        print(tabulate(table, headers, tablefmt="grid", floatfmt=".4f"))

        results.append([Ms[i], P_Pstar, T_Tstar, Po_Postar, rho_rhostar])

    return  results