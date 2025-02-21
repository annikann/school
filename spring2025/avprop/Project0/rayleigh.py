# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         Rayleigh Flow        
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
from scipy.optimize import fsolve
from tabulate import tabulate

def phi_Msqd_Rayleigh(M, gamma):
    phi_Msqd = ((M**2)*(1 + ((gamma - 1)/2)*(M**2)))/((1 + gamma*(M**2))**2)
    return phi_Msqd

def P_Pstar_Rayleigh(M, gamma):
    P_Pstar = (1 + gamma) / (1 + gamma*M**2)
    return P_Pstar

def T_Tstar_Rayleigh(M, gamma):
    T_Tstar = ((M**2)*(1 + gamma)**2) / ((1 + gamma*(M**2))**2)
    return T_Tstar

def rho_rhostar_Rayleigh(M, gamma):
    rho_rhostar = (1/(M**2))*((1 + gamma*(M**2))/(1 + gamma))
    return rho_rhostar

def Po_Postar_Rayleigh(M, gamma):
    P_Pstar = P_Pstar_Rayleigh(M, gamma)
    Po_Postar = P_Pstar*((1+(M**2)*(gamma - 1)/2)**(gamma/(gamma - 1)))/(1 + (gamma - 1)/2)**(gamma/(gamma - 1))
    return Po_Postar

def To_Tostar_Rayleigh(M, gamma):
    To_Tostar = ((M**2)*(1 + gamma)**2 / (1 + gamma*(M**2))**2) * (1 + (M**2)*(gamma - 1)/2)/(1 + (gamma - 1)/2)
    return To_Tostar 

def V_Vstar_Rayleigh(M, gamma):
    P_Pstar = P_Pstar_Rayleigh(M, gamma)
    T_Tstar = T_Tstar_Rayleigh(M, gamma)
    V_Vstar = T_Tstar/P_Pstar 
    return V_Vstar

def Mach_Rayleigh(M1, gamma, To_Tostar):

    def equation(M, To_Tostar, gamma):
        return ((((gamma + 1) * (M**2)) / ((1 + gamma * (M**2))**2)) * (2 + (gamma - 1) * (M**2))) - To_Tostar

    def solve_M(To_Tostar, gamma):
        M_subsonic_guess = 0.1
        M_supersonic_guess = 1.01
        M_subsonic = fsolve(equation, M_subsonic_guess, args=(To_Tostar, gamma))[0]
        M_supersonic = fsolve(equation, M_supersonic_guess, args=(To_Tostar, gamma))[0]
        return np.array([M_subsonic, M_supersonic])
    
    return solve_M(To_Tostar, gamma)[0] if M1 < 1.0 else solve_M(To_Tostar, gamma)[1]

def RayleighFlow(M1, gamma, To_Tostar_2):

    # Compute Mach number or To_Tostar ratio based on inputs
    if To_Tostar_2 is None:
        M = M1
        To_Tostar = To_Tostar_Rayleigh(M, gamma)
    else:
        M = Mach_Rayleigh(M1, gamma, To_Tostar_2)
        To_Tostar = To_Tostar_2

    # Compute all ratios
    phi_Msqd = phi_Msqd_Rayleigh(M, gamma)
    P_Pstar = P_Pstar_Rayleigh(M, gamma)
    T_Tstar = T_Tstar_Rayleigh(M, gamma)
    rho_rhostar = rho_rhostar_Rayleigh(M, gamma)
    Po_Postar = Po_Postar_Rayleigh(M, gamma)
    V_Vstar = V_Vstar_Rayleigh(M, gamma)

    # Create a table with labeled values
    headers = ["Parameter", "Value"]
    table = [
        ["Mach Number (M)", f"{M:.4f}"],
        ["𝜙(M^2)", f"{phi_Msqd:.4f}"],
        ["P / P*", f"{P_Pstar:.4f}"],
        #["ρ / ρ*", f"{rho_rhostar:.4f}"],
        ["Po / Po*", f"{Po_Postar:.4f}"],
        ["T / T*", f"{T_Tstar:.4f}"],
        ["To / To*", f"{To_Tostar:.4f}"],
        ["V / V*", f"{V_Vstar:.4f}"],
    ]

    print("\nRayleigh Flow Results")
    print(tabulate(table, headers, tablefmt="grid", floatfmt=".4f"))

    results = [M, P_Pstar, T_Tstar, rho_rhostar, Po_Postar, To_Tostar, V_Vstar]

    return results



