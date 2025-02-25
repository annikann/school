import rayleigh
import fanno

M1 = 0.3
To1 = 500.0
Po1 = 600.0
qdot = 500.0
y_a = 1.4
cp_a = 1.004
y_b = 1.325
cp_b = 1.171

To2a = To1 + (qdot/cp_a)
results1a = rayleigh.RayleighFlow(M1, y_a, None)
To_Tostar_2a = (results1a[6])*(To2a/To1)
results2a = rayleigh.RayleighFlow(M1, y_a, To_Tostar_2a)

# To2b = To1 + (qdot/cp_b)
# results1b = rayleigh.RayleighFlow(M1, gamma_b, None)
# To_Tostar_2b = (results1b[6])*(To2b/To1)
# results2b = rayleigh.RayleighFlow(M1, gamma_b, To_Tostar_2b)

# M1 = 0.6        # inlet Mach number
# P1 = 10.0       # inlet static pressure (psia)
# T1 = 500.0      # inlet static temperature (degR)
# L = 8.0         # duct length (ft)
# D = 1.0         # duct height (ft)
# cf = 0.004      # coefficient of friction
# y = 1.4

# results = fanno.FannoFlow(M1, cf, L, D, y)


# headers = ["Parameter", "Value"]
# table = [
#     ["Mach Number (M)", f"{Ms[i]:.4f}"],
#     ["4cfL* / D", f"{fLmax_D}"],
#     ["I / I*", f"{I_Istar:.4f}"],
#     ["P / P*", f"{P_Pstar:.4f}"],
#     ["Po / Po*", f"{Po_Postar:.4f}"],
#     ["T / T*", f"{T_Tstar:.4f}"],
#     ["ρ / ρ*", f"{rho_rhostar:.4f}"],
# ]

# print("\nFanno Flow Results")
# print('Station ' + str(i + 1) + ':')
# print(tabulate(table, headers, tablefmt="grid", floatfmt=".4f"))

# headers = ["Parameter", "Value"]
# table = [
#     ["Mach Number (M)", f"{M:.4f}"],
#     ["𝜙(M^2)", f"{phiMsqd:.4f}"],
#     ["P / P*", f"{P_Pstar:.4f}"],
#     ["ρ / ρ*", f"{rho_rhostar:.4f}"],
#     ["Po / Po*", f"{Po_Postar:.4f}"],
#     ["T / T*", f"{T_Tstar:.4f}"],
#     ["To / To*", f"{To_Tostar:.4f}"],
#     ["V / V*", f"{V_Vstar:.4f}"],
# ]

# print("\nRayleigh Flow Results")
# print(tabulate(table, headers, tablefmt="grid", floatfmt=".4f"))