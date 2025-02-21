import rayleigh

M1 = 0.3
To1 = 500.0
Po1 = 600.0
qdot = 500.0
gamma_a = 1.4
cp_a = 1.004
gamma_b = 1.325
cp_b = 1.171

To2a = To1 + (qdot/cp_a)
results1a = rayleigh.RayleighFlow(M1, gamma_a, None)
To_Tostar_2a = (results1a[5])*(To2a/To1)
results2a = rayleigh.RayleighFlow(M1, gamma_a, To_Tostar_2a)

To2b = To1 + (qdot/cp_b)
results1b = rayleigh.RayleighFlow(M1, gamma_b, None)
To_Tostar_2b = (results1b[5])*(To2b/To1)
results2b = rayleigh.RayleighFlow(M1, gamma_b, To_Tostar_2b)