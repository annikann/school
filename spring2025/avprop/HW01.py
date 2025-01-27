# Annika Carlson
# Adv. Topics in AV Propulsion
# Module 1 Assignment - Code

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ~~~~~~~~~~~~~~
#   Question 4
# ~~~~~~~~~~~~~~

## Create list of Mach numbers + set known values
Ms = np.arange(0.3, 2.1, 0.1).round(1)
a = 40000       # altitude, ft
ns = [1, 4]     # load factors
gc = 32.174     # gravitational constant, ft*lbm/lbf*s^2

# HF-1 Aircraft Data 
Wmax = 40000    # max gross TOW, lbf
W = 0.9*Wmax    # 90% MaxGTOW
S = 720         # wing area, sq ft

# From Altitude tables
y = 1.4                 # specific heat ratio
d = 0.1858              # pressure ratio at 40kft, P/Pstd
Pstd = 2116.2           # std day reference pressure, lbf/ft^2
th = 0.7519             # std day temperature ratio at 40kft
aspeed = th*1116        # speed of sound at 40kft, ft/s
rho = 0.07647*(d/th)    # density at 40kft, lbm/ft^3

## Set up given data from table 1.4 to interpolate K1 and CD0 values
Mnums = np.array([0.0, 0.8, 1.2, 1.4, 2.0])
K1vals = np.array([0.2, 0.2, 0.2, 0.25, 0.4])
CD0vals = np.array([0.012, 0.012, 0.02267, 0.028, 0.027])

# Create interpolating functions
K1_interp = interp1d(Mnums, K1vals, kind='linear')
CD0_interp = interp1d(Mnums, CD0vals, kind='linear')

# # Interpolate K1 and CD0 values
# interpolated_K1 = K1_interp(Ms)
# interpolated_CD0 = CD0_interp(Ms)

# print(interpolated_CD0)

## Loop through each Mach number to calculate:
#   dynamic pressure, drag coefficient, lift coefficient, L, D and L/D
qs = []
K1s = []
CD0s = []
CDs = []
CLs = []
Ls = [[],[]]
Ds = [[],[]]
LDs = [[],[]]

for j in range(len(ns)):
    for i in range(len(Ms)):
        q = 0.5*y*d*Pstd*(Ms[i])**2
        qs.append(q)

        if Ms[i] not in Mnums:
            # Interpolate K1 and CD0 values
            interpolated_K1 = K1_interp(Ms[i])
            interpolated_CD0 = CD0_interp(Ms[i])
            K1s.append(interpolated_K1.item())
            CD0s.append(interpolated_CD0.item())
        else:
            if Ms[i] <= 0.8:
                K1s.append(0.2)
                CD0s.append(0.012)
            elif Ms[i] == 1.2:
                K1s.append(0.2)
                CD0s.append(0.02267)
            elif Ms[i] == 1.4:
                K1s.append(0.25)
                CD0s.append(0.028)
            elif Ms[i] == 2.0:
                K1s.append(0.4)
                CD0s.append(0.027)
    
        L = ns[j]*W
        Ls[j].append(L)

        CL = (ns[j]*W)/(q*S)
        # CL = L/(q*S*gc)
        CD = K1s[i]*CL**2 + CD0s[i]
        CDs.append(CD)

        # D = CDs[i]*q*S
        D = CD*q*S
        Ds[j].append(D)

        LDs[j].append(L/D)

plt.figure('Lift vs Mach Number')
plt.plot(Ms, Ls[0], label='Load Factor = 1')
plt.plot(Ms, Ls[1], label='Load Factor = 4')
plt.xlabel("Mach Number")
plt.ylabel("Lift, L (lbf)")
plt.title("Lift vs Mach Number")
plt.tight_layout(); plt.grid(); plt.legend()

plt.figure('Drag vs Mach Number')
plt.plot(Ms, Ds[0], label='Load Factor = 1')
plt.xlabel("Mach Number")
plt.ylabel("Drag, D (lbf)")
plt.title("HF-1 Level Flight (n=1) Drag")
plt.tight_layout(); plt.grid(); plt.legend()

# plt.figure('Drag vs Mach Number')
# plt.plot(Ms, Ds[1], label='Load Factor = 4')
# plt.xlabel("Mach Number")
# plt.ylabel("Drag, D (lbf)")
# plt.title("HF-1 (n=4) Drag")
# plt.tight_layout(); plt.grid(); plt.legend()

plt.figure('L/D vs Mach Number')
plt.plot(Ms, LDs[0], label='Load Factor = 1')
plt.plot(Ms, LDs[1], label='Load Factor = 4')
plt.xlabel("Mach Number")
plt.ylabel("L/D")
plt.title("Lift to Drag Ratio vs Mach Number")
plt.tight_layout(); plt.grid(); plt.legend()

plt.show()

# print(qs)
# print(K1s,CD0s)
