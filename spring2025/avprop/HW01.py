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

## Set up given data from table 1.4 to interpolate K1 and CD0 values
Mnums = np.array([0.0, 0.8, 1.2, 1.4, 2.0])
K1vals = np.array([0.2, 0.2, 0.2, 0.25, 0.4])
CD0vals = np.array([0.012, 0.012, 0.02267, 0.028, 0.027])

# Create interpolating functions
K1_interp = interp1d(Mnums, K1vals, kind='linear')
CD0_interp = interp1d(Mnums, CD0vals, kind='linear')

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

plt.figure("Fig 1")
plt.plot(Ms, Ls[0], label='Lift')
plt.plot(Ms, Ds[0], label='Drag')
plt.xlabel("Mach Number"); plt.xlim(0.3, 2.0), plt.xticks(np.arange(0.3, 2.1, 0.1))
plt.ylabel("Force (lbf)"); plt.ylim(0.0, 40000.0)
plt.title("HF-1 Lift & Drag at 40kft, 90% Max GTOW (n=1)")
plt.tight_layout(); plt.grid(); plt.legend()

plt.figure("Fig 2")
plt.plot(Ms, Ls[1], label='Lift')
plt.plot(Ms, Ds[1], label='Drag')
plt.xlabel("Mach Number"); plt.xlim(0.3, 2.0), plt.xticks(np.arange(0.3, 2.1, 0.1))
plt.ylabel("Force (lbf)"); plt.ylim(0.0, 250000.0)
plt.title("HF-1 Lift & Drag at 40kft, 90% Max GTOW (n=4)")
plt.tight_layout(); plt.grid(); plt.legend()

plt.figure("Fig 3")
plt.plot(Ms, LDs[0], label='Load Factor = 1')
plt.plot(Ms, LDs[1], label='Load Factor = 4')
plt.xlabel("Mach Number"); plt.xlim(0.3, 2.0), plt.xticks(np.arange(0.3, 2.1, 0.1))
plt.ylabel("L/D"); plt.ylim(0.0, 12.0)
plt.title("HF-1 Lift to Drag Ratio at 40kft, 90% Max GTOW")
plt.tight_layout(); plt.grid(); plt.legend()

# ~~~~~~~~~~~~~~
#   Question 5
# ~~~~~~~~~~~~~~

## Set known values
M = 0.8
a = 40000       # altitude, ft
n = 1           # load factors
gc = 32.174     # gravitational constant, ft*lbm/lbf*s^2
g0 = gc

# HF-1 Aircraft Data 
Wmax = 40000    # max gross TOW, lbf
Wi = 0.9*Wmax   # initial weight, 90% MaxGTOW, lbf
S = 720         # wing area, sq ft

# From Altitude tables
y = 1.4                 # specific heat ratio
d = 0.1858              # pressure ratio at 40kft, P/Pstd
Pstd = 2116.2           # std day reference pressure, lbf/ft^2
th = 0.7519             # std day temperature ratio at 40kft
a = 1116*np.sqrt(th)    # speed of sound at 40kft, ft/s
rho = 0.07647*(d/th)    # density at 40kft, lbm/ft^3

## First, calculate CL and CD
q = 0.5*y*d*Pstd*(M**2)

# from table 1.4 for Mach 0.8:
K1 = 0.20
K2 = 0.0
CD0 = 0.012

CL = (n*Wi)/(q*S)
CD = K1*(CL**2) + (K2*CL) + CD0

## Calculate TSFC from equation 1.36b:
TSFC = (1.0 + 0.35*M)*np.sqrt(th)     # (lbm/h)/lbf

## Calculate the Range Factor
V = M*a*(3600/6080)             # nm/h
RF = (CL/CD)*(V/TSFC)*(gc/g0)   # nm

## Calculate s for a range of weight fractions for both cruise climb and level cruise
Wfrac = np.arange(0.1, 1.05, 0.05).round(2)
Wfrac_CC = []
Wfrac_LC = []
s_CC = []
s_LC = []
s05 = []

for i in range(len(Wfrac)):
    s_CC_val = -RF*np.log(Wfrac[i])
    s_LC_val = 2*RF*(1 - (np.sqrt(Wfrac[i])))

    if Wfrac[i] == 0.5:
        s05.append(s_CC_val)
        s05.append(s_LC_val)

    s_CC.append(s_CC_val)
    s_LC.append(s_LC_val)

plt.figure("Fig 4")
plt.plot(Wfrac, s_CC, label='Cruise Climb')
plt.scatter(0.5, s05[0], label='Cruise Climb Range @ $W_f$/$W_i$ = 0.5')
plt.plot(Wfrac, s_LC, label='Level Cruise')
plt.scatter(0.5, s05[1], label='Level Cruise Range @ $W_f$/$W_i$ = 0.5')
plt.xlabel("Weight Fraction, $W_f$/$W_i$"); plt.xlim(0.1, 1.0), plt.xticks(np.arange(0.1, 1.1, 0.1))
plt.ylabel("Range, s (ft)"); plt.ylim(0.0, 10000.0)
plt.title("HF-1 Range at 40kft, 90% Max GTOW")
plt.tight_layout(); plt.grid(); plt.legend()

plt.show()

