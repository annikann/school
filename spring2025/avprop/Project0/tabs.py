# Simple One-Dimensional Flow Functions
# Slade Brooks & Annika Carlson
# brooksl@mail.uc.edu | carlsoai@mail.uc.edu

import numpy as np
import os
import glob
from tabulate import tabulate
from compflow import compflow
from normshock import normshock
from rayleigh import rayleigh
from fanno import fanno

current_folder = os.path.dirname(os.path.abspath(__file__))
results_folder = os.path.join(current_folder, "simpleflowtabs")
os.makedirs(results_folder, exist_ok=True)
files = glob.glob(os.path.join(results_folder, "*"))
for file in files:
    if os.path.isfile(file):
        os.remove(file)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       Compressible Flow      
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# range of mach numbers and gammas
Ms = np.arange(0., 4.001, 0.01)
ys = np.array([1.4, 1.33, 1.3])
T_Tts = np.empty((len(ys), len(Ms)))
P_Pts = np.empty((len(ys), len(Ms)))
rho_rhots = np.empty((len(ys), len(Ms)))
A_Astars = np.empty((len(ys), len(Ms)))
MFPsRgcs = np.empty((len(ys), len(Ms)))
mus = np.empty((len(ys), len(Ms)))
vs = np.empty((len(ys), len(Ms)))

# loop through and calculate parameters
for i, y in enumerate(ys):
    for k, M in enumerate(Ms):
        T_Tts[i, k], P_Pts[i, k], rho_rhots[i, k], A_Astars[i, k], MFPsRgcs[i, k], mus[i, k], vs[i, k] = compflow(M, y)

    # output tables
    data = list(zip(Ms, T_Tts[i, :], P_Pts[i, :], rho_rhots[i, :], A_Astars[i, :], MFPsRgcs[i, :], mus[i, :], vs[i, :]))
    headers = ["M", "T/Tt", "P/Pt", "ρ/ρt", "A/A*", "MFP√(R/gc)", "μ", "v"]
    units = ["-", "-", "-", "-", "-", "W-s/m-√T", "deg", "deg"]
    header_units = [f"{h}\n({u})" for h, u in zip(headers, units)]
    fulldata = [header_units] + data
    colalign = ['center', 'center', 'center', 'center', 'center', 'center', 'center', 'center']
    output_path = os.path.join(results_folder, f"Compflow_y{ys[i]:.2f}_table.txt") 
    with open(output_path, "w", encoding="utf-8") as file:
        file.write(f"-----===== Compressible Flow Tables (gamma={ys[i]:.2f}) =====-----\n")
        file.write(tabulate(fulldata, headers="firstrow", colalign=colalign, tablefmt="fancy_grid", \
                            floatfmt=(".2f", ".6f", ".6f", ".6f", ".6f", ".6f", ".6f", ".6f")))    

# ~~~~~~~~~~~~~~~~~~~~~~~~
#       Normal Shock      
# ~~~~~~~~~~~~~~~~~~~~~~~~

# range of mach numbers and gammas
Mxs = np.arange(1., 4.001, 0.01)
yxs = np.array([1.4, 1.33, 1.3])
Mys = np.empty((len(ys), len(Ms)))
Pty_Ptxs = np.empty((len(ys), len(Ms)))
Py_Pxs = np.empty((len(ys), len(Ms)))
rhoy_rhoxs = np.empty((len(ys), len(Ms)))
Ty_Txs = np.empty((len(ys), len(Ms)))

# loop through and calculate parameters
for i, yx in enumerate(yxs):
    for k, Mx in enumerate(Mxs):
        Mys[i, k], Pty_Ptxs[i, k], Py_Pxs[i, k], rhoy_rhoxs[i, k], Ty_Txs[i, k] = normshock(Mx, yx)

    # output tables
    data = list(zip(Mxs, Mys[i, :], Pty_Ptxs[i, :], Py_Pxs[i, :], rhoy_rhoxs[i, :], Ty_Txs[i, :]))
    headers = ["Mx", "My", "Pty/Ptx", "Py/Px", "ρy/ρx", "Ty/Tx"]
    fulldata = [headers] + data
    colalign = ['center', 'center', 'center', 'center', 'center', 'center']
    output_path = os.path.join(results_folder, f"Normshock_y{ys[i]:.2f}_table.txt") 
    with open(output_path, "w", encoding="utf-8") as file:
        file.write(f"-----===== Normal Shock Tables (gamma={yxs[i]:.2f}) =====-----\n")
        file.write(tabulate(fulldata, headers="firstrow", colalign=colalign, tablefmt="fancy_grid", \
                            floatfmt=(".2f", ".2f", ".6f", ".6f", ".6f", ".6f")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~
#       Rayleigh Flow      
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# range of mach numbers and gammas
Ms = np.arange(0., 4.001, 0.01)
ys = np.array([1.4])
phiMsqds = np.empty((len(ys), len(Ms)))
Tt_Ttstar = np.empty((len(ys), len(Ms)))
T_Tstar = np.empty((len(ys), len(Ms)))
Pt_Ptstar = np.empty((len(ys), len(Ms)))
P_Pstar = np.empty((len(ys), len(Ms)))
rho_rhostar = np.empty((len(ys), len(Ms)))
V_Vstar = np.empty((len(ys), len(Ms)))

# loop through and calculate parameters
for i, y in enumerate(ys):
    for k, M in enumerate(Ms):
        _, phiMsqds[i, k], Tt_Ttstar[i, k], T_Tstar[i, k], Pt_Ptstar[i, k], P_Pstar[i, k], _, _ = rayleigh(M, y, None)

    # output tables
    data = list(zip(Ms, phiMsqds[i, :], Tt_Ttstar[i, :], T_Tstar[i, :], Pt_Ptstar[i, :], P_Pstar[i, :]))
    headers = ["M", "𝜙(M^2)", "Tt/Ttstar", "T/Tstar", "Pt/Ptstar", "P/Pstar"]
    fulldata = [headers] + data
    colalign = ['center', 'center', 'center', 'center', 'center', 'center']
    output_path = os.path.join(results_folder, f"Rayleighflow_y{ys[i]:.2f}_table.txt") 
    with open(output_path, "w", encoding="utf-8") as file:
        file.write(f"-----===== Rayleigh Flow Tables (gamma={ys[i]:.2f}) =====-----\n")
        file.write(tabulate(fulldata, headers="firstrow", colalign=colalign, tablefmt="fancy_grid", \
                            floatfmt=(".2f", ".6f", ".6f", ".6f", ".6f", ".6f")))  

# ~~~~~~~~~~~~~~~~~~~~~~
#       Fanno Flow      
# ~~~~~~~~~~~~~~~~~~~~~~

# range of mach numbers and gammas
Ms = np.arange(0., 4.001, 0.01)
ys = np.array([1.4])
fLmax_D = np.empty((len(ys), len(Ms)))
I_Istar = np.empty((len(ys), len(Ms)))
T_Tstar = np.empty((len(ys), len(Ms)))
Po_Postar = np.empty((len(ys), len(Ms)))
P_Pstar = np.empty((len(ys), len(Ms)))

# loop through and calculate parameters
for i, y in enumerate(ys):
    for k, M in enumerate(Ms):
        _, fLmax_D[i, k], I_Istar[i, k], T_Tstar[i, k], Po_Postar[i, k], P_Pstar[i, k], _ = fanno(M, y, None, None, None)[0]

    # output tables
    data = list(zip(Ms, fLmax_D[i, :], I_Istar[i, :], T_Tstar[i, :], Po_Postar[i, :], P_Pstar[i, :]))
    headers = ["M", "4cfL*/D", "I/I*", "T/T*", "Po/Po*", "P/P*"]
    fulldata = [headers] + data
    colalign = ['center', 'center', 'center', 'center', 'center', 'center']
    output_path = os.path.join(results_folder, f"Fannoflow_y{ys[i]:.2f}_table.txt") 
    with open(output_path, "w", encoding="utf-8") as file:
        file.write(f"-----===== Fanno Flow Tables (gamma={ys[i]:.2f}) =====-----\n")
        file.write(tabulate(fulldata, headers="firstrow", colalign=colalign, tablefmt="fancy_grid", \
                            floatfmt=(".2f", ".6f", ".6f", ".6f", ".6f", ".6f")))  