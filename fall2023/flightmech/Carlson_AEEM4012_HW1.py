# Annika Carlson
# Flight Mechanics, AEEM4012
# Assignment 1 - UAV Animation

import sys
sys.path.append('/Users/annikacarlson/Documents/flightmech/library')
import numpy as np
from math import cos, sin, tan
import scipy.linalg as linalg
import library.simulation_parameters as SIM
from library.cube_animation import cube_animation
from library.signalGenerator import signalGenerator
import keyboard as key
from library.sampleUAV_animation import UAV_animation

state=np.array([[0], [0], [-1], [0], [0], [0], [0], [0], [0], [0], [0], [0]])
UAV_anim = UAV_animation(state, scale=5, multfigs=True)

# initialize the simulation and signal generator
sim_time = SIM.start_time
angles = signalGenerator(np.deg2rad(40), 1)
translations = signalGenerator(10, 1)

print("Press Q to exit simulation")
pn=state[0,0]
pe=state[1,0]
pd=state[2,0]
phi=state[6,0]
theta=state[7,0]
psi=state[8,0]

# create subplots
angles_plot = UAV_anim.fig.add_subplot(2, 2, 2)
trans_plot = UAV_anim.fig.add_subplot(2, 2, 4)
sim_times = []
phis = []
thetas = []
psis = []
pns = []
pes = []
pds = []

while sim_time < SIM.end_time:

    # rotations
    if sim_time <= 1:
        phi = 0
    elif sim_time <= 2:
        phi = angles.sin(sim_time)
    elif sim_time <= 3:
        theta = -angles.sin(sim_time)
    elif sim_time <= 4:
        psi = angles.sin(sim_time)
    # translations
    elif sim_time <= 5:
        pn = translations.sin(sim_time)
    elif sim_time <= 6:
        pe = translations.sin(sim_time)
    elif sim_time <= 7:
        pd = translations.sin(sim_time)

    # make lists for plotting
    sim_times.append(sim_time)
    phis.append(phi)
    thetas.append(theta)
    psis.append(psi)
    pns.append(pn)
    pes.append(pe)
    pds.append(pd)

    # create plots
    angles_plot.clear(); trans_plot.clear()
    angles_plot.plot(sim_times, np.rad2deg(phis), color="c", label="$\phi$")
    angles_plot.plot(sim_times, -np.rad2deg(thetas), color="m", label="$\\theta$")
    angles_plot.plot(sim_times, np.rad2deg(psis), color="pink", label="$\psi$")
    trans_plot.plot(sim_times, pns, color="c", label="North")
    trans_plot.plot(sim_times, pes, color="m", label="East")
    trans_plot.plot(sim_times, pds, color="pink", label="Height")
    angles_plot.legend(loc="upper left"); angles_plot.grid(); angles_plot.set_ylim(-60,60)
    trans_plot.legend(loc="upper left"); trans_plot.grid(); trans_plot.set_ylim(-15,15)   

    # call simulation
    UAV_anim.update(pn, pe, pd, phi, theta, psi)

    # -------increment time-------------
    sim_time += SIM.ts_simulation

    # end simulation
    if key.is_pressed("q"): break
