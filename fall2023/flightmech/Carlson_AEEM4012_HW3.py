# Annika Carlson
# Flight Mechanics, AEEM4012
# Assignment 3 - MAV Animation with Aerodynamics

import sys
sys.path.append('/Users/annikacarlson/Documents/flightmech/library')
import numpy as np
from math import cos, sin, tan
import scipy.linalg as linalg
from library.signalGenerator import signalGenerator
import keyboard
from library.MAV_animation import MAV_animation
from library.mavDynamics import mavDynamics
from library.mavAero import mavAero
from library.wind import wind
import library.aerosonde_parameters as P
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore", category=UserWarning, module="control")

state = P.states0
MAV = mavDynamics()
Aero = mavAero()
MAV_anim = MAV_animation(limits=10, multfigs=True)
Vs = np.array([[5.],[3.],[2.]])
Wind = wind(Vs)

deltaa = 0.
deltae = 0.
deltar = 0.
deltat = 1.
Va = np.sqrt(state[3][0]**2 + state[4][0]**2 + state[5][0]**2)

pn=state[0][0]
pe=state[1][0]
pd=state[2][0]
u=state[3][0]
v=state[4][0]
w=state[5][0]
phi=state[6][0]
theta=state[7][0]
psi=state[8][0]
p=state[9][0]
q=state[10][0]
r=state[11][0]

# create subplots
throttle = plt.figure(1).add_subplot(3, 15, 1)
force_plot = MAV_anim.fig.add_subplot(5, 2, 7)
axpos = force_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.02; force_plot.set_position(axpos)
moment_plot = MAV_anim.fig.add_subplot(5, 2, 9)
axpos = moment_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.02; moment_plot.set_position(axpos)
deltas_plot = MAV_anim.fig.add_subplot(5, 2, 10)
axpos = deltas_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.02; deltas_plot.set_position(axpos)
trans_plot = MAV_anim.fig.add_subplot(5, 2, 2)
axpos = trans_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.03; trans_plot.set_position(axpos)
vel_plot = MAV_anim.fig.add_subplot(5, 2, 4)
axpos = vel_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.03; vel_plot.set_position(axpos)
angles_plot = MAV_anim.fig.add_subplot(5, 2, 6)
axpos = angles_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.03; angles_plot.set_position(axpos)
angleRs_plot = MAV_anim.fig.add_subplot(5, 2, 8)
axpos = angleRs_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.03; angleRs_plot.set_position(axpos)
sim_times = []
fxs = []; fys = []; fzs = []
ls = []; ms = []; ns = []
das = []; des = []; drs = []
pns = []; pes = []; pds = []
us = []; vs = []; ws = []
phis = []; thetas = []; psis = []
ps = []; qs = []; rs = []

print("Press Q to exit simulation")
sim_time = P.start_time
t_next_plot = sim_time + P.ts_plotting
while sim_time < P.end_time:

    # # keyboard inputs
    # if keyboard.is_pressed("down arrow"): deltae -= np.deg2rad(1)
    # if keyboard.is_pressed("up arrow"): deltae += np.deg2rad(1)
    # if keyboard.is_pressed("right arrow"): deltaa += np.deg2rad(0.5); deltar -= np.deg2rad(0.25)
    # if keyboard.is_pressed("left arrow"): deltaa -= np.deg2rad(0.5); deltar += np.deg2rad(0.25)
    # if keyboard.is_pressed("space"): deltae = 0.; deltaa = 0.; deltar = 0.
    # if keyboard.is_pressed("x"):
    #     if deltat < 1: deltat += 0.1
    # if keyboard.is_pressed("z"):
    #     if deltat > 0: deltat -= 0.1

    if sim_time <= 0.5:
        deltat = 1.0
    elif sim_time <= 1:
        deltat = 0.6
    elif sim_time <= 2:
        deltae -= np.deg2rad(3)
    elif sim_time <= 3:
        deltae += np.deg2rad(3)
    elif sim_time <= 3.5:
        deltae = 0.
    elif sim_time <= 4:
        deltaa += np.deg2rad(0.25)
        deltar -= np.deg2rad(0.01)
    elif sim_time <= 5:
        deltaa -= np.deg2rad(0.25)
        deltar += np.deg2rad(0.01)
    elif sim_time <= 6:
        deltaa -= np.deg2rad(0.25)
        deltar += np.deg2rad(0.01)
    elif sim_time <= 7:
        deltaa += np.deg2rad(0.25)
        deltar -= np.deg2rad(0.01)
    elif sim_time <= 7.5:
        deltat = 1.0
    elif sim_time <= 9:
        deltae -= np.deg2rad(3)

    Va, alpha, beta = Wind.windout(state, Va, sim_time)
    fx, fy, fz = Aero.forces(state, alpha, beta, deltaa, deltae, deltar, deltat, Va)
    l, m, n, = Aero.moments(state, alpha, beta, deltaa, deltae, deltar, deltat, Va)
    state = MAV.update(fx, fy, fz, l, m, n)
    pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = state.flatten()
    MAV_anim.update(pn, pe, pd, phi, theta, psi)

    # make lists for plotting
    sim_times.append(sim_time)

    fxs.append(fx); fys.append(fy); fzs.append(fz)
    ls.append(l); ms.append(m); ns.append(n)
    das.append(deltaa); des.append(deltae); drs.append(deltar)
    pns.append(pn); pes.append(pe); pds.append(pd)
    us.append(u); vs.append(v); ws.append(w)
    phis.append(phi); thetas.append(theta); psis.append(psi)
    ps.append(p); qs.append(q); rs.append(r)

    # plotting
    if sim_time > t_next_plot:
        throttle.clear(); throttle.bar(0, deltat); throttle.set_ylim(0, 1); throttle.set_xticklabels(""); throttle.set_xticks([])
        t_next_plot = sim_time + P.ts_plotting 

    # create plots
    force_plot.clear(); moment_plot.clear()
    deltas_plot.clear()
    trans_plot.clear(); vel_plot.clear(); angles_plot.clear(); angleRs_plot.clear()

    force_plot.plot(sim_times, fxs, color='c', label='fx')
    force_plot.plot(sim_times, fys, color='m', label='fy')
    force_plot.plot(sim_times, fzs, color='pink', label='fz')
    moment_plot.plot(sim_times, ls, color='c', label='l')
    moment_plot.plot(sim_times, ms, color='m', label='m')
    moment_plot.plot(sim_times, ns, color='pink', label='n')

    deltas_plot.plot(sim_times, das, color='c', label='$\delta_a$')
    deltas_plot.plot(sim_times, des, color='m', label='$\delta_e$')
    deltas_plot.plot(sim_times, drs, color='pink', label='$\delta_r$')

    trans_plot.plot(sim_times, pns, color='c', label='North')
    trans_plot.plot(sim_times, pes, color='m', label='East')
    trans_plot.plot(sim_times, pds, color='pink', label='Height')
    vel_plot.plot(sim_times, us, color='c', label='u')
    vel_plot.plot(sim_times, vs, color='m', label='v')
    vel_plot.plot(sim_times, ws, color='pink', label='w')
    angles_plot.plot(sim_times, np.rad2deg(phis), color='c', label='$\phi$')
    angles_plot.plot(sim_times, -np.rad2deg(thetas), color='m', label='$\\theta$')
    angles_plot.plot(sim_times, np.rad2deg(psis), color='pink', label='$\psi$')
    angleRs_plot.plot(sim_times, ps, color='c', label='p')
    angleRs_plot.plot(sim_times, qs, color='m', label='q')
    angleRs_plot.plot(sim_times, rs, color='pink', label='r')

    # Set plot titles and labels
    throttle.set_title('Throttle')
    force_plot.legend(loc='upper left'); force_plot.grid(); force_plot.set_title('Forces (N/m)')
    moment_plot.legend(loc='upper left'); moment_plot.grid(); moment_plot.set_title('Moments (Nm)')
    deltas_plot.legend(loc='upper left'); deltas_plot.grid(); deltas_plot.set_title('Surface Deflections (deg)')
    trans_plot.legend(loc="upper left"); trans_plot.grid(); trans_plot.set_title('Position (m)')
    vel_plot.legend(loc="upper left"); vel_plot.grid(); vel_plot.set_title('Velocity (m/s)')
    angles_plot.legend(loc="upper left"); angles_plot.grid(); angles_plot.set_title('Angle (deg)')
    angleRs_plot.legend(loc="upper left"); angleRs_plot.grid(); angleRs_plot.set_title('Angle Rate (deg/s)')

    # -------increment time-------------
    sim_time += P.ts_simulation

    # end simulation
    if keyboard.is_pressed("q"): break
