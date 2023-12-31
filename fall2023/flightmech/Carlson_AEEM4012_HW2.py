# Annika Carlson
# Flight Mechanics, AEEM4012
# Assignment 2 - MAV Animation with Dynamics

import sys
sys.path.append('/Users/annikacarlson/Documents/flightmech/library')
import numpy as np
from math import cos, sin, tan
import scipy.linalg as linalg
from library.signalGenerator import signalGenerator
import keyboard as key
from library.MAV_animation import MAV_animation
from library.mavDynamics import mavDynamics
import library.aerosonde_parameters as P

state = P.states0
MAV = mavDynamics()
MAV_anim = MAV_animation(state, scale=5, multfigs=True)

# initialize the simulation and signal generator
sim_time = P.start_time
forces = signalGenerator(1000., 1)
moments = signalGenerator(15., 1)

print("Press Q to exit simulation")
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
force_plot = MAV_anim.fig.add_subplot(4, 2, 5)
axpos = force_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.02; force_plot.set_position(axpos)
moment_plot = MAV_anim.fig.add_subplot(4, 2, 7)
axpos = moment_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.02; moment_plot.set_position(axpos)
trans_plot = MAV_anim.fig.add_subplot(4, 2, 2)
axpos = trans_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.03; trans_plot.set_position(axpos)
vel_plot = MAV_anim.fig.add_subplot(4, 2, 4)
axpos = vel_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.03; vel_plot.set_position(axpos)
angles_plot = MAV_anim.fig.add_subplot(4, 2, 6)
axpos = angles_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.03; angles_plot.set_position(axpos)
angleRs_plot = MAV_anim.fig.add_subplot(4, 2, 8)
axpos = angleRs_plot.get_position(); axpos.x0 += 0.0; axpos.y0 += 0.03; angleRs_plot.set_position(axpos)
sim_times = []
fxs = []; fys = []; fzs = []
ls = []; ms = []; ns = []
pns = []; pes = []; pds = []
us = []; vs = []; ws = []
phis = []; thetas = []; psis = []
ps = []; qs = []; rs = []

sim_time = P.start_time
while sim_time < P.end_time:
    # forces
    if sim_time <= 1:
        fx = 0.
        fy = 0.
        fz = 0.
        l = 0.
        m = 0.
        n = 0.
    elif sim_time <= 2:
        fx = forces.sin(sim_time)
        fy = 0.
        fz = 0.
        l = 0.
        m = 0.
        n = 0.
    elif sim_time <= 3:
        fx = 0.
        fy = forces.sin(sim_time)
        fz = 0.
        l = 0.
        m = 0.
        n = 0.
    elif sim_time <= 4:
        fx = 0.
        fy = 0.
        fz = forces.sin(sim_time)
        l = 0.
        m = 0.
        n = 0.
    # moments
    elif sim_time <= 5:
        fx = 0.
        fy = 0.
        fz = 0.
        l = moments.sin(sim_time)
        m = 0.
        n = 0.
    elif sim_time <= 6:
        fx = 0.
        fy = 0.
        fz = 0.
        l = 0.
        m = moments.sin(sim_time)
        n = 0.
    elif sim_time <= 7:
        fx = 0.
        fy = 0.
        fz = 0.
        l = 0.
        m = 0.
        n = moments.sin(sim_time)

    # call simulation
    y = MAV.update(fx, fy, fz, l, m, n)
    MAV_anim.update(y[0][0], y[1][0], y[2][0], y[6][0], y[7][0], y[8][0])

    # make lists for plotting
    sim_times.append(sim_time)

    fxs.append(fx); fys.append(fy), fzs.append(fz)
    ls.append(l); ms.append(m), ns.append(n)
    pns.append(y[0][0]); pes.append(y[1][0]); pds.append(y[2][0])
    us.append(y[3][0]); vs.append(y[4][0]); ws.append(y[5][0])
    phis.append(y[6][0]); thetas.append(y[7][0]); psis.append(y[8][0])
    ps.append(y[9][0]); qs.append(y[10][0]); rs.append(y[11][0])

    # create plots
    force_plot.clear(); moment_plot.clear()
    trans_plot.clear(); vel_plot.clear(); angles_plot.clear(); angleRs_plot.clear()

    force_plot.plot(sim_times, fxs, color='c', label='fx')
    force_plot.plot(sim_times, fys, color='m', label='fy')
    force_plot.plot(sim_times, fzs, color='pink', label='fz')
    moment_plot.plot(sim_times, ls, color='c', label='l')
    moment_plot.plot(sim_times, ms, color='m', label='m')
    moment_plot.plot(sim_times, ns, color='pink', label='n')

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

    # SET PLOT TITLES AND LABELS AS WELL
    force_plot.legend(loc='upper left'); force_plot.grid(); force_plot.set_title('Forces (N/m)')
    moment_plot.legend(loc='upper left'); moment_plot.grid(); moment_plot.set_title('Moments (Nm)')
    trans_plot.legend(loc="upper left"); trans_plot.grid(); trans_plot.set_title('Position (m)')
    vel_plot.legend(loc="upper left"); vel_plot.grid(); vel_plot.set_title('Velocity (m/s)')
    angles_plot.legend(loc="upper left"); angles_plot.grid(); angles_plot.set_title('Angle (deg)')
    angleRs_plot.legend(loc="upper left"); angleRs_plot.grid(); angleRs_plot.set_title('Angle Rate (deg/s)')

    # -------increment time-------------
    sim_time += P.ts_simulation

    # end simulation
    if key.is_pressed("q"): break
