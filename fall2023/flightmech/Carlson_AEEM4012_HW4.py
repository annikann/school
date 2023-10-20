# Annika Carlson
# Flight Mechanics, AEEM4012
# Assignment 4 - MAV Trim Conditions

import sys
sys.path.append('/Users/annikacarlson/Documents/flightmech/library')
import numpy as np
import matplotlib.pyplot as plt
from numpy import cos, sin, tan
import scipy.linalg as linalg
from library.signalGenerator import signalGenerator
import keyboard as key
from library.MAV_animation import MAV_animation
from library.mavDynamics import mavDynamics
from library.compute_trim import ComputeTrim
from library.mavAero import mavAero
from library.wind import wind
import library.aerosonde_parameters as P

state = P.states0
MAV_anim = MAV_animation(limits=10, scale=0.25, multfigs=True)
MAV = mavDynamics()
Trim = ComputeTrim()
Aero = mavAero()
Vs = np.array([[0.],[0.],[0.]])
Wind = wind(Vs)

# initialize the simulation and signal generator
sim_time = P.start_time

# pn = state[0][0]
# pe = state[1][0]
# pd = state[2][0]
# u = state[3][0]
# v = state[4][0]
# w = state[5][0]
# phi = state[6][0]
# theta = state[7][0]
# psi = state[8][0]
# p = state[9][0]
# q = state[10][0]
# r = state[11][0]

# create subplots
throttle = plt.figure(1).add_subplot(1, 60, 1)
axpos = throttle.get_position(); axpos.x0 -= 0.05; axpos.x1 -= 0.03; throttle.set_position(axpos)
force_plot = MAV_anim.fig.add_subplot(4, 3, 1)
axpos = force_plot.get_position(); axpos.x0 += 0.02; axpos.x1 += 0.02; axpos.y0 -= 0.05; axpos.y1 -= 0.05; force_plot.set_position(axpos)
moment_plot = MAV_anim.fig.add_subplot(4, 3, 4)
axpos = moment_plot.get_position(); axpos.x0 += 0.02; axpos.x1 += 0.02; axpos.y0 -= 0.1; axpos.y1 -= 0.1; moment_plot.set_position(axpos)
deltas_plot = MAV_anim.fig.add_subplot(4, 3, 7)
axpos = deltas_plot.get_position(); axpos.x0 += 0.02; axpos.x1 += 0.02; axpos.y0 -=0.15; axpos.y1 -=0.15; deltas_plot.set_position(axpos)
trans_plot = MAV_anim.fig.add_subplot(4, 3, 3)
axpos = trans_plot.get_position(); axpos.x0 += 0.02; axpos.x1 += 0.02; axpos.y0 += 0.03; trans_plot.set_position(axpos)
vel_plot = MAV_anim.fig.add_subplot(4, 3, 6)
axpos = vel_plot.get_position(); axpos.x0 += 0.02; axpos.x1 += 0.02; axpos.y0 += 0.03; vel_plot.set_position(axpos)
angles_plot = MAV_anim.fig.add_subplot(4, 3, 9)
axpos = angles_plot.get_position(); axpos.x0 += 0.02; axpos.x1 += 0.02; axpos.y0 += 0.03; angles_plot.set_position(axpos)
angleRs_plot = MAV_anim.fig.add_subplot(4, 3, 12)
axpos = angleRs_plot.get_position(); axpos.x0 += 0.02; axpos.x1 += 0.02; axpos.y0 += 0.03; angleRs_plot.set_position(axpos)

sim_times = []
fxs = []; fys = []; fzs = []
ls = []; ms = []; ns = []
pns = []; pes = []; pds = []
us = []; vs = []; ws = []
phis = []; thetas = []; psis = []
ps = []; qs = []; rs = []
deltaas = []; deltaes = []; deltars = []; deltats = []

Va = 35.
Y = np.deg2rad(7)
R = np.inf
# deltat = 0.5
alpha = 0
beta = 0

x_trim, u_trim = Trim.compute_trim(Va, Y, R)
deltae, deltat, deltaa, deltar = u_trim.flatten()
print("--- Trim Conditions ---")
print(f"E: {np.rad2deg(deltae):.2f} deg")
print(f"T: {deltat*100:.2f} %")
print(f"A: {np.rad2deg(deltaa):.2f} deg")
print(f"R: {np.rad2deg(deltar):.2f} deg")

pn = 0
pe = 0
pd = -100
u = x_trim.item(3)
v = x_trim.item(4)
w = x_trim.item(5)
phi = x_trim.item(6)
theta = x_trim.item(7)
psi = x_trim.item(8)
p = x_trim.item(9)
q = x_trim.item(10)
r = x_trim.item(11)

states = np.array([pn, pe, pd, u, v, w, phi, theta, psi, p, q, r])
state0 = np.array([[pn], [pe], [pd], [u], [v], [w], [phi], [theta], [psi], [p], [q], [r]])
MAV.state = np.ndarray.copy(state0)

Va_actual = np.sqrt(u**2 + v**2 + w**2)

print("Press Q to exit simulation")
sim_time = P.start_time
while sim_time < P.end_time:
    t_next_plot = sim_time + P.ts_plotting
    while sim_time < t_next_plot:
        # call simulation
        Va_actual, alpha, beta = Wind.windout(MAV.state, Va_actual, sim_time)
        fx, fy, fz = Aero.forces(MAV.state, alpha, beta, deltaa, deltae, deltar, deltat, Va_actual)
        l, m, n = Aero.moments(MAV.state, alpha, beta, deltaa, deltae, deltar, deltat, Va_actual)
        y = MAV.update(fx, fy, fz, l, m, n)

        MAV_anim.update(y[0][0], y[1][0], y[2][0], y[6][0], y[7][0], y[8][0])

        # make lists for plotting
        sim_times.append(sim_time)

        fxs.append(fx); fys.append(fy), fzs.append(fz)
        ls.append(l); ms.append(m), ns.append(n)
        pns.append(y[0][0]); pes.append(y[1][0]); pds.append(-1*y[2][0])
        us.append(y[3][0]); vs.append(y[4][0]); ws.append(y[5][0])
        phis.append(y[6][0]); thetas.append(y[7][0]); psis.append(y[8][0])
        ps.append(y[9][0]); qs.append(y[10][0]); rs.append(y[11][0])
        deltaas.append(deltaa); deltaes.append(deltae); deltars.append(deltar)

        # create plots
        throttle.clear(); throttle.bar(0, deltat); throttle.set_ylim(0, 1); throttle.set_xticklabels(""); throttle.set_xticks([])

        force_plot.clear(); moment_plot.clear(); deltas_plot.clear()
        trans_plot.clear(); vel_plot.clear(); angles_plot.clear(); angleRs_plot.clear()

        force_plot.plot(sim_times, fxs, color='c', label='fx')
        force_plot.plot(sim_times, fys, color='m', label='fy')
        force_plot.plot(sim_times, fzs, color='pink', label='fz')
        moment_plot.plot(sim_times, ls, color='c', label='l')
        moment_plot.plot(sim_times, ms, color='m', label='m')
        moment_plot.plot(sim_times, ns, color='pink', label='n')
        deltas_plot.plot(sim_times, deltaas, color='c', label='a')
        deltas_plot.plot(sim_times, deltaes, color='m', label='e')
        deltas_plot.plot(sim_times, deltars, color='pink', label='r')

        trans_plot.plot(sim_times, pns, color='c', label='North')
        trans_plot.plot(sim_times, pes, color='m', label='East')
        trans_plot.plot(sim_times, pds, color='pink', label='Height')
        vel_plot.plot(sim_times, us, color='c', label='u')
        vel_plot.plot(sim_times, vs, color='m', label='v')
        vel_plot.plot(sim_times, ws, color='pink', label='w')
        angles_plot.plot(sim_times, np.rad2deg(phis), color='c', label='$\phi$')
        angles_plot.plot(sim_times, np.rad2deg(thetas), color='m', label='$\\theta$')
        angles_plot.plot(sim_times, np.rad2deg(psis), color='pink', label='$\psi$')
        angleRs_plot.plot(sim_times, ps, color='c', label='p')
        angleRs_plot.plot(sim_times, qs, color='m', label='q')
        angleRs_plot.plot(sim_times, rs, color='pink', label='r')

        # SET PLOT TITLES AND LABELS AS WELL
        throttle.set_title('Throttle')
        force_plot.legend(loc='upper left'); force_plot.grid(); force_plot.set_title('Forces (N/m)')
        moment_plot.legend(loc='upper left'); moment_plot.grid(); moment_plot.set_title('Moments (Nm)')
        trans_plot.legend(loc="upper left"); trans_plot.grid(); trans_plot.set_title('Position (m)')
        vel_plot.legend(loc="upper left"); vel_plot.grid(); vel_plot.set_title('Velocity (m/s)')
        angles_plot.legend(loc="upper left"); angles_plot.grid(); angles_plot.set_title('Angle (deg)'); angles_plot.set_ylim(-100,100)
        angleRs_plot.legend(loc="upper left"); angleRs_plot.grid(); angleRs_plot.set_title('Angle Rate (deg/s)')
        deltas_plot.legend(loc="upper left"); deltas_plot.grid(); deltas_plot.set_title('Deflections (rads)')

    # -------increment time-------------
        sim_time += P.ts_simulation
    # plt.pause(0.0001)

    # end simulation
    if key.is_pressed("q"): break
