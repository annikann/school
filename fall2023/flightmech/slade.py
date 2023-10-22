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
from library.sladeAnim import MAV_animation
from library.mavDynamics import mavDynamics
from library.compute_trim import ComputeTrim
from library.mavAero import mavAero
from library.wind import wind
from library.sladegains import ComputeGains
import library.aerosonde_parameters as P
import warnings; warnings.filterwarnings("ignore", category=UserWarning, module="Control")

state = P.states0
MAV_anim = MAV_animation(limits=10, scale=0.25, multfigs=False)
MAV = mavDynamics()
Trim = ComputeTrim()
Aero = mavAero()
Vs = np.array([[0.],[0.],[0.]])
Wind = wind(Vs)
gains = ComputeGains()

# initialize the simulation and signal generator
sim_time = P.start_time

# create subplots
throttle = plt.figure(1).add_subplot(1, 20, 1)
axpos = throttle.get_position(); axpos.x0 -= 0.1; axpos.x1 -= 0.1; throttle.set_position(axpos)
force_plot = MAV_anim.fig.add_subplot(334)
axpos = force_plot.get_position(); axpos.x0 -= 0.015; axpos.x1 -= 0.05; force_plot.set_position(axpos)
moment_plot = MAV_anim.fig.add_subplot(337)
axpos = moment_plot.get_position(); axpos.x0 -= 0.015; axpos.x1 -= 0.05; moment_plot.set_position(axpos)
deltas_plot = MAV_anim.fig.add_subplot(331)
axpos = deltas_plot.get_position(); axpos.x0 -= 0.015; axpos.x1 -= 0.05; deltas_plot.set_position(axpos)
trans_plot = MAV_anim.fig.add_subplot(433)
axpos = trans_plot.get_position(); axpos.x0 += 0.125; axpos.x1 += 0.09; trans_plot.set_position(axpos)
vel_plot = MAV_anim.fig.add_subplot(436)
axpos = vel_plot.get_position(); axpos.x0 += 0.125; axpos.x1 += 0.09; vel_plot.set_position(axpos)
angles_plot = MAV_anim.fig.add_subplot(439)
axpos = angles_plot.get_position(); axpos.x0 += 0.125; axpos.x1 += 0.09; angles_plot.set_position(axpos)
angleRs_plot = MAV_anim.fig.add_subplot(4, 3, 12)
axpos = angleRs_plot.get_position(); axpos.x0 += 0.125; axpos.x1 += 0.09; angleRs_plot.set_position(axpos)

# set up data arrays
buffer = int(1/P.ts_plotting)
sim_times = np.zeros(buffer)
deltaas = np.zeros(buffer)
deltaes = np.zeros(buffer)
deltars = np.zeros(buffer)
fxs = np.zeros(buffer); fys = np.zeros(buffer); fzs = np.zeros(buffer)
ls = np.zeros(buffer); ms = np.zeros(buffer); ns = np.zeros(buffer)
pns = np.zeros(buffer); pes = np.zeros(buffer); pds = np.zeros(buffer)
us = np.zeros(buffer); vs = np.zeros(buffer); ws = np.zeros(buffer)
phis = np.zeros(buffer); thetas = np.zeros(buffer); psis = np.zeros(buffer)
ps = np.zeros(buffer); qs = np.zeros(buffer); rs = np.zeros(buffer)

Va = 35.
Y = np.deg2rad(10)
R = 100.
alpha = 0
beta = 0

x_trim, u_trim = Trim.compute_trim(Va, Y, R)
deltae, deltat, deltaa, deltar = u_trim.flatten()
print("\n--- Trim Conditions ---")
print(f"Elevator: {np.rad2deg(deltae):.2f} deg")
print(f"Throttle: {deltat*100:.2f} %")
print(f"Aileron:  {np.rad2deg(deltaa):.2f} deg")
print(f"Rudder:   {np.rad2deg(deltar):.2f} deg")

# Compute transfer functions and state spaces
T_phi_delta_a, T_chi_phi, T_theta_delta_e, T_h_theta, T_h_Va, T_Va_delta_t, T_Va_theta, T_Va_theta, T_beta_delta_r = gains.compute_tfs(x_trim, u_trim)
Alat, Blat, elatvalue, elongvalue = gains.statespace(x_trim, u_trim)

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

print("\nPress Q to exit simulation")
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
        sim_time += P.ts_simulation

    # make lists for plotting
    sim_times = np.concatenate((sim_times[1:], [sim_time]))
    deltaas = np.concatenate((deltaas[1:], [np.degrees(deltaa)]))
    deltaes = np.concatenate((deltaes[1:], [np.degrees(deltae)]))
    deltars = np.concatenate((deltars[1:], [np.degrees(deltar)]))
    fxs = np.concatenate((fxs[1:], [fx]))
    fys = np.concatenate((fys[1:], [fy]))
    fzs = np.concatenate((fzs[1:], [fz]))
    ls = np.concatenate((ls[1:], [l]))
    ms = np.concatenate((ms[1:], [m]))
    ns = np.concatenate((ns[1:], [n]))
    pns = np.concatenate((pns[1:], [y[0][0]]))
    pes = np.concatenate((pes[1:], [y[1][0]]))
    pds = np.concatenate((pds[1:], [-y[2][0]]))
    us = np.concatenate((us[1:], [y[3][0]]))
    vs = np.concatenate((vs[1:], [y[4][0]]))
    ws = np.concatenate((ws[1:], [y[5][0]]))
    phis = np.concatenate((phis[1:], [np.degrees(y[6][0])]))
    thetas = np.concatenate((thetas[1:], [np.degrees(y[7][0])]))
    psis = np.concatenate((psis[1:], [np.degrees(y[8][0])]))
    ps = np.concatenate((ps[1:], [np.degrees(y[9][0])]))
    qs = np.concatenate((qs[1:], [np.degrees(y[10][0])]))
    rs = np.concatenate((rs[1:], [np.degrees(y[11][0])]))

    # create plots
    throttle.clear(); throttle.bar(0, deltat); throttle.set_ylim(0, 1); throttle.set_xticklabels(""); throttle.set_xticks([])

    force_plot.clear(); moment_plot.clear(); deltas_plot.clear()
    trans_plot.clear(); vel_plot.clear(); angles_plot.clear(); angleRs_plot.clear()

    force_plot.plot(sim_times, fxs, label='fx')
    force_plot.plot(sim_times, fys, label='fy')
    force_plot.plot(sim_times, fzs, label='fz')
    moment_plot.plot(sim_times, ls, label='l')
    moment_plot.plot(sim_times, ms, label='m')
    moment_plot.plot(sim_times, ns, label='n')
    deltas_plot.plot(sim_times, deltaas, label='a')
    deltas_plot.plot(sim_times, deltaes, label='e')
    deltas_plot.plot(sim_times, deltars, label='r')

    trans_plot.plot(sim_times, pns, label='North')
    trans_plot.plot(sim_times, pes, label='East')
    trans_plot.plot(sim_times, pds, label='Height')
    vel_plot.plot(sim_times, us, label='u')
    vel_plot.plot(sim_times, vs, label='v')
    vel_plot.plot(sim_times, ws, label='w')
    angles_plot.plot(sim_times, phis, label='$\phi$')
    angles_plot.plot(sim_times, thetas, label='$\\theta$')
    angles_plot.plot(sim_times, psis, label='$\psi$')
    angleRs_plot.plot(sim_times, ps, label='p')
    angleRs_plot.plot(sim_times, qs, label='q')
    angleRs_plot.plot(sim_times, rs, label='r')

    # SET PLOT TITLES AND LABELS AS WELL
    throttle.set_title('Throttle')
    force_plot.legend(loc='upper right'); force_plot.grid(); force_plot.set_ylabel('Force (N)')
    force_plot.set_ylim((-1000, 1000))
    moment_plot.legend(loc='upper right'); moment_plot.grid(); moment_plot.set_ylabel('Moment (Nm)')
    moment_plot.set_ylim((-40, 40))
    trans_plot.legend(loc="upper right"); trans_plot.grid(); trans_plot.set_ylabel('Position (m)')
    vel_plot.legend(loc="upper right"); vel_plot.grid(); vel_plot.set_ylabel('Velocity (m/s)')
    vel_plot.set_ylim((-25, 100))
    angles_plot.legend(loc="upper right"); angles_plot.grid(); angles_plot.set_ylabel('Angle (deg)'); angles_plot.set_ylim((-180,180))
    angleRs_plot.legend(loc="upper right"); angleRs_plot.grid(); angleRs_plot.set_ylabel('Angle Rate (deg/s)')
    angleRs_plot.set_ylim((-360, 360))
    deltas_plot.legend(loc="upper right"); deltas_plot.grid(); deltas_plot.set_ylabel('Deflections (deg)')
    deltas_plot.set_ylim((-45, 45))

    plt.pause(0.1)

    # end simulation
    if key.is_pressed("q"): break
