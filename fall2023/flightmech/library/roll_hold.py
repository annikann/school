"""
Created on Mon Mar 21 22:55:46 2022

@author: raj
"""

def roll_hold(phi_c, phi, p, flag, dt):
    global roll_integrator
    global roll_differentiator
    global roll_error_d1
    
    limit1 = 0.7854
    limit2 = -0.7854
    
    kp = 0.1
    kd = 0.01
    ki = 0
    
    if flag == 1:
        roll_integrator = 0
        roll_differentiator = 0
        roll_error_d1 = 0
    
    error = phi_c - phi
    roll_integrator = roll_integrator + (dt/2)*(error + roll_error_d1)
    roll_differentiator = p
    roll_error_d1 = error
    
    u = kp*error + ki*roll_integrator + kd*roll_differentiator
    u_sat = sat(u, limit1, limit2)

    if ki != 0:
        roll_integrator = roll_integrator + dt/ki*(u_sat- u)
        
    return u_sat

def sat(inn, up_limit, low_limit):
    if inn > up_limit:
        out = up_limit
    elif inn < low_limit:
        out = low_limit
    else:
        out = inn
    return out
