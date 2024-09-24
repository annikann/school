import numpy as np
import matplotlib.pyplot as plt

# Define Variables

a = -0.00033617
B = -0.000064432
R = 200
u_c = -102.76
v_c = 95.93

s1 = np.linspace(0, (3/2)*np.pi*R, 100)
s2 = np.linspace(0, R, 100)
s3 = np.linspace(0, np.pi*R, 100)
s4 = np.linspace(0, (np.pi/2)*R, 100)

def f12(s): 
    return -u_c*s - (R**2)*np.sin(s/R)
def f23(s): 
    return f12(s1[-1]) - u_c*s
def f34(s): 
    return f23(s2[-1]) - u_c*s + 4*(R**2)*np.cos(s/(2*R))-4*R**2
def f45(s): 
    return f34(s3[-1]) - 3*R*s/2 - u_c*s - (R**2)/4*np.sin(2*s/R)
def g12(s): 
    return -v_c*s + (R**2)*(np.cos(s/R) - 1)
def g23(s): 
    return g12(s1[-1]) + R*s - v_c*s + (s**2)/2
def g34(s): 
    return g23(s2[-1]) -v_c*s + 4*(R**2)*np.sin(s/(2*R))
def g45(s): 
    return g34(s3[-1]) -v_c*s - (R**2)/4*(np.cos(2*s/R) - 1)

q12 = a*f12(s1) + B*g12(s1)
q23 = a*f23(s2) + B*g23(s2)
q34 = a*f34(s3) + B*g34(s3)
q45 = a*f45(s4) + B*g45(s4)

print(q12[-1], q23[-1], q34[-1])

max = np.max(q12)
min = np.min(q34)
print(max, min)

plt.figure()
plt.plot(s1, q12, label="q_s1", linewidth=2)
plt.plot(s1[-1] + s2, q23, label="q_s2", linewidth=2)
plt.plot(s1[-1] + s2[-1] + s3, q34, label="q_s3", linewidth=2)
plt.plot(s1[-1] + s2[-1] + s3[-1] + s4, q45, label="q_s4", linewidth=2)

plt.ylabel("Shear Flow $q(s)$ (N/mm)"); plt.xlabel("Arc-Length Coordinate $s$ (mm)")
plt.grid(); plt.legend(loc="lower left")

plt.show()

# print(q_s1)
