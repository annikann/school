# project 2 test file

from ramp import ramp

M0 = 0.9
M2 = 0.5
A2 = 1749.209
A1 = 716.999
theta1 = 10.634
theta2 = 14.535
pi_d = 0.98

results = ramp(M0, M2, theta1, theta2, A1, A2, pi_d)
print(results)