# Quick plotting for HW 5 part d

import numpy as np
import matplotlib.pyplot as plt

def t_0(n): return 1/(1 + n/(4*np.pi))
def t_R(n): return np.sqrt(1 + n/(12*np.pi))/(1 + n/(4*np.pi))

n = np.linspace(1, 10000, 1000)

t_small = t_0(n)
t_large = t_R(n)

plt.plot(n, t_small, label=r't $\rightarrow$ 0')
plt.plot(n, t_large, label=r't $\rightarrow$ R')
plt.ylabel('A\'/A'); plt.xlabel('Number of fins (n)')
# plt.title()
plt.grid(); plt.legend(loc='upper right')
plt.show()