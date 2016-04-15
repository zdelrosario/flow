# Standard libraries
import numpy as np
from numpy import exp
import matplotlib.pyplot as plt
from scipy.optimize import newton

# Functions
H = 0.5
Nt  = 34
Nbl = 22
yfm = 5.4e-4

r = 0.2 / (Nbl-2)
dy_fm = r*yfm

# Fine mesh
f = lambda k: r - (exp(k/(Nbl-2))-1)/(exp(k)-1)
kfm = newton(f,2)
print(kfm)

# Coarse mesh
dy_cm = yfm*( (exp(kfm*(Nbl-1)/(Nbl-2))-1)/(exp(kfm)-1) -1)
f2= lambda k: dy_cm - (H-yfm)*(exp(k/(Nt-Nbl))-1)/(exp(k)-1)
kcm = newton(f2,8)
print(kcm)

# Plot
K = np.linspace(1e-3,10)

plt.figure()
plt.plot(K,f(K),'b')
plt.plot(K,f2(K),'r')

plt.show()