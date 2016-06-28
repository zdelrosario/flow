# Constants
D   = 1.0       # Length scale, [m]
nu  = 1.544e-5  # Kinematic viscosity of air, [m^2/s]
gam = 1.4       # Adiabatic coefficient [-]
rho = 1.1462    # Density of air, [kg/m^3]
# Chosen dimensionless parameters
Re = 5e7
# M  = 0.9
M = 2.0
# Calculated parameters
U_inf = Re * nu / D
P_inf = rho/gam * (U_inf/M)**2
# State values
w0 = rho
w1 = rho*U_inf
w2 = 0
w3 = P_inf/(gam-1) + 0.5*rho*U_inf**2

e_inf  = w3 / rho

print("Parameters:")
print("U_inf = {:.4e}".format(U_inf))
print("P_inf = {:.4e}".format(P_inf))
print("e_inf = {:.4e}".format(e_inf))
print("State Values:")
print("w0 = {:.4e}".format(w0))
print("w1 = {:.4e}".format(w1))
print("w2 = {:.4e}".format(w2))
print("w3 = {:.4e}".format(w3))