import sys
import numpy as np
import matplotlib
# choose the backend to handle window geometry
matplotlib.use("Qt4Agg")
# Import pyplot
import matplotlib.pyplot as plt

from solve_blasius import blasius_solution, eta_fcn, u_fcn

from math import sqrt

# Plot settings
offset = [(0,0),(700,0),(1400,0)] # Plot window locations
d = 0.01
e = 1e-5

U_inf = 68.93 # m/s
nu = 15.11e-6 # m^2/s

# Command line argument form
# if len(sys.argv) < 2:
#     print('Usage:')
#     print('    python {} [grid file] [solution file]'.format(sys.argv[0]))
#     exit()

# grid_file = sys.argv[1]
# sol_file  = sys.argv[2]

# DEBUG -- fixed input arguments
grid_file = "solution.grid.dat"
sol_file  = "solution.val.dat"

# grid_file = "restart.grid.dat"
# sol_file  = "restart.val.dat"

##################################################
# Post-process data
##################################################
# Import grid
f = open(grid_file,'r')
# First line is grid dimensions
pair = f.readline(); val = pair.split(",")
n = int(val[0]); m = int(val[1])
X = []; Y = []
for line in f:
    pairs = line.split(";") # Split into x,y pairs
    for pair in pairs:
        val = pair.split(",")
        X.append(float(val[0]))
        Y.append(float(val[1]))
# Reshape arrays
Xm = np.reshape(X,(n,m))
Ym = np.reshape(Y,(n,m))

# Find midpoints of cells
Xs = []; Ys = []
for ii in range(n-1):
    for jj in range(m-1):
        Xs.append( np.mean([Xm[ii][jj],Xm[ii+1][jj],Xm[ii][jj+1],Xm[ii+1][jj+1]]))
        Ys.append( np.mean([Ym[ii][jj],Ym[ii+1][jj],Ym[ii][jj+1],Ym[ii+1][jj+1]]))
Xs = np.reshape(np.array(Xs),(-1,m-1))
Ys = np.reshape(np.array(Ys),(-1,m-1))

# Import solution
f = open(sol_file,'r')
W1 = []; W2 = []; W3 = []; W4 = []
for line in f:
    state = line.split(";") # Split into solution vectors
    for element in state:
        val = element.split(",")
        W1.append(float(val[0]))
        W2.append(float(val[1]))
        W3.append(float(val[2]))
        W4.append(float(val[3]))
# Velocities
U = np.reshape(np.array([W2[i]/W1[i] for i in range(len(W1))]),(-1,m-1))
V = np.reshape(np.array([W3[i]/W1[i] for i in range(len(W1))]),(-1,m-1))

W1 = np.reshape(np.array(W1),(-1,m-1))
W2 = np.reshape(np.array(W2),(-1,m-1))
W3 = np.reshape(np.array(W3),(-1,m-1))
W4 = np.reshape(np.array(W4),(-1,m-1))

# Non-dimensionalized results
ind = 20    # Horizontal station index
U_c = U[:,ind]
X_c = Xs[:,ind]
Y_c = Ys[:,ind]
Eta_c = np.array([eta_fcn(X_c[i],Y_c[i],nu,U_inf) for i in range(len(U_c))])
Utl_c = U_c / U_inf

##################################################
# Blasius solution
##################################################
Eta_b = np.linspace(0,12,int(1e4))
F = blasius_solution(Eta_b)
f = F[:,0]; fp = F[:,1]
Utl_b = fp

##################################################
# Plotting
##################################################

# Global View
fig = plt.figure()
# Velocity contour
cs = plt.contourf(Xs,Ys,W2)
plt.colorbar(cs)
# Axis limits
plt.xlim([min(X)-d,max(X)+d])
plt.ylim([min(Y)-d,max(Y)+d])
# Set plot location on screen
manager = plt.get_current_fig_manager()
x,y,dx,dy = manager.window.geometry().getRect()
manager.window.setGeometry(offset[0][0],offset[0][1],dx,dy)

# Velocity profile
fig = plt.figure()
plt.plot(Utl_c,Eta_c,'*')
plt.plot(Utl_b,Eta_b)
plt.ylim([Eta_b[0],Eta_b[-1]])
# Set plot location on screen
manager = plt.get_current_fig_manager()
x,y,dx,dy = manager.window.geometry().getRect()
manager.window.setGeometry(offset[1][0],offset[1][1],dx,dy)

# Show all plots
plt.show()