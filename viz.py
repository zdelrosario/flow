import sys
import numpy as np

import matplotlib
# choose the backend to handle window geometry
matplotlib.use("Qt4Agg")
# Import pyplot
import matplotlib.pyplot as plt

offset = [(1400,0),(0,500),(700,500),(1400,500)] # Plot window locations

d = 0.01
e = 1e-5

# Command line argument form
# if len(sys.argv) < 2:
#     print('Usage:')
#     print('    python {} [grid file] [solution file]'.format(sys.argv[0]))
#     exit()

# grid_file = sys.argv[1]
# sol_file  = sys.argv[2]

# DEBUG -- solver final output
grid_file = "solution.grid.dat"
sol_file  = "solution.val.dat"
# DEBUG -- solver restart file
# grid_file = "restart.grid.dat"
# sol_file  = "restart.val.dat"

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
U = [W2[i]/W1[i] for i in range(len(W1))]
V = [W3[i]/W1[i] for i in range(len(W1))]

W1 = np.reshape(np.array(W1),(-1,m-1))
W2 = np.reshape(np.array(W2),(-1,m-1))
W3 = np.reshape(np.array(W3),(-1,m-1))
W4 = np.reshape(np.array(W4),(-1,m-1))

# Plot Gridpoints
# fig = plt.figure()
# plt.plot(X,Y,'.k')
# plt.xlim([min(X)-d,max(X)+d])
# plt.ylim([min(Y)-d,max(Y)+d])
# plt.show()

# Global View
# fig = plt.figure()
# # Interior grid lines
# # for ii in range(n-1):
# #     for jj in range(m-1):
# #         plt.plot([Xm[ii][jj],Xm[ii+1][jj]],[Ym[ii][jj],Ym[ii+1][jj]],'k-')
# #         plt.plot([Xm[ii][jj],Xm[ii][jj+1]],[Ym[ii][jj],Ym[ii][jj+1]],'k-')
# # # Right and bottom grid line
# # plt.plot([Xm[0][m-1],Xm[n-1][m-1]],[Ym[0][m-1],Ym[n-1][m-1]],'k-')
# # plt.plot([Xm[n-1][0],Xm[n-1][m-1]],[Ym[n-1][0],Ym[n-1][m-1]],'k-')
# # Velocity contour
# cs = plt.contourf(Xs,Ys,W1); plt.title('Density')
# # cs = plt.contourf(Xs,Ys,W2); plt.title('Horizontal Momentum')
# plt.colorbar(cs)
# # Axis limits
# plt.xlim([min(X)-d,max(X)+d])
# plt.ylim([min(Y)-d,max(Y)+d])

# # Short View
# fig = plt.figure()
# # Interior grid lines
# for ii in range(n-1):
#     for jj in range(m-1):
#         plt.plot([Xm[ii][jj],Xm[ii+1][jj]],[Ym[ii][jj],Ym[ii+1][jj]],'k-')
#         plt.plot([Xm[ii][jj],Xm[ii][jj+1]],[Ym[ii][jj],Ym[ii][jj+1]],'k-')
# # Right and bottom grid line
# plt.plot([Xm[0][m-1],Xm[n-1][m-1]],[Ym[0][m-1],Ym[n-1][m-1]],'k-')
# plt.plot([Xm[n-1][0],Xm[n-1][m-1]],[Ym[n-1][0],Ym[n-1][m-1]],'k-')
# # Velocity contour
# cs = plt.contourf(Xs,Ys,W2,shading='flat')
# plt.colorbar(cs)
# # Axis limits
# plt.xlim([min(X)-d,max(X)+d])
# plt.ylim([-1e-5,e])

# Simple matrix plot
cs = plt.matshow( W1 ); plt.title('Density')
plt.colorbar(cs)

# cs = plt.matshow( W2 ); plt.title('Horizontal Momentum')
# plt.colorbar(cs)

##################################################
# Density
##################################################
fig = plt.figure()
# Velocity contour
cs = plt.contourf(Xs,Ys,W1)
plt.colorbar(cs)
plt.title("Density")
# Axis limits
plt.xlim([min(X)-d,max(X)+d])
plt.ylim([min(Y)-d,max(Y)+d])
# Set plot location on screen
manager = plt.get_current_fig_manager()
x,y,dx,dy = manager.window.geometry().getRect()
manager.window.setGeometry(offset[0][0],offset[0][1],dx,dy)

##################################################
# Velocity
##################################################
fig = plt.figure()
# Velocity contour
cs = plt.contourf(Xs,Ys,(W2**2+W3**2)/W1**2)
plt.colorbar(cs)
plt.title("Velocity Magnitude")
# Axis limits
plt.xlim([min(X)-d,max(X)+d])
plt.ylim([min(Y)-d,max(Y)+d])
# Set plot location on screen
manager = plt.get_current_fig_manager()
x,y,dx,dy = manager.window.geometry().getRect()
manager.window.setGeometry(offset[1][0],offset[1][1],dx,dy)

##################################################
# Velocity Quiver
##################################################
fig = plt.figure()
plt.quiver(Xs,Ys,U,V)
plt.title("Velocity Quiver")
# Set plot location on screen
manager = plt.get_current_fig_manager()
x,y,dx,dy = manager.window.geometry().getRect()
manager.window.setGeometry(offset[2][0],offset[2][1],dx,dy)

##################################################
# Pressure
##################################################
fig = plt.figure()
# Velocity contour
cs = plt.contourf(Xs,Ys,(0.4)*(W4-0.5*(W2**2+W3**2)))
plt.colorbar(cs)
plt.title("Pressure")
# Axis limits
plt.xlim([min(X)-d,max(X)+d])
plt.ylim([min(Y)-d,max(Y)+d])
# Set plot location on screen
manager = plt.get_current_fig_manager()
x,y,dx,dy = manager.window.geometry().getRect()
manager.window.setGeometry(offset[3][0],offset[3][1],dx,dy)

# Show all plots
plt.show()