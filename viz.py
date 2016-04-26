import sys
import matplotlib.pyplot as plt
import numpy as np

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

d = X[1]-X[0]

# Find midpoints of cells
Xs = []; Ys = []
for ii in range(n-1):
    for jj in range(m-1):
        Xs.append( np.mean([Xm[ii][jj],Xm[ii+1][jj],Xm[ii][jj+1],Xm[ii+1][jj+1]]))
        Ys.append( np.mean([Ym[ii][jj],Ym[ii+1][jj],Ym[ii][jj+1],Ym[ii+1][jj+1]]))

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

# Post process
U = [W2[i]/W1[i] for i in range(len(W1))]
V = [W3[i]/W1[i] for i in range(len(W1))]

# Plot Gridpoints
# fig = plt.figure()
# plt.plot(X,Y,'.k')
# plt.xlim([min(X)-d,max(X)+d])
# plt.ylim([min(Y)-d,max(Y)+d])
# plt.show()

# Plot Cells
fig = plt.figure()
# Interior grid lines
for ii in range(n-1):
    for jj in range(m-1):
        plt.plot([Xm[ii][jj],Xm[ii+1][jj]],[Ym[ii][jj],Ym[ii+1][jj]],'k-')
        plt.plot([Xm[ii][jj],Xm[ii][jj+1]],[Ym[ii][jj],Ym[ii][jj+1]],'k-')
# Right and bottom grid line
plt.plot([Xm[0][m-1],Xm[n-1][m-1]],[Ym[0][m-1],Ym[n-1][m-1]],'k-')
plt.plot([Xm[n-1][0],Xm[n-1][m-1]],[Ym[n-1][0],Ym[n-1][m-1]],'k-')
# Axis limits
plt.xlim([min(X)-d,max(X)+d])
plt.ylim([min(Y)-d,max(Y)+d])
plt.show()

# Plot velocity quiver
# fig = plt.figure()
# plt.quiver(Xs,Ys,U,V)
# plt.show()