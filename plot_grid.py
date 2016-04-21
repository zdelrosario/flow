import matplotlib.pyplot as plt

# Import grid
f = open('grid.dat','r')
X = []; Y = []
for line in f:
    pairs = line.split(";") # Split into x,y pairs
    for pair in pairs:
        val = pair.split(",")
        X.append(float(val[0]))
        Y.append(float(val[1]))

d = X[1]-X[0]

# Plot
fig = plt.figure()
plt.plot(X,Y,'.k')
plt.xlim([min(X)-d,max(X)+d])
plt.ylim([min(Y)-d,max(Y)+d])
plt.show()