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

# Plot
fig = plt.figure()
plt.plot(X,Y,'.k')
plt.xlim([min(X)-1,max(X)+1])
plt.ylim([min(Y)-1,max(Y)+1])
plt.show()