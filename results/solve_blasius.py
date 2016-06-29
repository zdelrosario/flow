import numpy as np
from scipy.integrate import ode
from scipy.optimize import bisect
import matplotlib.pyplot as plt
from math import sqrt

def blasius_solution(T):

    # Define Blasius equation as system of 1st order ODE
    def F(t,y):
        return [y[1],y[2],-0.5*y[0]*y[2]]
    # Set solver parameters
    backend = 'zvode'
    method = 'bdf'
    # Define shooting function
    def shoot(alpha):
        y0 = [0,0,alpha]
        Y = [y0]
        # Define solver object
        solver = ode(F).set_integrator(backend,method=method)
        solver.set_initial_value(y0,T[0])
        # Run the solver
        for i in range(1,len(T)):
            Y.append(solver.integrate(T[i]))
            if solver.successful() != True:
                raise ValueError("Solver unsuccessful, exiting...")
        # Return solution
        return np.real(np.array(Y))
    # Run rootfinding
    alpha = bisect(lambda alp: shoot(alp)[-1,1]-1, 2e-1, 5e-1)
    return shoot(alpha)

def eta_fcn(x,y,nu,U):
    return y * sqrt( U/2./nu/x )

def u_fcn(fp,U):
    return U*fp

if __name__ == "__main__":

    T = np.linspace(0,12,int(1e4))
    Y = blasius_solution(T)

    fig = plt.figure()
    plt.plot(T,Y[:,1])
    plt.show()