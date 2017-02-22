# ##############################################
# # Kriek van der Meulen
# # Stellar Structures
# # 2017
# ##############################################
#
from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# Function for plotting analytical solution
def anaplot(analist,color):
    ax1.plot(xlist,analist,linestyle="dashed",linewidth ="3",color = color,label="analytic n = "+str(n))
    anafinalrho = np.array(analist)**n
    ax2.plot(xlist,anafinalrho,linestyle="dashed",linewidth ="3",color = color,label="analytic n = "+str(n))

# Function for ode integration
def Ode(solreal,ksireal,n,dksi,begin,end,initialtheta,initialphi):

    # define function used to integrate
    def Stellar(y,ksi,n,a):
        theta,phi = y
        dydt = phi/(ksi**2),-a*(ksi**2 * (np.power(theta,n)))
        return dydt

    # Starting values are defined here
    y0 = [initialtheta,initialphi]
    a=1
    ksi = np.arange(begin,end,dksi)

    # perform function
    sol = odeint(Stellar,y0,ksi,args=(n,a))

    # fetch new initial values to pass to the function
    initialtheta = sol[:,0][-1]
    initialphi = sol[:,1][-1]

    # concatenate the arrays
    ksireal = np.append(ksireal,ksi)
    solreal = np.append(solreal,sol[:,0])

    # Check whether we found a solution and return

    for i in range(0,len(solreal)):
        if solreal[i] <=0 or isnan(solreal[i]):
            return ksireal[:i],solreal[:i]

    # call Ode on the stack with new initial values and return when solution is found
    begin = end
    end +=end
    return Ode(solreal,ksireal,n,dksi,begin,end,initialtheta,initialphi)


# define figure
fig,(ax1,ax2) = plt.subplots(2,1,figsize=[20,20])
ax1.set_ylim([-0.1,1.1])
ax2.set_ylim([-0.1,1.1])
ax1.set_xlabel(r'$\xi$',fontsize=18)
ax1.set_ylabel(r'$\theta$',fontsize=18,rotation=360)
ax2.set_xlabel(r'$\xi$',fontsize=18)
ax2.set_ylabel(r'$\frac{\rho}{\rho_c}$',fontsize=22,rotation=360)

# define nvalues to loop over and check whether using arrayspace
nvalues = [0,1,3./2,2,3,4]
# nvalues = [0]
arrayspace = True
for n in nvalues:
    # define starting values
    initialtheta = 1.0
    initialphi = 0.
    dksi = 0.01
    end = 10
    begin = dksi
    solreal = []
    ksireal = []

    # perform calculations
    xlist,ylist = Ode(solreal,ksireal,n,dksi,begin,end,initialtheta,initialphi)

    if n == 0:
        analist = (1-(1./6)*xlist**2)
        color = "red"
        anaplot(analist,color)

    if n == 1:
        analist = np.exp((-1./6)*xlist**2)
        color = "blue"
        anaplot(analist,color)

    # plot everything!
    ax1.plot(xlist,ylist,label="n = "+str(n) )
    ax1.axhline(y=0,linestyle="dotted")
    finalrho = np.array(ylist)**n
    ax2.plot(xlist,finalrho,label="n = "+str(n))
    ax2.axhline(y=0,linestyle="dotted")

ax1.legend()
ax2.legend()
plt.show()
