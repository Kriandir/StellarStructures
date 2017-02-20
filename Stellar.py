##############################################
# Kriek van der Meulen
# Stellar Structures
# 2017
##############################################

from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# theta = 1
# phi = 0
# ksi = 0
# dksi = 0.01
# n = 1
# result = []

# Function for ode integration
def Ode(n,ksi,dksi,arrayspace):
    # define function used to integrate
    def Stellar(y,ksi,n,a):
        theta,phi = y
        dydt = phi/(ksi**2),-a*(ksi**2 * (np.power(theta,n)))
        return dydt
    # initial values, theta = 1 for xi = 0 and phi = 0 for xi = 0
    y0 = [1.,0.]
    a=1

    # for doing a non-np.array integration of function stellar
    def solfunc(ksi):
        sol = odeint(Stellar,y0,ksi,args=(n,a))
        return sol
    # check if doing np.arry integration
    if arrayspace:
        ksi = np.arange(0.001,20,0.001)
        sol = odeint(Stellar,y0,ksi,args=(n,a))
        return ksi,sol[:,0]

    # perform non-np.array integration(warning a metric ton slower but as a
    #  benefit has no limit on xi)
    else:
        sollist = []
        i=0
        print "hoi"
        while True:

            sollist = solfunc(ksi)
            if sollist[i][0] <= 0 or isnan(sollist[i][0]):
                break
            c = ksi[i] + dksi
            ksi.append(c)
            i+=1
        return ksi,sollist[:,0]



#
# def Euler(n,dksi):
#     theta = 1
#     phi = 0
#     ksi = 0
#     ylist = []
#     xlist = []
#     while theta > 0:
#         ksi+=dksi
#         theta += phi/(ksi**2) * dksi
#         phi += -1*(ksi**2 * (theta**n))*dksi
#         result.append([phi,theta,ksi])
#         print ksi
#         xlist.append(ksi)
#         ylist.append(theta)
#     return xlist,ylist

# print result
# Euler(5,0.1)
# plt.plot(xlist,ylist)
# plt.show()

# define figure
fig,(ax1,ax2) = plt.subplots(2,1,figsize=[20,20])
ax1.set_xlabel(r'$\xi$',fontsize=18)
ax1.set_ylabel(r'$\theta$',fontsize=18,rotation=360)
ax2.set_xlabel(r'$\xi$',fontsize=18)
ax2.set_ylabel(r'$\frac{\rho}{\rho_c}$',fontsize=22,rotation=360)

# define nvalues to loop over and check whether using arrayspace
nvalues = [0,1,3./2,2,3,4]
arrayspace = True
for n in nvalues:
    if arrayspace:
        ksi = [0.01]
        dksi = 0.01
    # perform calculations
    xlist,ylist = Ode(n,ksi,dksi,arrayspace)

    #  if using arrayspace cut-off the ylist at theta <= 0 to remove
    # unnecessary entries for Xi
    if arrayspace:
        finaly=[]
        finalx=[]
        for i in range(0,len(ylist)):
            finaly.append(ylist[i])
            finalx.append(xlist[i])
            if ylist[i] <=0:
                break
    # plot everything!
    ax1.plot(finalx,finaly,label="n = "+str(n) )
    finalrho = np.array(finaly)**n
    ax2.plot(finalx,finalrho,label="n = "+str(n))
ax1.legend()
ax2.legend()
plt.show()
