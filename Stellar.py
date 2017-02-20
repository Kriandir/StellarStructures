##############################################
# Kriek van der Meulen
# Stellar Structures
# 2017
##############################################

from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

theta = 1
phi = 0
ksi = 0
dksi = 0.01
n = 1
result = []


def Ode(n,ksi,dksi,linspace):
    print n
    def Stellar(y,ksi,n,a):

        theta,phi = y
        dydt = phi/(ksi**2),-a*(ksi**2 * (np.power(theta,n)))
        return dydt
    y0 = [1.,0.]
    a=1

    def solfunc(ksi):
        sol = odeint(Stellar,y0,ksi,args=(n,a))
        return sol

    if linspace:
        ksi = np.arange(0.001,20,0.001)
        sol = odeint(Stellar,y0,ksi,args=(n,a))
        return ksi,sol[:,0]

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




def Euler(n,dksi):
    theta = 1
    phi = 0
    ksi = 0
    ylist = []
    xlist = []
    while theta > 0:
        ksi+=dksi
        theta += phi/(ksi**2) * dksi
        phi += -1*(ksi**2 * (theta**n))*dksi
        result.append([phi,theta,ksi])
        print ksi
        xlist.append(ksi)
        ylist.append(theta)
    return xlist,ylist


# print result
# Euler(5,0.1)
# plt.plot(xlist,ylist)
# plt.show()
fig,(ax1,ax2) = plt.subplots(2,1,figsize=[20,20])
ax1.set_xlabel(r'$\xi$',fontsize=18)
ax1.set_ylabel(r'$\theta$',fontsize=18,rotation=360)
ax2.set_xlabel(r'$\xi$',fontsize=18)
ax2.set_ylabel(r'$\frac{\rho}{\rho_c}$',fontsize=22,rotation=360)
#
# ax1.plot(Stars.instances[i].x,Stars.instances[i].my_model(Stars.instances[i].x,*Stars.instances[i].iniresult[j]),color= "red",linestyle="dotted",lw=3)
# ax2.plot(Stars.instances[i].x, (Stars.instances[i].y - Stars.instances[i].my_model(Stars.instances[i].x,*Stars.instances[i].iniresult[j]))/Stars.instances[i].y, 'o')

nvalues = [0,1,3./2,2,3,4]
for n in nvalues:
    ksi = [0.01]
    dksi = 0.01
    xlist,ylist = Ode(n,ksi,dksi,True)
    finaly=[]
    finalx=[]
    for i in range(0,len(ylist)):
        finaly.append(ylist[i])
        finalx.append(xlist[i])
        if ylist[i] <=0:
            break
    ax1.plot(finalx,finaly,label="n = "+str(n) )
    finalrho = np.array(finaly)**n
    ax2.plot(finalx,finalrho,label="n = "+str(n))
ax1.legend()
ax2.legend()
plt.show()
