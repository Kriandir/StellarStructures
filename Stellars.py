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
from scipy import stats

# Make class structure for performing calculations for the white dwarf and neutron star

class Stellar:
    def __init__(self,name):
        self.c = 2.99792458*10**10
        self.G = 6.67259*10**(-8)
        self.hbar = 1.05457266*10**(-27)
        self.mp = 1.6726231*10**(-24)
        self.K = 0
        self.xlist = np.array
        self.ylist = np.array
        self.r = np.array
        self.rmax = 0
        self.mass = 0
        self.rho = np.array
        self.phi = 0
        self.rhocenter = 0

        if name == "dwarf":
            self.Dwarf()
        if name == "neutron":
            self.Neutron()

    # appenddate to our objects
    def appenddata(self,xlist,ylist,phi):
        self.xlist = xlist
        self.ylist = ylist
        self.phi = phi

    # define usefull values for Dwarfstar
    def Dwarf(self):
        self.K = (3**(1./3)*np.pi**(2./3)*self.hbar*self.c)/(2**(4./3)*4*self.mp**(4./3))
        self.n = 3

    # define usefull values for Neutronstar
    def Neutron(self):
        self.n = 1
        self.mass = 1.4*1.99*10**33
        self.rmax = 1*10**6

    # calculate rhocenter for neutronstar
    def calcrhocenter(self):
        rhoaverage = self.mass/((4./3)*np.pi*self.rmax**3)
        self.rhocenter = (-1./3/self.phi*xlist[-1]**3)*rhoaverage

    # calculate mass for white dwarf
    def calcmass(self):
        print self.phi
        self.mass = (-4./(np.sqrt(np.pi)))*(self.K/self.G)**(3./2)*self.phi




#-----------------------------------------
# define starting values
initialtheta = 1.0
initialphi = 0.
dksi = 0.001
end = 5
arraystep = 5
begin = dksi
Neutron = Stellar("neutron")
Dwarf = Stellar("dwarf")
# ----------------------------------------



# Function for plotting analytical solution
def anaplot(analist,color):
    ax1.plot(xlist,analist,linestyle="dashed",linewidth ="3",color = color,label="analytic n = "+str(n))
    anafinalrho = np.array(analist)**n
    ax2.plot(xlist,anafinalrho,linestyle="dashed",linewidth ="3",color = color,label="analytic n = "+str(n))

    ax3.plot(xlist,ylist-analist,linestyle="dashed",linewidth ="3",color = color,label="numeric n - analytic n for n = "+str(n))
    ax4.plot(xlist,finalrho-anafinalrho,linestyle="dashed",linewidth ="3",color = color,label="numeric n - analytic n for n = "+str(n))
    print stats.ttest_ind(ylist,analist)

# Function for ode integration
def Ode(solreal,ksireal,phireal,n,dksi,begin,end,initialtheta,initialphi):

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
    phireal = np.append(phireal,sol[:,1])


    # Check whether we found a solution and return

    for i in range(0,len(solreal)):
        if solreal[i] <=0 or isnan(solreal[i]):
            return ksireal[:i],solreal[:i],phireal[i]

    # call Ode on the stack with new initial values and return when solution is found
    begin = end
    end = end + arraystep
    return Ode(solreal,ksireal,phireal,n,dksi,begin,end,initialtheta,initialphi)


# define figure


fig,(ax1,ax2) = plt.subplots(2,1,figsize=[20,20])
ax1.set_ylim([-0.1,1.1])
ax2.set_ylim([-0.1,1.1])
ax1.set_xlabel(r'$\xi$',fontsize=18)
ax1.set_ylabel(r'$\theta$',fontsize=18,rotation=360)
ax2.set_xlabel(r'$\xi$',fontsize=18)
ax2.set_ylabel(r'$\frac{\rho}{\rho_c}$',fontsize=22,rotation=360)

fig2,(ax3,ax4) = plt.subplots(2,1,figsize=[20,20])
ax3.set_ylim([-0.000006,0.000006])
ax4.set_ylim([-0.000006,0.000006])
ax3.set_xlabel(r'$\xi$',fontsize=18)
ax3.set_ylabel(r'$difference$ $in$ $\theta$',fontsize=18)
ax4.set_xlabel(r'$\xi$',fontsize=18)
ax4.set_ylabel(r'$difference$ $in$ $\frac{\rho}{\rho_c}$',fontsize=22)

# define nvalues to loop over and check whether using arrayspace
nvalues = [1,3]
arrayspace = True

for n in nvalues:
    solreal = []
    ksireal = []
    phireal = []

    # perform calculations
    xlist,ylist,phi = Ode(solreal,ksireal,phireal,n,dksi,begin,end,initialtheta,initialphi)
    finalrho = np.array(ylist)**n


    if n == 0:
        analist = (1-(1./6)*xlist**2)
        color = "red"
        anaplot(analist,color)

    if n == 1:
        analist = np.sin(xlist)/xlist
        color = "blue"
        anaplot(analist,color)
        Neutron.appenddata(xlist,ylist,phi)

    if n == 3:
        Dwarf.appenddata(xlist,ylist,phi)

    # plot everything!
    ax1.plot(xlist,ylist,label="n = "+str(n) )
    ax1.axhline(y=0,linestyle="dotted")
    ax2.plot(xlist,finalrho,label="n = "+str(n))
    ax2.axhline(y=0,linestyle="dotted")


# Perform calculations

Neutron.calcrhocenter()
print Neutron.rhocenter
Dwarf.calcmass()
print Dwarf.mass/(1.99*10**33)




ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()
plt.show()
