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
import scipy.stats as ss
import Starobjects as so
from sklearn.metrics import cohen_kappa_score


#-----------------------------------------
# define starting values and initialize objects
nvalues = [0,1,(3./2),2,3,4]
initialtheta = 1.0
initialphi = 0.
dksi = 0.000001
end = 5
arraystep = 5
begin = dksi
Neutron = so.Stellar("neutron")
Dwarf = so.Stellar("dwarf")
# ----------------------------------------


# Function for plotting analytical solution
def anaplot(analist,color):
    ax1.plot(xlist,analist,linestyle="dashed",linewidth ="3",color = color,label="analytic n = "+str(n))
    ax3.plot(xlist,ylist-analist,linestyle="dashed",linewidth ="3",color = color,label="deviation for analytic and numeric n = "+str(n))
    # print "A sigma confidence level of: " + str(p_to_sigmas(1-ss.ttest_ind(ylist,analist)[1])) + "for n = " + str(n)

    # finding out how much the analytical solutions deviate from the numerical solutions
    ctr001 = 0
    ctr0001 = 0
    ctr00001 = 0
    ctr000001 = 0
    ctr0000001 = 0
    for i in range(len(ylist)):
        if abs(ylist[i] - analist[i])<= 0.001:
            ctr001 +=1
        if abs(ylist[i] -analist[i]) <= 0.0001:
            ctr0001 +=1
        if abs(ylist[i] - analist[i]) <= 0.00001:
            ctr00001 +=1
        if abs(ylist[i] - analist[i]) <= 0.000001:
            ctr000001 +=1
        if abs(ylist[i] - analist[i]) <= 0.0000001:
            ctr0000001 +=1

    print "------------------------------------"
    print "for n = " + str(n) + " and stepsize "+ str(dksi)
    print "------------------------------------"
    print " we have " +str(float(ctr001)/len(ylist)*100) + '%'+" of all the values of the analytical and numerical solution within 0.001"
    print " we have " +str(float(ctr0001)/len(ylist)*100) + '%'+" of all the values of the analytical and numerical solution within 0.0001"
    print " we have " +str(float(ctr00001)/len(ylist)*100) + '%'+" of all the values of the analytical and numerical solution within 0.00001"
    print " we have " +str(float(ctr000001)/len(ylist)*100) + '%'+" of all the values of the analytical and numerical solution within 0.000001"
    print " we have " +str(float(ctr0000001)/len(ylist)*100) + '%'+" of all the values of the analytical and numerical solution within 0.0000001"
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
fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=[20,20])
fig.suptitle("Plotted Graphs using a stepsize for "r'$\xi$'+" of " + str(dksi),fontsize = 30 )
fig.subplots_adjust(hspace=.5)
ax1.set_ylim([-0.1,1.1])
ttl = ax1.title
ttl.set_position([.5, 1.05])
ttl2 = ax2.title
ttl2.set_position([.5, 1.05])
ttl3 = ax2.title
ttl3.set_position([.5, 1.05])
ax1.set_title('Graph depicting '+ r'$\theta$'+ " dependence on " r'$\xi$',fontsize=20)
ax2.set_ylim([-0.1,1.1])
ax2.set_title('Graph depicting '+ r'$\frac{\rho}{\rho_c}$'+ " dependence on " r'$\xi$',fontsize=20)
ax3.set_title('Graph depicting the deviation of '+ r'$\theta$'+ "\'s dependence on " r'$\xi$',fontsize=20)
ax1.set_ylabel(r'$\theta$',fontsize=20,rotation=360)
ax2.set_ylabel(r'$\frac{\rho}{\rho_c}$',fontsize=25,rotation=360)
ax3.set_xlabel(r'$\xi$',fontsize=20)
ax3.set_ylabel(r'$\Delta\theta$',rotation=360,fontsize=20)
ax1.axhline(y=0,linestyle="dotted")
ax2.axhline(y=0,linestyle="dotted")
ax3.axhline(y=0,linestyle="dotted")
# loop through n values and integrate
for n in nvalues:
    solreal = []
    ksireal = []
    phireal = []

    # perform calculations
    xlist,ylist,phi = Ode(solreal,ksireal,phireal,n,dksi,begin,end,initialtheta,initialphi)
    finalrho = np.array(ylist)**n

    # if analytical or dwarf or neutronstar perform additional calculations.
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

    ax2.plot(xlist,finalrho,label="n = "+str(n))



# Perform calculations

Neutron.calcrhocenter()
print "---------------------------------------------------------------"
print "The density at the center of the neutronstar = " + str(Neutron.rhocenter)
Dwarf.calcmass()
print "The mass of the Dwarf star in solar mass = " + str(Dwarf.mass/(1.99*10**33))




ax1.legend(prop={'size':10})
ax2.legend(prop={'size':10})
ax3.legend(prop={'size':10})
plt.show()
