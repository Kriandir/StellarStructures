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
        self.rho = np.array
        self.phi = 0
        self.rhocenter = 0
        self.mass = 0

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
        self.rhocenter = (-1./3/self.phi*self.xlist[-1]**3)*rhoaverage

    # calculate mass for white dwarf
    def calcmass(self):
        self.mass = (-4./(np.sqrt(np.pi)))*(self.K/self.G)**(3./2)*self.phi
