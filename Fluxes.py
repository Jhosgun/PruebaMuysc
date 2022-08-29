#!/usr/bin/env python
# coding: utf-8

# In[28]:


# -*- encoding: latin1 -*-

"""-------------------------------------
MUOGRAPHY SIMULATION CODE,               

 _      _    ___  _ ____  ____    _      ____  ____  _____ _    
/ \__/|/ \ /\\  \/// ___\/   _\  / \__/|/  _ \/  _ \/  __// \      
| |\/||| | || \  / |    \|  /    | |\/||| / \|| | \||  \  | |   
| |  ||| \_/| / /  \___ ||  \_   | |  ||| \_/|| |_/||  /_ | |_/\
\_/  \|\____//_/   \____/\____/  \_/  \|\____/\____/\____\\____/
                                                                
 _____ _     _    ___  _ ____ 
/    // \   / \ /\\  \/// ___\
|  __\| |   | | || \  / |    \
| |   | |_/\| \_/| /  \ \___ |
\_/   \____/\____//__/\\\____/
                              



The Flux models  
  Jorge Jaimes, Jesus Peña  2021             
-------------------------------------

Updated: Frebrary 28, 2021
"""


"""
Full documentation at:
"""


import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from scipy import signal
from scipy.optimize import minimize
from scipy.interpolate import UnivariateSpline
import matplotlib as mpl

 
class Flux(object):
    
    def __init__(self,Elevation,matrixDatos,distances,h):
        
        
        self.cenith = matrixDatos[:,0]
        self.azimuth = matrixDatos[:,1]
        self.cenith_matrix = np.reshape(self.cenith,(50,50))
        self.z = np.abs(distances)
        self.h = h
        self.show()
        
    def int_flux(self,Lm, theta_in):
        """
        Aislada
        """

        # Minimum muon energy estimation

        E = np.linspace(1,1e4,100000) # GeV
        y = np.log10(E)
        l0 = 0.2549
        l1 = 0.0801
        l2 = 0.0368
        l3 = -0.0461
        l4 = 0.0154

        M = 10000
        Eu = 105.6 # muon mass [eV/c2]
        dE = (1e4-1)/M
        
        dEdp = -10**((l4*y**3.65) + (l3*y**3)  + (l2*y**2) + l1*y + l0)

        Op = -E*1e3/dEdp

        L = Lm*1e2 # Lenght to cm
        rho = 2.65 # Standard rock density g/cm3
        p = L*rho # Opacity g/cm2

        Eminp = E[np.argmin((Op - p)**2)] - Eu/1e3 

        # Muon flux model

        cenith = theta_in
        theta = cenith*np.pi/180.0 # cenith angle / redians
        E0 = np.linspace(1,1e4,10000) # Muon energy / GeV

        c = 1 # Speed of light
        m0 = 0.1056 # Muon mass in GeV/c^2
        p = np.sqrt((E0**2-m0**2*c**4)/c**2) # Bugaev model

        y = np.log10(p*np.cos(theta)) # Bugaev/Reyna model
        AB = 0.00253
        a0 = 0.2455
        a1 = 1.288
        a2 = -0.2555
        a3 = 0.0209
        
        #
        p = np.sqrt((E0)**2-((Eu/1e3)**2))
        h0 = 4900 + 750*p 
        Phi_Bugaev_Reyna = AB*(p**(-(a3*y**3 + a2*y**2 + a1*y + a0)))*(np.cos(theta))**3
        Phi_Bugaev_Reyna = np.exp(self.h/h0)*Phi_Bugaev_Reyna
        # Integrated flux estimation

        N = len(Phi_Bugaev_Reyna)
        Int_flux = 0
        Open_Sky = 0

        for i in range(N):
            Open_Sky = Open_Sky + Phi_Bugaev_Reyna[i]  # Open sky flux
            if E0[i] >= Eminp:
                Int_flux = Int_flux + Phi_Bugaev_Reyna[i] # Traversing flux

        Open_Sky = Open_Sky*dE
        Int_flux = Int_flux*dE 
        if math.isnan(Int_flux):
            print("Theta:",theta,"-",cenith)
        return  Int_flux*86400, Open_Sky*86400
    
    
    def Paso2(self,z,Elevation,cenith_matrix):
        Eu = 105.6 # muon mass [eV/c2]
        M = 10000
        dE = (1e4-1)/M
        N = 50
        Trav_Flux = np.zeros((N,N))
        Open_Sky_Flux = np.zeros((N,N))
        cenithal = np.zeros((N,N))

        for i in range(N):
            for j in range(N):
                cenithal[i,j] = (np.pi/2 - (Elevation + cenith_matrix[i,j]))*180/np.pi  # Cenith angle in grades
                if cenithal[i,j]>90:
                    cenithal[i,j]= 90

                A, B = self.int_flux(z[i,j], cenithal[i,j])


                Trav_Flux[i,j] = A
                Open_Sky_Flux[i,j] = B 
        
        self.cenithal = cenithal
        return Trav_Flux, Open_Sky_Flux, cenithal
    
    def show(self):
        
        azimuth = self.azimuth
        #cenithal = self.cenithal
        Trav_Flux, Open_Sky_Flux, cenithal = self.Paso2(self.z,np.radians(1.2217907712098335),self.cenith_matrix)
        FlujoIntegrado = Trav_Flux
        FlujoIntegrado[self.z==0] = 0

        fig = plt.figure(figsize=(20, 5))
        extent = (min(azimuth)*180/np.pi, max(azimuth)*180/np.pi, np.max(cenithal), np.min(cenithal))
        plt.imshow(FlujoIntegrado, interpolation='nearest', extent=extent, origin='upper',norm=mpl.colors.LogNorm(), cmap='jet')
        plt.xlabel("Azimuth [degree]", fontsize = 25)
        plt.ylabel("Zenith [degree]", fontsize = 25)


        # Color bar

        clb = plt.colorbar()
        clb.set_label('Integrated Flux [$cm^{-2}sr^{-1}day^{-1}$]', fontsize = 25)
        clb.ax.tick_params(labelsize = 20)

        labelsx = np.round(np.linspace(min(azimuth)*180/np.pi, max(azimuth)*180/np.pi, 11),0)
        labelsy = np.round(np.linspace(np.max(cenithal), np.min(cenithal),  11),0)

        plt.xticks(labelsx, fontsize = 15)
        plt.yticks(labelsy, fontsize = 15)

        plt.show()
        
        self.show2(FlujoIntegrado)
    
    
    def show2(self,FlujoIntegrado):
        
        fig = plt.figure(figsize=(20, 5))    
        plt.imshow(FlujoIntegrado, interpolation='nearest', origin='upper',norm=mpl.colors.LogNorm(), cmap='jet')

        # Color bar

        clb = plt.colorbar()
        clb.set_label('Integrated Flux [$cm^{-2}sr^{-1}day^{-1}$]', fontsize = 25)
        clb.ax.tick_params(labelsize = 20)
        print(FlujoIntegrado.max())
        #plt.xticks(labelsx, fontsize = 15)
        #plt.yticks(labelsy, fontsize = 15)

        plt.show()


