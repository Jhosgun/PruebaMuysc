#!/usr/bin/env python
# coding: utf-8

# In[45]:



# -*- encoding: latin1 -*-

"""-------------------------------------
MUOGRAPHY SIMULATION CODE,               


  \  | |  |\ \  / __|  __|
 |\/ | |  | \  /\__ \ (   
_|  _|\__/   _| ____/\___|


The structre of MUTE    
  Jorge Jaimes, Jesus Peña  2021             
-------------------------------------

Updated: Frebrary 28, 2021
"""


"""
Full documentation at:
"""
#import import_ipynb

import os
import requests
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from scipy import signal
import matplotlib as mpl
import Fluxes

from srtm import Srtm1HeightMapCollection

srtm1_data = Srtm1HeightMapCollection()


# This class will make end-to-end simulation of the muography 
class Mute(object):

    
    # Input variables are expected to be in the following units:
    # Region        = [latitude1, latitude2, longitude1, longitude2, name]
    # srtm1_data    = Module 
    # Points        = Number of resolution
    
    def __init__(self,region,points,srtm1_data,cmap):
        """Initializing global variables 
        """
        
        self.lat1    = region[0] # First coordinates of latitude
        self.lat2    = region[1] # First coordinates of longitude
        self.lon1    = region[2] # Second coordinates of latitude
        self.lon2    = region[3] # Second coordinates of longitude
        self.name    = region[4] # Name of the study object
        
        
        self.srtm1_data = srtm1_data
        self.cmap = cmap
        self.points = points 
        self.regionPoints = region
        self.downloadData()
        self.elevation()

    def downloadData(self):
        # Path where topography data will be storage
        path = 'C:/Users/Jorge/Desktop/MUYSC0.1/TopographyData'
        os.chdir(path=path)


        def downloadData(regionPoints):
            """ This function download the data
            """
            if regionPoints[0] < 0:
                A = 'S'
            else:
                A = 'N'

            if regionPoints[1] < 0:
                B = 'S'
            else:
                B = 'N'
            
            if regionPoints[2] < 0:
                C = 'W'
            else:
                C = 'E'

            if regionPoints[3] < 0:
                D = 'W'
            else:
                D = 'E'

            # Calculate of decimals regions
            Lat1 = str(abs(int(regionPoints[0]))).zfill(2)
            Lat2 = str(abs(int(regionPoints[1]))).zfill(2)
            Log1 = str(abs(int(regionPoints[2]))).zfill(3)
            Log2 = str(abs(int(regionPoints[3]))).zfill(3)


            # Conditions of files
            if A+Lat1 == B+Lat2 and C+Log1 == D+Log2:
                topography_file = A + Lat1 + C + Log1 + '.SRTMGL1.hgt.zip'
            elif A+Lat1 == B+Lat2 and C+Log1 != D+Log2:
                topography_file1 = A + Lat1 + C + Log1 + '.SRTMGL1.hgt.zip'
                topography_file2 = A + Lat1 + D + Log2 + '.SRTMGL1.hgt.zip'
            elif A+Lat1 != B+Lat2 and C+Log1 == D+Log2:
                topography_file1 = A + Lat1 + C + Log1 + '.SRTMGL1.hgt.zip'
                topography_file2 = B + Lat2 + C + Log1 + '.SRTMGL1.hgt.zip'
            elif A+Lat1 != B+Lat2 and C+Log1 != D+Log2:
                topography_file1 = A + Lat1 + C + Log1 + '.SRTMGL1.hgt.zip'
                topography_file2 = B + Lat2 + D + Log2 + '.SRTMGL1.hgt.zip'
            
            # if there is differents coordinates
            if 'topography_file1' and 'topography_file2' in locals():
                return [topography_file1, topography_file2]
            else:
                return [topography_file]



        def download_file(downloadUrl):
            req = requests.get(downloadUrl)
            filename = req.url[downloadUrl.rfind('/')+1:]

            with open(filename, 'wb') as f:
                for chunk in req.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)



        for i in downloadData(self.regionPoints):
            downloadLink = 'http://step.esa.int/auxdata/dem/SRTMGL1/'+i+'.SRTMGL1.hgt.zip'
            download_file(downloadLink)   
        


    def measure(self, lat1, lon1, lat2, lon2):  
        """Converting GMS coordinates to local meters.
        """
        
        R = 6378.137                                                      # Radius of earth in KM
        dLat = lat2 * np.pi / 180 - lat1 * np.pi / 180                    # Distance of the Latitude coordinates
        dLon = lon2 * np.pi / 180 - lon1 * np.pi / 180                    # Distance of the Longitude coordinates 
        a = np.sin(dLat/2) * np.sin(dLat/2) + np.cos(lat1 * np.pi / 180)  * np.cos(lat2 * np.pi / 180) * np.sin(dLon/2) * np.sin(dLon/2)
        c = 2 * math.atan2(np.sqrt(a), np.sqrt(1-a))
        d = R * c                                                         # Total distance
        return d * 1000 # meters

    
    def elevation(self):
        """Calculating the elevation of the geological structure.
        """
        
        op=0
        # List of latitude and longitude
        lisLat = np.linspace(self.lat1, self.lat2, num=self.points)
        lisLong = np.linspace(self.lon1, self.lon2, num=self.points)
        
        # List of latitude and longitude in meteres
        meterspointLat = self.measure(self.lat1, self.lon1, self.lat2, self.lon1)
        meterspointLong = self.measure(self.lat1, self.lon1, self.lat1, self.lon2)
        
        if self.lon1 and self.lon2 <0:
            meterspointLong = -meterspointLong
        # List from the initial point ot     
        lisLatM = np.linspace(0,meterspointLat, num=self.points) 
        lisLongM = np.linspace(0,meterspointLong, num=self.points) 

        if op==0:
            self.X, self.Y = np.meshgrid(lisLong, lisLat) # Creation grid in point with meters
        if op==1:
            self.X, self.Y = np.meshgrid(lisLongM, lisLatM)

        
        # Calculated the elevation
        lisZ=np.zeros([self.points,self.points])
        for i in range(self.points):
            for j in range(self.points):
                lisZ[i][j]=srtm1_data.get_altitude(latitude=lisLat[i], longitude=lisLong[j])
        self.lisZ = lisZ
        

        
    def plot_structure(self): 
        """Plot the geological structure
        """
        
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(20,17)) # Geological structure
        surf = ax.plot_surface(self.X, self.Y, self.lisZ, edgecolors='grey', alpha=0.5,  cmap=self.cmap,
                            linewidth=0, antialiased=False)
        
        
        # Labels                                            
        ax.set_xlabel("Longitud", fontsize=22)
        ax.set_ylabel("Latitud", fontsize=22)
        ax.set_zlabel("Elevación", fontsize=22)
        ax.set_title(self.name, fontsize=25)
        plt.show()
        
        
    def pointView(self, ObsPoint, RefPoint):
        """ Set the observation point and reference point
        """
        
        # Observation point
        self.obsPX = ObsPoint[1]
        self.obsPY = ObsPoint[0]
        self.obsPZ = srtm1_data.get_altitude(latitude=self.obsPY, longitude=self.obsPX)
        
        # Reference point
        if RefPoint=="max":
            self.maxpoint = np.where(self.lisZ==self.lisZ.max())
            self.RefPX = self.X[self.maxpoint[0],self.maxpoint[1]]
            self.RefPY = self.Y[self.maxpoint[0],self.maxpoint[1]]
            self.RefPZ = self.lisZ.max()
        else:
            self.RefPX=RefPoint[1]
            self.RefPY=RefPoint[0]
            self.RefPZ=srtm1_data.get_altitude(latitude=self.RefPY, longitude=self.RefPX)
        
        
        ca = self.measure(self.obsPY,self.obsPX,self.RefPY,self.RefPX)
        co = self.RefPZ
        self.angle = np.arctan(co/ca)
   
        print("El angulo: ",self.angle)
       
        
        
        
    def plot_lines(self, cenit, azimut): 
        """ Plot the geological structure with the moun's lines
        """
        
        cenitF = cenit[0]
        cenitS = cenit[1]
        cenitP = cenit[2]
        
        azimutF = azimut[0]
        azimutS = azimut[1]
        azimutP = azimut[2]
        
        # Equation creation
        t = 4
        self.PjecX = self.obsPX+(self.RefPX-self.obsPX)*t
        self.PjecY = self.obsPY+(self.RefPY-self.obsPY)*t
        self.PjecZ = self.obsPZ+(self.RefPZ-self.obsPZ)*t
        
        
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(20,17)) 
        surf = ax.plot_surface(self.X, self.Y, self.lisZ, edgecolors='grey', alpha=0.5,  cmap=self.cmap,
                            linewidth=0, antialiased=False, label='Cerro machín')
        
        # Plot points
        ax.scatter(self.RefPX,self.RefPY,self.RefPZ, c="r")          # Reference point
        ax.scatter(self.obsPX, self.obsPY, self.obsPZ, c="r")        # Telescope point
        ax.plot([self.obsPX,self.PjecX],[self.obsPY,self.PjecY],     # Line projection                                                  
                [self.obsPZ,self.PjecZ], color="red")                      
        
        
        # Plot lines of muons

        cenit = np.linspace(cenitS,cenitF,cenitP)
        azimut = np.linspace(azimutF,azimutS,azimutP)
        self.azimut = azimut
        self.cenit = cenit # corregir variable / creada para el funcionamiento de h flujo
        self.distances = np.zeros((cenitP,azimutP))
        t = np.linspace(0,1,1000)
        for k,i in enumerate(azimut):
            for m,j in enumerate(cenit):
                
                self.pointfX = self.obsPX + math.cos(i * np.pi / 180) * (self.PjecX - self.obsPX) - math.sin(i * np.pi / 180) * (self.PjecY - self.obsPY)
                self.pointfY = self.obsPY + math.sin(i * np.pi / 180) * (self.PjecX - self.obsPX) + math.cos(i * np.pi / 180) * (self.PjecY - self.obsPY)
                self.pointfZ = self.PjecZ + math.sin(j * np.pi / 180) * self.PjecZ 
                #=======================================================
                self.equationX = self.obsPX+(self.pointfX-self.obsPX)*t
                self.equationY = self.obsPY+(self.pointfY-self.obsPY)*t
                self.equationZ = self.obsPZ+(self.pointfZ-self.obsPZ)*t
                
                # Call function calculate distance and add to matrix distances 
                self.distances[m][k]= self.calculate_distance(self.equationX,self.equationY,self.equationZ)
                ax.plot(self.equationX,self.equationY,self.equationZ, color="black") # Line equation
             
        # Labels                                            
        ax.set_xlabel("Longitud", fontsize=22)
        ax.set_ylabel("Latitud", fontsize=22)
        ax.set_zlabel("Elevación", fontsize=22)
        ax.set_title("TOPOGRAFÍA "+self.name, fontsize=25)
        plt.show()
        
        
        
                
        Lmcenit = []
        Lmazimut = []
        for i in self.cenit:
            for j in self.azimut:
               Lmcenit.append(i)
               Lmazimut.append(j)

        matrixDatos = np.empty((self.distances.shape[0]*self.distances.shape[1],3))        
        matrixDatos[:,0] = np.radians(Lmcenit)
        matrixDatos[:,1] = np.radians(Lmazimut)
        matrixDatos[:,2] = self.distances.flatten()
        
        self.matrixDatos = matrixDatos
        
        return self.distances,cenit,azimut,[self.obsPX,self.obsPY,self.obsPZ],[self.RefPX,self.RefPY,self.RefPZ]
    
    
    
    def show_Convolution(self):

        E0 = np.linspace(0.01, int(1E4), 10000000)
        self.Densidad = self.distances * 2.65 * 100
        
        Model = "ReynaB"
  
        a = np.zeros(shape=(self.distances.shape[0],self.distances.shape[1]))
        a[:,:] = 5.237486
        cenitTrans = 90 - (self.cenit + self.angle)
        cenitTrans[cenitTrans > 90] = 90
        
        
        matrixCenit = np.zeros(shape=(len(cenitTrans),len(self.azimut)))
        cenitmat = np.zeros(shape=(len(cenitTrans),len(self.azimut)))
        for i in range(self.distances.shape[1]):
            matrixCenit[:,i] = cenitTrans
            cenitmat[:,i] = self.cenit
        

        
        
        
        Flujo = ModelsFlux.Fluxs(E0,Model,self.obsPZ,matrixCenit,self.distances)
        A = Flujo.IntegratedFlux()
    

        
        fig = plt.figure(figsize=(20,17))
        ax= fig.add_subplot(311)
        
        # Labels
        ax.set_xlabel('  'r'$\alpha$[Acimutal]', fontsize=22)
        ax.set_ylabel('  'r'$\phi$[Cenital]', fontsize=22)
        ax.set_title("DISTANCIA RECORRIDA", fontsize=25)
        
        plt.imshow(A, interpolation='nearest', origin='upper',norm=mpl.colors.LogNorm())
        # Color bar
        clb = plt.colorbar()
        clb.set_label('d [km]',fontsize = 25)
        clb.ax.tick_params(labelsize = 20)
        
        
        # Center location
        plt.axvline(x=0, color='k', lw=1, linestyle='--')
        plt.axhline(y=0, color='k', lw=1, linestyle='--')
        
        # Increse tick labels
        ax.tick_params(axis='both', which='major',labelsize=15)
        
        # Save the image
        # plt.savefig()
        plt.show()
        
        
        return cenitmat, matrixCenit
        
        

    
    
    def calculate_distance(self,equationX,equationY,equationZ):
        """ This function calculate the distance for each point
        """
        d=0
        indexlist = []
        
        self.projection=np.zeros(len(equationX))
        for i in range(len(equationX)):
            self.projection[i] =srtm1_data.get_altitude(latitude=equationY[i], longitude=equationX[i])
    
        coor = equationZ<self.projection
        
        for j,i in enumerate(coor):
            if i==True:
                if j==0:
                    if coor[j+1]==True:
                        indexlist.append(j)
                else: 
                    if coor[j+1]==False and coor[j-1]==True:
                        indexlist.append(j)
                    elif coor[j-1]==False and coor[j+1]==False:
                        dtemp = self.measure(equationY[j-1],equationX[j-1],equationY[j+1],equationX[j+1])
                        d+=dtemp/3
                    elif coor[j-1]==False:
                        indexlist.append(j)
                        
        x,y = equationX[indexlist],equationY[indexlist]
        listp = list(zip(indexlist[0::2], indexlist[1::2]))
        
        for i,j in enumerate(listp):
            dtemp = self.measure(equationY[listp[i][0]],equationX[listp[i][0]],equationY[listp[i][1]],equationX[listp[i][1]])
            d+=dtemp
        return d
    
    
    
    def point_show(self):
        """ This function is on progress to show multiple points
        """
        
        #a = self.lisZ==self.obsPZ
        zmin = self.obsPZ-10
        zmax = self.obsPZ+10
        indexX = []
        indexY = []
        

        # List of latitude and longitude
        lisLat = np.linspace(self.lat1, self.lat2, num=self.points)
        lisLong = np.linspace(self.lon1, self.lon2, num=self.points)
        
        # List of latitude and longitude in meteres
        meterspointLat = self.measure(self.lat1, self.lon1, self.lat2, self.lon1)
        meterspointLong = self.measure(self.lat1, self.lon1, self.lat1, self.lon2)
        
        
        self.X, self.Y = np.meshgrid(lisLong, lisLat) # Creation grid in point with meters
        
        # Calculated the elevation
        lisZ=np.zeros([self.points,self.points])
        for i in range(self.points):
            for j in range(self.points):
                lisZ[i][j]=srtm1_data.get_altitude(latitude=lisLat[i], longitude=lisLong[j])
                if self.lisZ[i][j]>zmin and self.lisZ[i][j]<zmax:
                    indexX.append(lisLat[i])
                    indexY.append(lisLong[j])
        return indexX,indexY
        
        
    def save_data(self,name):
        """ This function save .dat data 
        """
        
        np.savetxt(name+'.dat', self.matrixDatos)
        
        
        # Add first 4 lines
        
        text_list = ["# Op " + str(self.obsPX) + " " + str(self.obsPY) + " " + str(self.obsPZ) + "\n", 
                     "# Pp " + str(self.RefPX) + " " + str(self.RefPY) + " " + str(self.RefPZ) + "\n",
                     "# A "  + str(self.lat1)  + " " + str(self.lat2)  + " " + str(self.lon1) + " " + str(self.lon2) + "\n",
                     "# cenith[rad] azimuth[rad] distance[km]\n"]
      

        
        
        
        
        my_file = open(name+".dat")
        string_list = my_file.readlines()
        my_file.close()

        newstring_list = text_list + string_list

        my_file = open(name+".dat", "w")
        new_file_contents = "".join(newstring_list)


        my_file.write(new_file_contents)
        my_file.close()

        
        
    def section(self):
        """ This function plot the section of the geolical structure traversed
        """
        
        fig, ax = plt.subplots(figsize=(20,17))
        ax.fill_between(self.equationX, self.projection, y2=0, alpha=0.5, color="navy")
        ax.plot(self.equationX,self.equationZ, color="red")
        
        # Labels
        #ax.set_ylim([2200,3000])
        ax.set_xlabel("Latitud", fontsize=25)
        ax.set_ylabel("Elevación", fontsize=25)
        ax.set_title('CERRO MACHÍN SECCIÓN', fontsize=25)                                         
        plt.show()


        
        
    def show_distances(self):
        """ This function plot the distance traversed in geological structured
        """
        fig = plt.figure(figsize=(20,17))
        extent = (min(self.azimut), max(self.azimut), min(self.cenit),max(self.cenit))
        plt.imshow(self.distances, interpolation='nearest', extent=extent, origin='upper', cmap=self.cmap)
        plt.xlabel("Azimuth [degree]", fontsize = 25)
        plt.ylabel("Zenith [degree]", fontsize = 25)
        plt.title("DISTANCIA RECORRIDA "+self.name, fontsize = 25)

        # Color bar
        clb = plt.colorbar()
        clb.set_label('d [m]', fontsize = 25)
        clb.ax.tick_params(labelsize = 20)

        labelsx = np.round(np.linspace(min(self.azimut), max(self.azimut), 11),0)
        labelsy = np.round(np.linspace(min(self.cenit), max(self.cenit), 11),0)

        plt.xticks(labelsx, fontsize = 15)
        plt.yticks(labelsy, fontsize = 15)

        plt.show()
        
     
   

    def IntegratedFlux(self):
        t = Fluxes.Flux(self.angle,self.matrixDatos,self.distances,self.obsPZ)
        



