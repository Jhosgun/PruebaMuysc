import Muysc
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from scipy import signal
import matplotlib as mpl
from srtm import Srtm1HeightMapCollection
srtm1_data = Srtm1HeightMapCollection()

#!export SRTM1_DIR=/home/jorge/Desktop/PruebaMuysc/TopographyData
 
path = '/home/jorge/Desktop/PruebaMuysc/TopographyData'

P1 = [4.492298, -75.381092]
RefPoint = [4.487717, -75.387880]
regionPoints = [4.466944, 4.500833,-75.370000,  -75.404720, "CERRO MACHÍN"]

cenit = [-10, 20,50]
azimut = [-10,22,50]

a = Muysc.Mute(regionPoints,80,srtm1_data,"jet")
a.downloadData(path)
a.pointView(P1,RefPoint)
a.plot_lines(cenit,azimut)         
a.show_distances()
a.section()
a.save_data("Prueba4")
a.IntegratedFlux()