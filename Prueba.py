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


ObsPoint = [37.744571, 14.997342]
RefPoint = [37.747857, 14.999084]


regionPoints = [37.738283, 37.781430,14.964277, 15.026140, "MONTE ETNA"]

cenit = [-5, 20,50]
azimut = [-25,15,50]

a = Muysc.Mute(regionPoints,80,srtm1_data,"jet")
a.pointView(ObsPoint,RefPoint)
a.plot_lines(cenit,azimut)         
a.show_distances()
a.section()
a.save_data("Prueba4")
a.IntegratedFlux()