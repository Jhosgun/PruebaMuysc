import requests
import numpy as np
import os


# Path where topography data will be storage
path = 'C:/Users/Jorge/Desktop/MUYSC0.1/TopographyData'
os.chdir(path=path)

# Lat 1, Lat 2 , Longitude 1, Longitude 2
regionPoints = [4.500833, 4.466944, -75.404720, -76.370000, "CERRO MACHÍN"]


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



for i in downloadData(regionPoints):
    downloadLink = 'http://step.esa.int/auxdata/dem/SRTMGL1/'+i+'.SRTMGL1.hgt.zip'
    download_file(downloadLink)   
