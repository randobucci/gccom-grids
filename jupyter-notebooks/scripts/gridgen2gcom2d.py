"""
Name: gridgen2gcom2d.py
Description: Python script to format gridgen output for use in GCCOM
    + Read mesh xy coordinates
    + Read bathy z coordinates
    + Output 2D mesh for input into GCCOM
* INPUT: ./output/output_bathy_CE, ./output/output_grid_CE
* OUTPUT: ./output/gcom_grid2d.dat, ./output/ProbSize_mont.dat
"""

import numpy as np
import urllib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors as c
from matplotlib import ticker
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import netCDF4
import os
import sys
import pytz
import datetime
from pytz import timezone
from datetime import datetime, timedelta


#- Define source files
src_dir = './output/'
fname_xy = src_dir+'/output_grid_CE'
fname_z = src_dir+'/output_bathy_CE'

#- Get xy-coords
i=0
x=[]
y=[]
f = open(fname_xy, 'r') 
for line in f:  
    line = line.rstrip()
    s = line
    if (i > 0):
        coords = s.split(' ')
        x.append(float(coords[0])) 
        y.append(float(coords[1])) 
    elif(i == 0):
        coords = s.split(' ')    #-- ## 50 x 50
        nx = int(coords[1])
        ny = int(coords[3])
    i = i + 1
f.close()

print ("\nNX = "+str(nx)+"; NY = "+str(ny))

#- Get z-coords
z=[]
f = open(fname_z, 'r') 
for line in f:  
    line = line.rstrip()
    z.append(float(line))
f.close()

#- Now loop through i-j indices and fill up x-y-z arrays
k = 0
xx = np.empty((nx,ny))
yy = np.empty((nx,ny))
zz = np.empty((nx,ny))
cutoff = -5.0
for j in range(ny):
    for i in range(nx):
        xx[i,j] = x[k]
        yy[i,j] = y[k]

        #- Fill up z-array with depth values, if
        #- less than cutoff value, fill with cutoff value
        if (z[k] < cutoff):
            zz[i,j] = z[k]
        else:
            zz[i,j] = cutoff
        k = k + 1

outFile = src_dir+'gcom_grid2d.dat'
f = open(outFile,'w')
for j in range(ny):
    for i in range(nx):
        #outData = [i+1,j+1,xx[i,j],yy[i,j],zz[i,j]]
        f.write("%4d %4d %15.7f %15.7f %15.7f\n" % (i+1,j+1,xx[i,j],yy[i,j],zz[i,j]))
f.close()
print('Created 2-D Grid File: '+outFile)


#- Now create a ProbSize File for input into GCCOM program
probSizeFile = src_dir+'ProbSize_mont.dat'
nk = 13
f = open(probSizeFile,'w')
f.write(str(nx)+"\n")
f.write(str(ny)+"\n")
f.write(str(nk)+"\n")
f.write("1\n1\n1\n")
f.close()

print('Created GCCOM ProbSize File: '+probSizeFile+'\n')
