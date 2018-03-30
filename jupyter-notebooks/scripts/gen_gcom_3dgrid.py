"""
Name: gen_gcom_3dgrid.py
Description: Python script to generate 3d mesh from 2d
    + Read mesh dimensions from ProbSize.dat file
    + Build 3d levels using spline
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
import pytz
import datetime
from pytz import timezone
from datetime import datetime, timedelta

#from gccom_plot import *

#- Define source files
src_dir = './output'
fname = src_dir+'/gcom_grid2d.dat'
prob_size_file = src_dir+'/ProbSize_mont.dat'
out_fname = src_dir+'/gcom_grid3d.dat'

#- Figure out dimensions from ProbSize file
i = 0
f = open(prob_size_file, 'r') 
for line in f:  
    line = line.rstrip()
    if (i == 0):
        im1=int(line)
    elif (i==1):
        jm1=int(line)
    elif (i==2):
        km1=int(line)
    i = i + 1
f.close()

print("Mesh Dimensions: "+str(im1)+"x"+str(jm1)+"x"+str(km1))

#- Initialize
i=0
im = im1 * 2 + 1
jm0 = jm1 + 2
km0 = km1 + 2
imz = im1 + 2
jm = jm0 - 1
km = km0 - 1

#- Set up x,y,z arrays to be size:
#-     x(-1:im,0:jm0,0:km0)
#-     y(-1:im,0:jm0,0:km0)
#-     z(-1:im,0:jm0,0:km0)

x = np.empty((im1,jm1,km))
y = np.empty((im1,jm1,km))
z = np.empty((im1,jm1,km))
xi = np.empty((im1,jm1))
yi = np.empty((im1,jm1))
zi = np.empty((im1,jm1))
s = []
f = open(fname, 'r') 
for line in f:  
    line = line.rstrip()
    s.append(line)
f.close()

#- First get base 2D level 
theI = 0
for j in range(jm1):
    for i in range(im1):
        #- split on whitespace
        coords = s[theI].split()    #---  20   98  478242.7866000  3604379.1357500    -41.83000
        xi[i,j] = float(coords[2])
        yi[i,j] = float(coords[3])
        zi[i,j] = float(coords[4])
        theI = theI + 1

#- Now start building first/bottom k-level
for j in range(jm1):
    for i in range(im1):
        x[i,j,0] = xi[i,j]
        y[i,j,0] = yi[i,j]
        z[i,j,0] = zi[i,j]

#- Now building top k-level
for j in range(jm1):
    for i in range(im1):
        x[i,j,km-1] = x[i,j,0]
        y[i,j,km-1] = y[i,j,0]
        z[i,j,km-1] = -0.50

#- Now building up from k-level=2
for j in range(jm1):
    for i in range(im1):
        for k in range(1,km1):
            #delta = (k-1)/km1          #- Compute diff from k-1 to ceiling
            delta = (k)/km1              #- Compute diff from k-1 to ceiling
            xf = x[i,j,km1] - x[i,j,0]  #- Compute diff of x from bottom to top 
            yf = y[i,j,km1] - y[i,j,0]  #- Compute diff of y from bottom to top 
            zf = z[i,j,km1] - z[i,j,0]  #- Compute diff of z from bottom to top 
            x[i,j,k] = x[i,j,0] + xf*delta
            y[i,j,k] = y[i,j,0] + yf*delta
            z[i,j,k] = z[i,j,0] + zf*delta

#out_fname = src_dir+'/gcom_grid3d.dat'
out_fname = src_dir+'/Grid_mont.dat'
f = open(out_fname, 'w') 
for k in range(km):
    for j in range(jm1):
        for i in range(im1):
            f.write("%15.7f,%15.7f,%9.6f\n" % (x[i,j,k],y[i,j,k],z[i,j,k]))
f.close() 

#- Now fix ProbSize file, there is actually k+1 levels in output
outFile = 'ProbSize_mont.dat'
fid = open(src_dir+'/'+outFile,'w')
fid.write("%d\n%d\n%d\n1\n1\n1\n" % (im1,jm1,km1+1))
fid.close()

'''
#- Now format the this file for use in GCCOM
#- Note: prob can simply output this instead of above coords
xtmp = []
ytmp = []
ztmp = []
ai = []
aj = []
ak = []
ni = im1
nj = jm1
nk = km1
IM=ni
JM=nj
KM=nk 
f = open(out_fname, 'r') 
for line in f:  
    line = line.rstrip()
    coords = line.split()    #---  20   98  478242.7866000  3604379.1357500    -41.83000
    ai.append(coords[0])
    aj.append(coords[1])
    ak.append(coords[2])
    xtmp.append(float(coords[3]))
    ytmp.append(float(coords[4]))
    ztmp.append(float(coords[5]))
f.close()

#- set up figure
fig = plt.figure()
ax = fig.add_subplot(211)
plt.pcolormesh(x[:,:,1],y[:,:,1],z[:,:,1])

outFile = 'sdbay_gcominput.dat'
fid = open(src_dir+'/'+outFile,'w')

#- Iterate through each k-level and write the entire
#- Level to file ordered by i-j
for k in range(nk):
    ii = xtmp[k:(ni*nk)*nj:nk]
    jj = ytmp[k:(ni*nk)*nj:nk]
    kk = ztmp[k:(ni*nk)*nj:nk]
    for i in range(len(ii)):
        fid.write("%15.7f,%15.7f,%9.6f\n" % (ii[i],jj[i],kk[i]))
fid.close()
'''

'''

#- Now re-read and re-order
x = []
y = []
z = []
fid = open(src_dir+'/'+outFile,'r')
for line in fid:  
    line = line.rstrip()
    coords = line.split(',')    #---  20   98  478242.7866000  3604379.1357500    -41.83000
    x.append(float(coords[0]))
    y.append(float(coords[1]))
    z.append(float(coords[2]))
fid.close()

X = np.reshape(x,(ni,nj,nk))
Y = np.reshape(y,(ni,nj,nk))
Z = np.reshape(z,(ni,nj,nk))

xx = []
yy = []
zz = []
theI = 0
for k in range(nk):
#    print(k)
    for j in range(nj,0,-1):
        for i in range(ni,0,-1):
            xx.append(X[i-1,j-1,k])
            yy.append(Y[i-1,j-1,k])
            zz.append(Z[i-1,j-1,k])
            theI = theI + 1

outFile = 'Grid_sdbay.dat'
fid = open(src_dir+'/'+outFile,'w')
for i in range(len(xx)):
    fid.write("%15.7f,%15.7f,%9.6f\n" % (xx[i],yy[i],zz[i]))
fid.close()


#- Now fix ProbSize file, there is actually k-1 levels in output
outFile = 'ProbSize.dat'
fid = open(src_dir+'/'+outFile,'w')
fid.write("%d\n%d\n%d\n1\n1\n1\n" % (im1,jm1,km1-1))
fid.close()



ax = fig.add_subplot(212)
plt.plot(x,y,'.c')
plt.plot(x[0:ni+5],y[0:ni+5],'.r')
plt.plot(x[0],y[0],'og')

plt.show()
'''
