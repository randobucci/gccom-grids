#-----------------------------------------------------------------------------------
#- gccom_ploy.py
#- DESCRIPTION: A collection of functions to plot GCCOM model/grid data
#-----------------------------------------------------------------------------------
import numpy as np
import sys
import netCDF4
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors as c
from matplotlib import ticker
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d

def read_prob_size(prob_size_file):
    #- Figure out dimensions from ProbSize file
    i = 0
    f = open(prob_size_file, 'r')
    for line in f:
        line = line.rstrip()
        if (i == 0):
            IMax=int(line)
        elif (i==1):
            JMax=int(line)
        elif (i==2):
            KMax=int(line)
        i = i + 1
    f.close()
    return IMax,JMax,KMax

#- Read Grid from ascii file, return x,y,z vectors
def read_grid_ascii(grid_file,prob_size_file):
    #- Figure out dimensions from ProbSize file
    i = 0
    IMax,JMax,KMax = read_prob_size(prob_size_file)

    #- Read in grid data
    f = open(grid_file, 'r')
    x,y,z = [],[],[]
    for line in f:
        line = line.strip()
        cols = line.split(',')
        x.append(float(cols[0]))
        y.append(float(cols[1]))
        z.append(float(cols[2]))
    f.close()

    #- Reshape to 3-D array
    X = np.reshape(x,(IMax,JMax,KMax),order='F')
    Y = np.reshape(y,(IMax,JMax,KMax),order='F')
    Z = np.reshape(z,(IMax,JMax,KMax),order='F')
    return X,Y,Z,IMax,JMax,KMax
#-------------------------------------------------------------------


#---PLOT_BATHY------------------------------------------------------
#- INPUT: 3-D variables X,Y,Z
#-
#-------------------------------------------------------------------
def plot_bathy(X,Y,Z,level=None,fname=None):
    if(level == None):
        level = 0
    im_width = 5
    im_height = 5
    [nk,nj,ni] = X.shape        #- (10, 23, 98)
    zmin = np.min(Z)
    zmax = np.max(Z)

    #- Set up figure
    fig = plt.figure(figsize=(im_width, im_height),dpi=72)
    ax = fig.add_subplot(111)

    ##- Plot wireframe
    #ax.plot_wireframe(X[0,:,:], Y[0,:,:],Z[0,:,:] )
    ##- Plot surface
    print(zmin)
    print(zmax)
    heatmap = plt.pcolor(X[:,:,level], Y[:,:,level],Z[:,:,level],cmap=cm.jet, vmin=zmin,vmax=-3.0)
    #- colorbar legend
    cbar = plt.colorbar(heatmap)
    cbar.ax.set_ylabel('Depth (m)')

    ##- Label axes
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Monterey Curvilinear Mesh: Bathymetry K-Level '+str(level))
    if (fname == None):
        plt.show() 
    else:
        plt.savefig(fname,dpi=96)
        print("    Plotting: "+fname)
 

#---PLOT_3D------------------------------------------------------
#- INPUT: 3-D variables X,Y,Z
#-
#-------------------------------------------------------------------
def plot_3d(X,Y,Z,level=None):
    if(level == None):
        level = 0
    [ni,nj,nk] = X.shape        #- (10, 23, 98)
    zmin = np.min(Z[:,:,level])
    zmax = np.max(Z[:,:,level])

    #- Set up figure
    #fig = plt.figure()
    fig = plt.figure(figsize=plt.figaspect(0.5)*1.5)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X[:,:,level], Y[:,:,level],Z[:,:,level],cmap=cm.jet, \
                    vmin=zmin,vmax=zmax, rstride=1,cstride=1)
    #ax.plot_surface(X[:,:,nk-1], Y[:,:,nk-1],Z[:,:,nk-1],cmap=cm.jet, \
    #                vmin=zmin,vmax=zmax, rstride=1,cstride=1)
    #ax.scatter(X[1,1,1], Y[1,1,1],Z[1,1,1],c='y',s=100)
    ##- Label axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Curvilinear Mesh Level: K='+str(level))#+' and K='+str(nk) )
    #ax.set_title('Seamount Test Case: Bathymetry')

    plt.show()


src_dir = './output/'
prob_size_file = src_dir+'ProbSize_mont.dat'
grid_file = src_dir+'Grid_mont.dat'
X,Y,Z,imax,jmax,kmax = read_grid_ascii(grid_file,prob_size_file)
plot_3d(X,Y,Z)
