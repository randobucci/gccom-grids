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

from gccom_utils import *


#---PLOT_VECTOR_POINTS------------------------------------------------------
#- INPUT: 2 vector variables x,y
#- OUTPUT: Simple points of bottom most level in 2D
#-------------------------------------------------------------------
def plot_vector_points(x,y,IM,JM,fname=None):
    im_width,im_height = 5.0,5.0
    [nx,ny] = len(x),len(y)
    fig = plt.figure(figsize=(im_width, im_height),dpi=72)
    fig.set_size_inches(im_width, im_width)
    ax = fig.add_subplot(111)
    #- Plot all the points
    plt.plot(x,y,'.k')
    #- Plot the first i=0:ni how they are ordered
    plt.plot(x[0:IM],y[0:IM],'b')
    plt.plot(x[0:IM],y[0:IM],'.b')
    #- Plot (0,0) point
    plt.plot(x[0],y[0],'og',markersize=15)
    #- Plot (ni,0) point
    plt.plot(x[IM-1],y[IM-1],'or',markersize=15)

    ##- Label axes
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('San Diego Bay Mesh 2D Points')
    if (fname == None):
        plt.show() 
    else:
        plt.savefig(fname,dpi=128)
        print("    Plotting: "+fname)

#---PLOT_INDICES------------------------------------------------------
#- INPUT: 2 mesh 3-D variables x,y
#-
#-------------------------------------------------------------------
def plot_indices(X,Y,fname=None):
    im_width,im_height = 10,10
    [nx,ny] = len(X),len(Y)
    fig = plt.figure(figsize=(im_width, im_height),dpi=72)
    #- First plot real world coordinates
    ax = fig.add_subplot(211)
    #- First plot all the points
    plt.plot(X[:,:,0],Y[:,:,0],'.c')
    #- Next plot i=0 points
    plt.plot(X[:,0,0],Y[:,0,0],'.r')
    #- Plot (0,0) point
    plt.plot(X[0,0,0],Y[0,0,0],'og')
    ##- Label axes
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('San Diego Bay Mesh 2D Points')

    #- Now plot computational space
    ax = fig.add_subplot(212)
    #- Plot origin
    xcolor = 'blue'
    ycolor = 'yellow'
    start_color = 'lime'
    end_color = 'red'

    xmin = 0-5
    xmax = nx+5
    ymin = 0-5
    ymax = ny+5
    compx = np.linspace(0,nx-1,nx)
    compy = np.linspace(0,ny-1,ny)
    yv,xv = np.meshgrid(compy,compx)
    zv = np.ones((nx,ny),order='F')
    plt.axis([xmin,xmax,ymin,ymax])
    plt.fill([xmin,xmin,xmax,xmax],[ymin,ymax,ymax,ymin],'gray')
    plt.pcolormesh(xv, yv, zv, facecolor='none', edgecolors='k', 
               alpha=0.8, cmap='gray',zorder=2)
    plt.plot(xv[0:nx,0],yv[0:nx,0],c=xcolor)
    plt.plot(xv[0:nx,0],yv[0:nx,0],'.',c=xcolor,markersize=12)
    plt.plot(xv[0,0:ny],yv[0,0:ny],c=ycolor)
    plt.plot(xv[0,0:ny],yv[0,0:ny],'.',c=ycolor,markersize=12)
    plt.plot(xv[0,0],yv[0,0],'o',c=start_color,markersize=14)
    plt.plot(xv[0,ny-1],yv[0,ny-1],'s',c=end_color,markersize=14)
    plt.plot(xv[nx-1,0],yv[nx-1,0],'s',c=end_color,markersize=14)
    plt.xlabel('I')
    plt.ylabel('J')


    if (fname == None):
        plt.show()
    else:
        plt.savefig(fname,dpi=96)
        print("    Plotting: "+fname)



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
    heatmap = plt.pcolor(X[:,:,level], Y[:,:,level],Z[:,:,level],cmap=cm.jet, vmin=zmin,vmax=zmax)
    #- colorbar legend
    cbar = plt.colorbar(heatmap)
    cbar.ax.set_ylabel('Depth (m)')

    ##- Label axes
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('La Jolla Curvilinear Mesh: Bathymetry K-Level '+str(level))
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
    ax.plot_surface(X[:,:,nk-1], Y[:,:,nk-1],Z[:,:,nk-1],cmap=cm.jet, \
                    vmin=zmin,vmax=zmax, rstride=1,cstride=1)
    #ax.scatter(X[1,1,1], Y[1,1,1],Z[1,1,1],c='y',s=100)
    ##- Label axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('San Diego Bay Curvilinear Mesh: Bathymetry')
    #ax.set_title('Seamount Test Case: Bathymetry')

    plt.show()



