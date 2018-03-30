#-------------------------------------------------------
#- File: create_grid.py
#- Description: Create curvilinear grid using PYGRIDGEN
#-    + Open file w/ geographic coords
#-        + Project geo coords to utm
#-        + Create utm file
#-    + Open vertex file and find corresponding lat/lon 
#-    + Create grid file with 4 vertices labeled 1
#- Author: Randy Bucciarelli
#- Date: June 15, 2017
#-------------------------------------------------------
import warnings
warnings.simplefilter('ignore')
import numpy as np
from matplotlib import pyplot as plt
import sys
import seaborn
clear_bkgd = {'axes.facecolor':'none', 'figure.facecolor':'none'}
seaborn.set(style='ticks', context='notebook', rc=clear_bkgd)
import pygridgen
from pyproj import Proj
from pyproj import Geod

#- Define function to plot results
def plot_grid(grid, ax):
    ax.plot(grid.x.flatten(), grid.y.flatten(), 'k.', label='Grid nodes', zorder=5)
    ax.set_aspect('equal')
    #ax.set_xlim([0, 4])
    #ax.set_ylim([0, 4])
    ax.plot(grid.xbry, grid.ybry, '-', color='0.5', zorder=0)
    pos = np.nonzero(grid.beta == 1)
    neg = np.nonzero(grid.beta == -1)
    #ax.plot(x[pos], y[pos], 'go', label='Positive', zorder=2, alpha=0.5)
    #ax.plot(x[neg], y[neg], 'rs', label='Negative', zorder=2, alpha=0.5)
    ax.legend(numpoints=1)

def plot_mesh(xi,yi,ax):
    #- Set up plot
    xcolor='blue'
    ycolor='yellow'
    start_color='lime'
    end_color='red'
    ni,nj = xi.shape
    zi = np.ones((ni,nj),order='F')
    ax.pcolormesh(xi, yi, zi, facecolor='none', edgecolors='k',
               alpha=0.2, cmap='gray',zorder=2)
    ax.plot(xi[0,0],yi[0,0],'o',c=start_color,markersize=14) #- Green Origin (0,0)
    ax.plot(xi[0,0:nj],yi[0,0:nj],c=ycolor)
    ax.plot(xi[0,nj-1],yi[0,nj-1],'s',c=end_color,markersize=14)  #- Y-axis end
    ax.plot(xi[0:ni,0],yi[0:ni,0],c=xcolor)
    ax.plot(xi[ni-1,0],yi[ni-1,0],'s',c=end_color,markersize=14)  # X-axis end


#- Set file for simple polygon in geographic coords
#- Define source files
input_args = sys.argv
if (len(sys.argv) < 2):
    src_dir = './input/'
    geo_file = src_dir+'mont-outline_geo.xy'
else:
    geo_file = input_args[1]


#- Initialize
lons,lats,verts=[],[],[]
x,y=[],[]
vert_lons,vert_lats,vert_ids = [],[],[]  #- vertices array

#- Read in vertex positions in mont-vertices.txt; should be 'lon lat'
vert_file = './input/mont-vertices.txt'
file = open(vert_file)
for line in file:  
    line = line.rstrip()
    coords = line.split()
    vert_lons.append(float(coords[0]))
    vert_lats.append(float(coords[1]))
    vert_ids.append(float(coords[2]))

#- Open kml file and read into memory
file = open(geo_file)
for line in file:  
    line = line.rstrip()
    coords = line.split()
    lon = float(coords[0])
    lat = float(coords[1])
    lons.append(lon)
    lats.append(lat)
    verts.append(0.0)
file.close()

#- Project data from geo to utm
myProj = Proj("+proj=utm +zone=10S, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
geod = Geod(ellps='WGS84') 
EarthRadius = 6371000 #- meters
x, y = myProj(lons, lats)

#- Loop through vertices and find closest lat/lon index
for v in range(4):     #- only 4 vertices
    vert_lon = vert_lons[v]
    vert_lat = vert_lats[v]
    distances = []
    for i in range(len(lats)):
        #angle1,angle2,distance = geod.inv(long1, lat1, long2, lat2)
        angle1,angle2,dist = geod.inv(vert_lon, vert_lat, lons[i], lats[i])
        distances.append(dist)
    min_index = distances.index(min(distances))
    verts[min_index] = 1
    
#- Write utm w/ correct vertex out to file
utm_file = geo_file.replace('geo','utm')
thefile = open(utm_file, 'w')
for i in range(len(x)):
    thefile.write("%13.4f\t%13.4f\t%d\n" % (x[i],y[i],verts[i]))
thefile.close()

#- Now create grid using pygridgen
x = np.array(x)
y = np.array(y)
beta = np.array(verts)

#- Call gridgen using pygridgen. Usage:
#- pygridgen.Gridgen(xbry, ybry, beta, shape, focus=None)
ni,nj = 50,25

#- Example 1- 10x50
grid = pygridgen.Gridgen(x, y, beta, shape=(nj, ni))


#- Need to reorder grid indices to match the ordering of digitized boundary
xi = grid.x.T    #- Take transpose of grid x coordinates
yi = grid.y.T    #- Take transpose of grid y coordinates
xi = np.fliplr(xi)#.flatten()
yi = np.fliplr(yi)#.flatten()
xi = np.flipud(xi)#.flatten()
yi = np.flipud(yi)#.flatten()

#- Do same for grid center points
xi_rho = grid.x_rho.T    #- Take transpose of grid x coordinates
yi_rho = grid.y_rho.T    #- Take transpose of grid y coordinates
xi_rho = np.fliplr(xi_rho)#.flatten()
yi_rho = np.fliplr(yi_rho)#.flatten()
xi_rho = np.flipud(xi_rho)#.flatten()
yi_rho = np.flipud(yi_rho)#.flatten()


#- Set up plot
xcolor='blue'
ycolor='yellow'
start_color='lime'
end_color='red'

fig, axes = plt.subplots(nrows=1,ncols=2,figsize=(7, 7))
ax = axes.flatten()[0]
plot_mesh(xi,yi, ax)
nx,ny = ni,nj
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.title.set_text('Curvilinear Grid: '+str(nx)+' x '+str(ny))

#- Now plot computational space
ax = axes.flatten()[1]
nx,ny = ni,nj
xmin = 0-2
xmax = nx+2
ymin = 0-2
ymax = ny+2
compx = np.linspace(0,nx-1,nx)
compy = np.linspace(0,ny-1,ny)
yv,xv = np.meshgrid(compy,compx)
plot_mesh(xv,yv, ax)
plt.axis([xmin,xmax,ymin,ymax])
plt.fill([xmin,xmin,xmax,xmax],[ymin,ymax,ymax,ymin],'gray')
ax.set_xlabel('I')
ax.set_ylabel('J')
ax.title.set_text('Computational Space: '+str(nx)+' x '+str(ny))

fname = './export/mont-grid.png'
plt.savefig(fname,dpi=100)
print("Created plot: "+fname) 




#- Export Grid vertices to file
out_file = './output/mont_grid.xy'
nx,ny = xi.shape
header = '## '+str(nx)+' x '+str(ny)
fid = open(out_file,'w')
fid.write('%s\n' % header)
for j in range(ny):
    for i in range(nx):
        fid.write('%.2f %.2f\n' % (xi[i,j],yi[i,j]))
fid.close()
print("Created grid corner nodes: "+out_file) 

#- Export Grid center points to file
out_file = './output/mont_grid_CE.xy'
nx,ny = xi_rho.shape
header = '## '+str(nx)+' x '+str(ny)
fid = open(out_file,'w')
fid.write('%s\n' % header)
for j in range(ny):
    for i in range(nx):
        fid.write('%.2f %.2f\n' % (xi_rho[i,j],yi_rho[i,j]))
fid.close()
print("Created grid center nodes: "+out_file) 


