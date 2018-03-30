'''
exec(open('./scripts/create_initial_conds.py').read())
% Pre-processing tool for La Jolla Region
% *Description:* Step by step for setting the input files required by GCCOM 
% after the mesh is created.
% 
% INPUT: 
% 
%      Grids fies:  Grid.dat and ProbSize.dat
% 
%  OUTPUT:
% 
%     Netcdf input file and Boundary Conditions
% 
%  IMPORTANT: This file needs to be edited manually depending to set up the 
% proper conditions for each problem. Please read carfully each step.
% 
% Info from ROMS
% The variables being extracted are: 
%    + temperature (temp), 
%    + salinity (salt), 
%    + u (zonal current speed), 
%    + v (meridional current speed)
%    + ssh (sea surface height)
% - They are all on the same 3D (lon/lat/depth) grid, except that ssh is defined at 
%   the surface only so there is no depth dimension.   
% - The lon/lat/depth variables are also extracted in the script.   
%     o There are 12 depth levels from the surface (= 0 meters) to a depth of 400 meters
%       (see depthdata for all 12 depth values).    
% - Concerning lon/lat, this script extracts the full domain (351x391 points, 3km resolution)
% - Also, note that for the ocean variables: 
%     o land values are given the 'missing' value of -9999
'''

import os
import sys
import shutil
import math
import netCDF4
import numpy as np
import datetime
#from datetime import *
import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append("./scripts")
from gccom_utils import *
from gccom_plot import *
from gis_utils import *

#-----------------------------------------------------------------------------------------------------
#- Function to convert ROMS time to mat
#-     INPUT: ncid = netcdf dataset resource
#      OUTPUT: dates,num_times,xlabel_time,time_sec
#-----------------------------------------------------------------------------------------------------
def ROMStime2mat(ncid):
    import datetime as dt
    nc = ncid#netCDF4.Dataset(ncid)
#    nc_vars = [var for var in nc.variables]       # list of nc variables
    times = nc.variables['time']
    timeunits = nc.variables['time'].units        #- 'hour since 2015-09-15 03:00:00'
    tmp = timeunits.split()
    timebase = tmp[2].split('-')             #- (YYYY,MM,DD)
    timebase2 = tmp[3].split(':')
    hour_origin = int(timebase2[0])
    time_origin = dt.datetime(int(timebase[0]),int(timebase[1]),int(timebase[2]),hour_origin,0,0);
    roms_times = []          #- time object list
    roms_strings = []        #- time string list
    xlabel_times = []
    delta_hour = dt.timedelta(hours=+1)    #- number of secs in hour
    roms_times.append(time_origin)
    roms_strings.append(time_origin.strftime('%Y-%m-%d %H:%M:%S'))
    for i in range(len(times)):
        delta_hour = dt.timedelta(hours=int(times[i]))
        new_time = time_origin + delta_hour
        roms_times.append(new_time)
        roms_strings.append(new_time.strftime('%Y-%m-%d %H:%M:%S'))
    start_time = time_origin
    final_time = roms_times[len(times)-1]
    num_times  = len(times)
    elapsed_time = final_time - start_time
    elapsed_sec = elapsed_time.total_seconds()
    return roms_times,roms_strings,elapsed_sec



#- Begin preprocessing script
#- Define source files, will default to lookd for 'Grid.dat' & 'ProbSize.dat'
#- User can input 'Grid_sdbay.dat' & script will automatically find 'ProbSize_sdbay.dat'
input_args = sys.argv
if (len(sys.argv) < 2):
    GridFileName = './LAJOLLA_40x24x14/Grid_40x24x14.dat'
    ProbSizeFileName='./LAJOLLA_40x24x14/ProbSize_40x24x14.dat'
else:
    GridFileName = input_args[1]       #- Grid_sdbay.dat
    ProbSizeFileName=GridFileName.replace('Grid','ProbSize')    #- ProbSize_sdbay.dat


if (os.path.exists(GridFileName) == False):
    print("Cannot find grid file: ./"+GridFileName)
    sys.exit()
if (os.path.exists(ProbSizeFileName) == False):
    print("Cannot find grid file: ./"+ProbSizeFileName)
    sys.exit()

#- Get ROMS netcdf dataset
#ROMSFileName = 'http://west.rssoffice.com:8080/thredds/dodsC/roms/CA3km-forecast/CA/ca_subCA_fcstair8_2015091503.nc'
ROMSFileName = 'http://west.rssoffice.com:8080/thredds/dodsC/roms/CA3km-forecast/CA/ca_subCA_fcst_2014120503.nc'
fname = ROMSFileName
print('Opening ROMS NETCDF file')
nc = netCDF4.Dataset(fname)


#- Grab start/end date from ROMS dataset
roms_times,roms_strings,elapsed_sec = ROMStime2mat(nc)

print('\nROMS start time %s '%roms_strings[0])
print('ROMS final time %s '%roms_strings[len(roms_strings)-1])
print('ROMS elapsed time in seconds: %d \n'%elapsed_sec)

#- 
#- Compute time since '1980-01-01 00:00:00'
start_time = roms_times[0]
ref_time = datetime.datetime(1980,1,1,0,0,0)
delta_time = start_time-ref_time              #- approx = 13041.041666666666 ;
#timer = (roms_strings[0],roms_strings[len(roms_strings)-1])
timer = roms_times[0].toordinal()


#-- Set the Bcs directory
BCDir = 'BCS_ROMS'
#-- Set the BC time filename
OpenBCtimesfile = BCDir+'/roms_bry_timesLJ.dat'
#-- Set Sub Directories for BCS
newdiru='u_roms_bdy_LJ'
newdirv='v_roms_bdy_LJ'
newdirw='w_roms_bdy_LJ'
newdirp='p_roms_bdy_LJ'

#- Check mesh to ensure it meets stability conditions
#checkmesh(GridFileName,ProbSizeFileName)

#- Read Grid and return as vectors IMax*JMax*KMax long
X,Y,Z,IMax,JMax,KMax = read_grid_ascii(GridFileName,ProbSizeFileName)

#- set output file name --> gcom_ini_SD_100x25x10_ROMSBC.nc
ini_file = 'gcom_ini_LJ_'+str(IMax)+'x'+str(JMax)+'x'+str(KMax)+'_ROMSBC.nc';

print('Creating NETCDF file: %s '% ini_file)
#- Initializing the gccom variables
gccom_write_netcdf(ini_file,timer,IMax,JMax,KMax,X,Y,Z )

#- Find Boundaries of GCCOM
#- Compute staggered grids
x_u,y_u,z_u = staggered_grid('u',X,Y,Z,IMax,JMax,KMax)
x_v,y_v,z_v = staggered_grid('v',X,Y,Z,IMax,JMax,KMax)
x_w,y_w,z_w = staggered_grid('w',X,Y,Z,IMax,JMax,KMax)

#- Find the boundary conditions at each grid
BCx_u,BCy_u,BCz_u = findBC(x_u,y_u,z_u)
BCx_v,BCy_v,BCz_v = findBC(x_v,y_v,z_v)

#- Transform from utm (meters) to geographic (degrees)
BCx_u,BCy_u=utm2ll(BCx_u,BCy_u,11)
BCx_v,BCy_v=utm2ll(BCx_v,BCy_v,11)

#- Plotting here ???
#lon,lat=utm2ll(np.array(x), np.array(y),11); ## for visualization porpose only

#- Read info from ROMS netcdf file
Rlat = nc.variables['lat'][:]
Rlon = nc.variables['lon'][:]
Rlon = Rlon[:] - 360
Rdepth = nc.variables['depth'][:]*(-1.0)   #- depths are positive in ROMS data
RIM,RJM,RKM = len(Rlon),len(Rlat),len(Rdepth)

#- Create 3-D matrix from ROMS vectors
[RXq,RYq,RZq] = np.meshgrid(Rlon,Rlat,Rdepth)

#-Find the region to extract the data from ROMS
GRegLon= [min(BCx_u[:]), max(BCx_u[:])]
GRegLat= [min(BCy_v[:]), max(BCy_v[:])]

#- find the closest point to the GCCOM Bcs  
#- Follow up xx,='np.where' with comma to remove tuple and 
#- assign xx to array
ind_lon_min, = np.where(Rlon<=GRegLon[0])
ind_lon_min = ind_lon_min[-1]
ind_lon_max, = np.where(Rlon<=GRegLon[1])
ind_lon_max = ind_lon_max[-1]
ind_lat_min, = np.where(Rlat<=GRegLat[0])
ind_lat_min = ind_lat_min[-1]
ind_lat_max, = np.where(Rlat<=GRegLat[1])
ind_lat_max = ind_lat_max[-1]

#- create a region where data will be extracted
#Regionlat=ind_lat_min-10:ind_lat_max+10
#Regionlon=ind_lon_min-5:ind_lon_max+5
 
#- Extract ROMS Latitude and Longitude in a buffered
#- region: LAT=+/-10 cells; LON +/- 5 cells 
#- Note: these are only vectors
dx = 5
Rlon_ind = list(np.linspace(ind_lon_min-dx,ind_lon_max+dx,dx*2+1))
Rlon_ind = [int(i) for i in Rlon_ind]
dy = 10 
Rlat_ind = list(np.linspace(ind_lat_min-dy,ind_lat_max+dy,dy*2+1))
Rlat_ind = [int(i) for i in Rlat_ind]
RSublon =Rlon[Rlon_ind]
RSublat =Rlat[Rlat_ind]

#- How long you want to run the simulation
#- The netcdf file has 73 files, frequency every hour
bry_times = 3600*np.arange(0,73)         #- setting for 72 hours.
Nrec = len(bry_times)

#- Initialize all the variables for BCS 
ubcw,ubce,ubcn,ubcs,pbct = gccom_ini_BC('u',Nrec,IMax,JMax,KMax)
vbcw,vbce,vbcn,vbcs,pbct  = gccom_ini_BC('v',Nrec,IMax,JMax,KMax)
wbcw,wbce,wbcn,wbcs,pbct = gccom_ini_BC('w',Nrec,IMax,JMax,KMax)
#- Note pbcw,pbce,pbcn,pbcs are dummies variables only BC at the top
#-  Represents the pressure at the top of water column
pbcw,pbce,pbcn,pbcs,pbct = gccom_ini_BC('p',Nrec,IMax,JMax,KMax)


#- Get ROMS Variables and subset the region using lat/lon bounds
#- Note: Python reads in NETCDF opposite order: 
#-   + Shape of ROMS 'u' = 351x391x14x73
#-   + Shape of subset ROMS 'u' = 11x21x14x73
#- However, will be read in Python as:
#-   + Shape of ROMS 'u' = 73x14x391x351
#-   + Shape of subset ROMS 'u' = 73x14x21x11
print('Reading in U and V variables from NETCDF')
u = nc.variables['u'][:,:,Rlat_ind,Rlon_ind]    #- u.shape = (73, 14, 21, 11)   
v = nc.variables['v'][:,:,Rlat_ind,Rlon_ind]    #- v.shape = (73, 14, 21, 11) 

#- Close ROMS file_id
nc.close()

#- Plot u-values
#fig = plt.figure()
#plt.pcolor(RSublon,RSublat,u[0,0,:,:])

#- Interpolate these variables

#- Loop through all time records (N=Nrec)
for rec in range(Nrec):
#for rec in range(2):
    print("Time Record: "+str(rec))
    #- Interpolate variables
    uone= np.squeeze(u[rec,:,:,:])  #- uone.shape is (w,v,u) 
    uone = uone.swapaxes(0,2)       #- Will reorder (u,v,w)
    vone= np.squeeze(v[rec,:,:,:])  #- vone.shape is (w,v,v)
    vone = vone.swapaxes(0,2)       #- will reorder (u,v,w)
    print('    Interpolating values ...')
    u_interp,u_real = interp_roms2gccom(RSublon,RSublat,Rdepth,uone,BCx_u,BCy_u,BCz_u)
    v_interp,v_real = interp_roms2gccom(RSublon,RSublat,Rdepth,vone,BCx_v,BCy_v,BCz_v)
    #- Manually: U-Mode forcing in U_direction
    u_in = np.zeros([JMax-1,KMax-1])
    v_in = np.zeros([JMax,KMax-1])
   
    #- BC for for U
    l=0
    for k in range(KMax-1):
        for j in range(JMax-2,-1,-1):  #- range is non-inclusive of endpt, but inclusive of start
            u_in[j,k]= u_interp[l]
            l=l+1


    #- BC for for U
    l=0
    for k in range(KMax-1):
        for j in range(JMax-1,-1,-1):  #- range is non-inclusive of endpt, but inclusive of start
            v_in[j,k]= v_interp[l]
            l=l+1

    # Settings ghost points for bcs for U. (Do not touch)
    print('    Setting Ghost Points ...')
    u_inall = gccom_ghost_BC('ubcw',u_in, IMax,JMax,KMax)
    v_inall = gccom_ghost_BC('vbcw',v_in, IMax,JMax,KMax);

    #- Simulating a kind of tidal wave. ---> NOT SURE ABOUT THIS
    ubcw[:,:,:,rec] = u_inall[:,:,:]
    vbcw[:,:,:,rec] = v_inall[:,:,:]


#- Make a directory called './BCS_ROMS'
if(os.path.exists(BCDir)):
    shutil.rmtree(BCDir)
if(not os.path.exists(BCDir)):
    os.makedirs(BCDir)

#- Call function to write the BCs to netcdf file
gccom_write_BC(BCDir,'u',bry_times,newdiru,ubcw,ubce,ubcs,ubcn,IMax, JMax,KMax)
gccom_write_BC(BCDir,'v',bry_times,newdirv,vbcw,vbce,vbcs,vbcn,IMax, JMax,KMax)
gccom_write_BC(BCDir,'w',bry_times,newdirw,wbcw,wbce,wbcs,wbcn,IMax, JMax,KMax)
gccom_write_BC(BCDir,'p',bry_times,newdirp,pbcw,pbce,pbcs,pbct,IMax, JMax,KMax) 
#- Note pbcw,pbce,pbcn,pbcs are dummies variables only BC at the top is uses.

#- Now write times to OpenBCtimesfile = 'roms_bry_timesSD.dat'
Nrec = len(bry_times); 
fileID = open(OpenBCtimesfile,'w')
for nt in range(Nrec):
    fileID.write('%12.16f\n' % bry_times[nt])
fileID.close()

print('*****************************************')
print('*** GCCOM inputs files ready!!***********')
print('*****************************************')


#- Now copy to final directory
grid_file = os.path.dirname(GridFileName)    #-  '../SDBAY_99x19x16/Grid_99x19x16.dat'
gccom_dir = os.path.relpath(grid_file)		#- '../SDBAY_99x19x16'
shutil.copytree(BCDir,gccom_dir+'/BCS_ROMS')
shutil.copy(ini_file,gccom_dir)


