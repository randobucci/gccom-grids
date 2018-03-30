#-----------------------------------------------------------------------------------
#- function check_mesh
#-----------------------------------------------------------------------------------

import numpy as np
import sys
import shutil
import os
import netCDF4

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

#-------------------------------------------------------------------
#- FUNCTION findBC
#- [BCx_u,BCy_u,BCz_u]=findBC(x_u, y_u, z_u);
def findBC(x_u, y_u, z_u):
    end = len(x_u[0,:,0]) - 1
    x_cells = x_u[0,end::-1,:]
    y_cells = y_u[0,end::-1,:]
    z_cells = z_u[0,end::-1,:]
    NI,NK = x_cells.shape
    xBC = np.reshape(x_cells,NI*NK,1)
    yBC = np.reshape(y_cells,NI*NK,1)
    zBC = np.reshape(z_cells,NI*NK,1)
    return xBC,yBC,zBC
#------------------------------------------------------------------


#-------------------------------------------------------------------
#- Input 3-Vars and return centers as vectors 
#-------------------------------------------------------------------
def calc_center(x,y,z):
    IM,JM,KM = x.shape 
    Xc = np.empty((IM,JM,KM))
    Yc = np.empty((IM,JM,KM))
    Zc = np.empty((IM,JM,KM))
    theI = 0
    for k in range(KM-1):
        for j in range(JM-1):
            for i in range(IM-1):
                theI = theI + 1
                Xc[i,j,k]=(x[i,j,k]+x[i+1,j,k]+x[i,j+1,k]+x[i+1,j+1,k]+ \
                       x[i,j,k+1]+x[i+1,j,k+1]+x[i,j+1,k+1]+x[i+1,j+1,k+1]);
                Yc[i,j,k]=(y[i,j,k]+y[i+1,j,k]+y[i,j+1,k]+y[i+1,j+1,k]+ \
                       y[i,j,k+1]+y[i+1,j,k+1]+y[i,j+1,k+1]+y[i+1,j+1,k+1]);
                Zc[i,j,k]=(z[i,j,k]+z[i+1,j,k]+z[i,j+1,k]+z[i+1,j+1,k]+ \
                       z[i,j,k+1]+z[i+1,j,k+1]+z[i,j+1,k+1]+z[i+1,j+1,k+1]);
    Xc = Xc/8.0
    Yc = Yc/8.0
    Zc = Zc/8.0
    return Xc,Yc,Zc
#-------------------------------------------------------------------
 


#-------------------------------------------------------------------
def check_mesh(grid_file,prob_size_file):

    X,Y,Z,IMax,JMax,KMax = read_grid_ascii(grid_file,prob_size_file)

    #- Calculate centers
    Xc,Yc,Zc = calc_center(X,Y,Z)

    #- Get indices of bottom k-level
    l = 1
#-------------------------------------------------------------------

#--STAGGERED_GRID-----------------------------------------------------------------
#- From P.Choboter matlab code StaggeredGrid.m
'''
function [x_out,y_out,z_out]=StaggeredGrid(g,x,y,z,IM,JM,KM)
%  This code takes x,y,z on the cell corners, and returns the 
%  x, y, and z coordinates on the other cell faces and center.  
%  x,y,z := coordinates on the cell corners, 3D arrays
%  g := what grid to output: 'u', 'v', 'w', or 'r' (r=rho,pressure,T,S)
%  IM,JM,KM := size of the grid.
%  Written by P. Choboter
%
%  Usage:  
%  [x_u, y_u, z_u] = StaggeredGrid('u',Z,y,z,IM,JM,KM);
%  [x_v, y_v, z_v] = StaggeredGrid('v',x,y,z,IM,JM,KM);
%  [x_r, y_r, z_r] = StaggeredGrid('r',x,y,z,IM,JM,KM);
'''
def staggered_grid(g,x,y,z,IMax,JMax,KMax):
    #- Adjust for python indexing starting at i=0, not i=1
    IM = IMax
    JM = JMax
    KM = KMax
    imm=IMax-1
    jmm=JMax-1
    kmm=KMax-1
    if (g == 'u'):
        x_out = 0.25*(x[0:IM,0:jmm,0:kmm]+x[0:IM,1:JM,0:kmm]+ \
               x[0:IM,0:jmm,1:KM]+x[0:IM,1:JM,1:KM])
        y_out = 0.25*(y[0:IM,0:jmm,0:kmm]+y[0:IM,1:JM,0:kmm]+ \
               y[0:IM,0:jmm,1:KM]+y[0:IM,1:JM,1:KM])
        z_out = 0.25*(z[0:IM,0:jmm,0:kmm]+z[0:IM,1:JM,0:kmm]+ \
               z[0:IM,0:jmm,1:KM]+z[0:IM,1:JM,1:KM])
    elif (g == 'v'):
        x_out = 0.25*(x[0:imm,0:JM,0:kmm]+x[1:IM,0:JM,0:kmm]+ \
               x[0:imm,0:JM,1:KM]+x[1:IM,0:JM,1:KM])
        y_out = 0.25*(y[0:imm,0:JM,0:kmm]+y[1:IM,0:JM,0:kmm]+ \
               y[0:imm,0:JM,1:KM]+y[1:IM,0:JM,1:KM])
        z_out = 0.25*(z[0:imm,0:JM,0:kmm]+z[1:IM,0:JM,0:kmm]+ \
               z[0:imm,0:JM,1:KM]+z[1:IM,0:JM,1:KM])
    elif (g == 'w'):
        x_out = 0.25*(x[0:imm,0:jmm,0:KM]+x[1:IM,0:jmm,0:KM]+ \
               x[0:imm,1:JM,0:KM]+x[1:IM,1:JM,0:KM])
        y_out = 0.25*(y[0:imm,0:jmm,0:KM]+y[1:IM,0:jmm,0:KM]+ \
               y[0:imm,1:JM,0:KM]+y[1:IM,1:JM,0:KM])
        z_out = 0.25*(z[0:imm,0:jmm,0:KM]+z[1:IM,0:jmm,0:KM]+ \
               z[0:imm,1:JM,0:KM]+z[1:IM,1:JM,0:KM])
    elif ( (g == 'r') | (g=='T') | (g=='S') | (g=='D') ):
        x_out = 0.125*( x[0:imm,0:jmm,0:kmm]+x[1:IM,0:jmm,0:kmm]+ \
               x[0:imm,1:JM,0:kmm]+x[1:IM,1:JM,0:kmm]+ \
               x[0:imm,0:jmm,1:KM]+x[1:IM,0:jmm,1:KM]+ \
               x[0:imm,1:JM,1:KM] + x[1:IM,1:JM,1:KM] )
        y_out = 0.125*( y[0:imm,0:jmm,0:kmm]+y[1:IM,0:jmm,0:kmm]+ \
               y[0:imm,1:JM,0:kmm]+y[1:IM,1:JM,0:kmm]+ \
               y[0:imm,0:jmm,1:KM]+y[1:IM,0:jmm,1:KM]+ \
               y[0:imm,1:JM,1:KM] + y[1:IM,1:JM,1:KM] )
        z_out = 0.125*( z[0:imm,0:jmm,0:kmm]+z[1:IM,0:jmm,0:kmm]+ \
               z[0:imm,1:JM,0:kmm]+z[1:IM,1:JM,0:kmm]+ \
               z[0:imm,0:jmm,1:KM]+z[1:IM,0:jmm,1:KM]+ \
               z[0:imm,1:JM,1:KM] + z[1:IM,1:JM,1:KM] )
    else:
        print('StaggeredGrid error! First input should be u, v, w, r, T, S, or D, in quotes.')
        sys.exit()
    return x_out,y_out,z_out
#-------------------------------------------------------------------

#----ZB_STAGGERED_GRID------------------------------------------------
#- From P.Choboter matlab code zb_StaggeredGrid.m
'''
%  This code takes z on the cell corners, and returns the 
%  value of z at the bottom of the domain in the various grids.
%  z := 3D z-coordinates on the cell corners
%  IM,JM,KM := size of the grid.  KM not used. 
%  zb   := 2D array below the cell corners
%  zb_u := 2D array below the u-grid
%  zb_v := 2D array below the v-grid
%  zb_r := 2D array below the cell centers
%  Written by P. Choboter
%
%  Usage:  
%  [zb,zb_u,zb_v,zb_r]=zb_Staggered(z,IM,JM);
'''
#-------------------------------------------------------------------
def zb_staggered_grid(z,IM,JM):
    imm = IM-1
    jmm = JM-1
    zb = np.squeeze(z[:,:,0])
    zb_u = 0.5*( zb[0:IM,0:jmm] + zb[0:IM,1:JM] )
    zb_v = 0.5*( zb[0:imm,0:JM] + zb[1:IM,0:JM] )
    zb_r =0.25*( zb[0:imm,0:jmm]+ zb[1:IM,0:jmm] + \
                 zb[0:imm,1:JM] + zb[1:IM,1:JM] )
    return zb, zb_u, zb_v, zb_r


#-------------------------------------------------------------------
def gccom_write_netcdf(ncfile_name,timer,IM,JM,KM,x,y,z):
    #nc_id = netCDF4.Dataset(ncfile_name, 'w', format='NETCDF4_CLASSIC')
    nc_id = netCDF4.Dataset(ncfile_name, 'w')

    #- Initialize a bunch of variables
    uic = np.empty((IM,JM-1,KM-1))
    vic = np.empty((IM-1,JM,KM-1))
    wic = np.empty((IM-1,JM-1,KM))
    Tic = np.empty((IM-1,JM-1,KM-1))
    Sic = np.empty((IM-1,JM-1,KM-1))
    Dic = np.empty((IM-1,JM-1,KM-1))
    Pic = np.empty((IM-1,JM-1,KM-1))
    Vortex_x = np.empty((IM,JM,KM))
    Vortex_y = np.empty((IM,JM,KM))
    Vortex_z = np.empty((IM,JM,KM))

    #- Compute staggered grids
    x_u,y_u,z_u = staggered_grid('u',x,y,z,IM,JM,KM)
    x_v,y_v,z_v = staggered_grid('v',x,y,z,IM,JM,KM)
    x_w,y_w,z_w = staggered_grid('w',x,y,z,IM,JM,KM)
    x_r,y_r,z_r = staggered_grid('r',x,y,z,IM,JM,KM)
    zb,zb_u,zb_v,zb_r = zb_staggered_grid(z,IM,JM)

    #- Initial conditions for temperature
    for k in range(KM-1):
        for j in range(JM-1):
            for i in range(IM-1):
                t1 = 12-0.6*(np.tanh((-1*z_r[i,j,k]-27)/5))
                t2 = -0.8-0.9*(np.tanh((-1*z_r[i,j,k]-70)/15))
                Tic[i,j,k]=t1+t2

 
   
    #- Declare some helpful variables 
    nxp1,nyp1,nzp1 = IM,JM,KM
    nx,ny,nz = IM-1,JM-1,KM-1
    nc_dims = ('time','nx','ny','nz','nxp1','nyp1','nzp1')
#    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    #- Begin writing to netcdf file
    nc_id.description = "Computational Science Research Center (CSRC) "
    nc_id.description += "Initialization File For General Curvilinear "
    nc_id.description += "Ocean Model (GCCOM)"
    #- Create Dimensions
    nc_id.createDimension('time',None) 
    nc_id.createDimension('nx',nx) 
    nc_id.createDimension('ny',ny) 
    nc_id.createDimension('nz',nz) 
    nc_id.createDimension('nxp1',nxp1) 
    nc_id.createDimension('nyp1',nyp1) 
    nc_id.createDimension('nzp1',nzp1) 
    #- Create Variables
    time=nc_id.createVariable('time', 'f8',('time',))
    ulon=nc_id.createVariable('ulon', 'f8',('nz','ny','nxp1'))
    ulat=nc_id.createVariable('ulat', 'f8',('nz','ny','nxp1'))
    ulev=nc_id.createVariable('ulev', 'f8',('nz','ny','nxp1'))
    vlon=nc_id.createVariable('vlon', 'f8',('nz','nyp1','nx'))
    vlat=nc_id.createVariable('vlat', 'f8',('nz','nyp1','nx'))
    vlev=nc_id.createVariable('vlev', 'f8',('nz','nyp1','nx'))
    wlon=nc_id.createVariable('wlon', 'f8',('nzp1','ny','nx'))
    wlat=nc_id.createVariable('wlat', 'f8',('nzp1','ny','nx'))
    wlev=nc_id.createVariable('wlev', 'f8',('nzp1','ny','nx'))
    lon=nc_id.createVariable('lon', 'f8',('nz','ny','nx'))
    lat=nc_id.createVariable('lat', 'f8',('nz','ny','nx'))
    lev=nc_id.createVariable('lev', 'f8',('nz','ny','nx'))
    u=nc_id.createVariable('u', 'f8',('nz','ny','nxp1'))
    v=nc_id.createVariable('v', 'f8',('nz','nyp1','nx'))
    w=nc_id.createVariable('w', 'f8',('nzp1','ny','nx'))
    p=nc_id.createVariable('p', 'f8',('nz','ny','nx'))
    T=nc_id.createVariable('T', 'f8',('nz','ny','nx'))
    S=nc_id.createVariable('S', 'f8',('nz','ny','nx'))
    D=nc_id.createVariable('D', 'f8',('nz','ny','nx'))
    vortx=nc_id.createVariable('Vortex_x', 'f8',('nzp1','nyp1','nxp1'))
    vorty=nc_id.createVariable('Vortex_y', 'f8',('nzp1','nyp1','nxp1'))
    vortz=nc_id.createVariable('Vortex_z', 'f8',('nzp1','nyp1','nxp1'))

    #- Create Attributes
    time.setncatts({'units':'days since 1980-01-01 00:00:00',\
                    'long_name':'time','calendar':'gregorian'})
    ulon.setncatts({'long_name':'longitude','units':'degrees_east','axis':'X',\
                    'comment':'longitude of the u component of velocity'})
    ulat.setncatts({'long_name':'latitude','units':'degrees_north','axis':'Y',\
                    'comment':'latitude of the u component of velocity'})
    ulev.setncatts({'long_name':'depth','units':'meters','axis':'Z','positive':'up',\
                    'comment':'depth of the u component of velocity'})

    vlon.setncatts({'long_name':'longitude','units':'degrees_east','axis':'X',\
                    'comment':'longitude of the v component of velocity'})
    vlat.setncatts({'long_name':'latitude','units':'degrees_north','axis':'Y',\
                    'comment':'latitude of the v component of velocity'})
    vlev.setncatts({'long_name':'depth','units':'meters','axis':'Z','positive':'up',\
                    'comment':'depth of the v component of velocity'})

    wlon.setncatts({'long_name':'longitude','units':'degrees_east','axis':'X',\
                    'comment':'longitude of the w component of velocity'})
    wlat.setncatts({'long_name':'latitude','units':'degrees_north','axis':'Y',\
                    'comment':'latitude of the w component of velocity'})
    wlev.setncatts({'long_name':'depth','units':'meters','axis':'Z','positive':'up',\
                    'comment':'depth of the w component of velocity'})

    lon.setncatts({'long_name':'longitude','units':'degrees_east','axis':'X',\
                    'comment':'longitude of the grid cell center'})
    lat.setncatts({'long_name':'latitude','units':'degrees_north','axis':'Y',\
                    'comment':'latitude of the grid cell center'})
    lev.setncatts({'long_name':'depth','units':'meters','axis':'Z','positive':'up',\
                    'comment':'depth of the grid cell center'})

    u.setncatts({'units':'m/s','long_name':'u component of the velocity',\
                 'coordinates':'ulon ulat ulev time'})
    v.setncatts({'units':'m/s','long_name':'v component of the velocity',\
                 'coordinates':'vlon vlat vlev time'})
    w.setncatts({'units':'m/s','long_name':'w component of the velocity',\
                 'coordinates':'wlon wlat wlev time'})
    p.setncatts({'units':'bar','long_name':'non-hydrostatic pressure',\
                 'coordinates':'lon lat lev time'})
    T.setncatts({'units':'degrees C','long_name':'potential temperature',\
                 'coordinates':'lon lat lev time'})
    S.setncatts({'units':'practical salinity units','long_name':'salinity',\
                 'coordinates':'lon lat lev time'})
    D.setncatts({'units':'g/cm3','long_name':'density',\
                 'coordinates':'lon lat lev time'})
    
    #- First write start time record
    nc_id.variables['time'][0] = timer

    #- Define a list of vars to iterate through and initialize  
    var_list = [var for var in nc_id.variables]
    var_list.pop(0)   #- remove time variable
    data_list = [x_u,y_u,z_u,x_v,y_v,z_v,x_w,y_w,z_w,x_r,y_r,z_r,
                 uic,vic,wic,Pic,Dic,Tic,Sic,Vortex_x,Vortex_y,Vortex_z]
    #- Loop through the variable list and set to vals in corres data_list 
    for i in range(len(var_list)):
        #- Need to swap the axes of these data i.e. --> (w,v,u) to (u,v,w)
        out_data = data_list[i].swapaxes(0,2)
        #nc_id.variables[var_list[i]][:] = data_list[i]
        nc_id.variables[var_list[i]][:] = out_data
    
    #- Close the netcdf file
    nc_id.close() 
#-------------------------------------------------------------------

#----GCCOM_INI_BC------------------------------------------------------
#- From GCCOM matlab code gccom_ini_BC.m
'''
   This code takes initializes all of the BC variables 
   INPUTS:
          z := string key letter: 'u','v','w','p'
          Nrec := Number of vlues
          IMax,JMax,KMax := size of the grid.
   OUTPUTS:
          bcw,bce,bcn,bcs,pbct
'''
#---------------------------------------------------------------------
def gccom_ini_BC(var,Nrec,IMax,JMax,KMax):
    if(var == 'u'):
        bcw = np.zeros((3, JMax+3, KMax+3, Nrec))
        bce = np.zeros((3, JMax+3, KMax+3, Nrec))
        bcn = np.zeros((IMax+4, 2, KMax+3, Nrec))
        bcs = np.zeros((IMax+4, 2, KMax+3, Nrec))
        pbct = None
    elif(var == 'v'):
        bcw = np.zeros((2, JMax+4, KMax+3, Nrec))
        bce = np.zeros((2, JMax+4, KMax+3, Nrec))
        bcn = np.zeros((IMax+3, 3, KMax+3, Nrec))
        bcs = np.zeros((IMax+3, 3, KMax+3, Nrec))
        pbct = None
    elif(var == 'w'):
        bcw = np.zeros((2, JMax+3, KMax+4, Nrec))
        bce = np.zeros((2, JMax+3, KMax+4, Nrec))
        bcn = np.zeros((IMax+3, 2, KMax+4, Nrec))
        bcs = np.zeros((IMax+3, 2, KMax+4, Nrec))
        pbct = None
    elif(var == 'p'):
        #- pbct := Pressure at top of water col
        pbct = np.zeros((IMax+3, JMax+3, Nrec))
        bcw = np.zeros((2, JMax+3, KMax+4, Nrec))
        bce = np.zeros((2, JMax+3, KMax+4, Nrec))
        bcn = np.zeros((IMax+3, 2, KMax+4, Nrec))
        bcs = np.zeros((IMax+3, 2, KMax+4, Nrec))
    return bcw,bce,bcn,bcs,pbct


#----GCCOM_WRITE_BC------------------------------------------------------
#- From GCCOM matlab code gccom_write_BC.m
'''
   This code takes writes all of the BC variables 
   INPUTS:
          z := string key letter: 'u','v','w','p'
          Nrec := Number of vlues
          IMax,JMax,KMax := size of the grid.
   OUTPUTS:
          bcw,bce,bcn,bcs,pbct
'''
#---------------------------------------------------------------------
def gccom_write_BC(bc_dir,variable,bry_times,newdirvar,bcw,bce,bcs,bcn,IMax, JMax,KMax):
    Nrec = len(bry_times)
    #- Switch variables
    if(variable == 'u'):
        newdir = bc_dir+'/'+newdirvar      #-  ./BCS_ROMS/u_roms_bdy_SD
        if(not os.path.exists(newdir)):
            os.makedirs(newdir)
        #- Loop through each times
        for rec in (range(Nrec)):
            if (len(str(rec)) == 1): IDno='00000'+str(rec)
            if (len(str(rec)) == 2): IDno='0000'+str(rec)
            if (len(str(rec)) == 3): IDno='000'+str(rec)
            filename = newdir+'/u_bry_'+IDno+'.dat'
            fileID = open(filename,'w')
            #-  Array size reminder
            #- ubcw(3, JMax+3, KMax+3, Nrec);
            #- ubce(3, JMax+3, KMax+3, Nrec);
            #- ubcn(IMax+4, 2, KMax+3, Nrec);
            #- ubcs(IMax+4, 2, KMax+3, Nrec);
            for k in range(KMax+3):
                for j in range(JMax+3):
                    for i in range(3):
                        fileID.write('%26.16f\n' % bcw[i,j,k,rec])
            for k in range(KMax+3):
                for j in range(JMax+3):
                    for i in range(3):
                        fileID.write('%26.16f\n' % bce[i,j,k,rec])
            for k in range(KMax+3):
                for j in range(2):
                    for i in range(IMax+4):
                        fileID.write('%26.16f\n' % bcs[i,j,k,rec])
            for k in range(KMax+3):
                for j in range(2):
                    for i in range(IMax+4):
                        fileID.write('%26.16f\n' % bcn[i,j,k,rec])
            fileID.close()
            print('*****************************************')
            print('***BC for u are ready:  '+filename)
            print('*****************************************')
    elif(variable == 'v'):
        newdir = bc_dir+'/'+newdirvar      #-  ./BCS_ROMS/u_roms_bdy_SD
        if(not os.path.exists(newdir)):
            os.makedirs(newdir)
        #- Loop through each times
        for rec in (range(Nrec)):
            if (len(str(rec)) == 1): IDno='00000'+str(rec)
            if (len(str(rec)) == 2): IDno='0000'+str(rec)
            if (len(str(rec)) == 3): IDno='000'+str(rec)
            filename = newdir+'/v_bry_'+IDno+'.dat'
            fileID = open(filename,'w')
            #-  Array size reminder
            #- vbcw(2, JMax+4, KMax+3, Nrec);
            #- vbce(2, JMax+4, KMax+3, Nrec);
            #- vbcn(IMax+3, 3, KMax+3, Nrec);
            #- vbcs(IMax+3, 3, KMax+3, Nrec);
            for k in range(KMax+3):
                for j in range(JMax+4):
                    for i in range(2):
                        fileID.write('%26.16f\n' % bcw[i,j,k,rec])
            for k in range(KMax+3):
                for j in range(JMax+4):
                    for i in range(2):
                        fileID.write('%26.16f\n' % bce[i,j,k,rec])
            for k in range(KMax+3):
                for j in range(3):
                    for i in range(IMax+3):
                        fileID.write('%26.16f\n' % bcs[i,j,k,rec])
            for k in range(KMax+3):
                for j in range(3):
                    for i in range(IMax+3):
                        fileID.write('%26.16f\n' % bcn[i,j,k,rec])
            fileID.close()
            print('*****************************************')
            print('***BC for v are ready:  '+filename)
            print('*****************************************')
    elif(variable == 'w'):
        newdir = bc_dir+'/'+newdirvar      #-  ./BCS_ROMS/u_roms_bdy_SD
        if(not os.path.exists(newdir)):
            os.makedirs(newdir)
        #- Loop through each times
        for rec in (range(Nrec)):
            if (len(str(rec)) == 1): IDno='00000'+str(rec)
            if (len(str(rec)) == 2): IDno='0000'+str(rec)
            if (len(str(rec)) == 3): IDno='000'+str(rec)
            filename = newdir+'/w_bry_'+IDno+'.dat'
            fileID = open(filename,'w')
            #-  Array size reminder
            #- vbcw(2, JMax+3, KMax+4, Nrec);
            #- vbce(2, JMax+3, KMax+4, Nrec);
            #- vbcn(IMax+3, 3, KMax+4, Nrec);
            #- vbcs(IMax+3, 3, KMax+4, Nrec);
            for k in range(KMax+4):
                for j in range(JMax+3):
                    for i in range(2):
                        fileID.write('%26.16f\n' % bcw[i,j,k,rec])
            for k in range(KMax+4):
                for j in range(JMax+3):
                    for i in range(2):
                        fileID.write('%26.16f\n' % bce[i,j,k,rec])
            for k in range(KMax+4):
                for j in range(2):
                    for i in range(IMax+3):
                        fileID.write('%26.16f\n' % bcs[i,j,k,rec])
            for k in range(KMax+4):
                for j in range(2):
                    for i in range(IMax+3):
                        fileID.write('%26.16f\n' % bcn[i,j,k,rec])
            fileID.close()
            print('*****************************************')
            print('***BC for w are ready:  '+filename)
            print('*****************************************')
    elif(variable == 'p'):
        newdir = bc_dir+'/'+newdirvar      #-  ./BCS_ROMS/u_roms_bdy_SD
        if(not os.path.exists(newdir)):
            os.makedirs(newdir)
        #- Loop through each times
        for rec in (range(Nrec)):
            if (len(str(rec)) == 1): IDno='00000'+str(rec)
            if (len(str(rec)) == 2): IDno='0000'+str(rec)
            if (len(str(rec)) == 3): IDno='000'+str(rec)
            filename = newdir+'/p_bry_'+IDno+'.dat'
            fileID = open(filename,'w')
            #-  Array size reminder
            #- vbcs(IMax+3, JMax+3, Nrec);
            for j in range(JMax+3):
                for i in range(IMax+3):
                    fileID.write('%26.16f\n' % bcn[i,j,rec])
            fileID.close()
            print('*****************************************')
            print('***BC for p are ready:  '+filename)
            print('*****************************************')



#----INTERP_ROMS2GCCOM -----------------------------------------------
#- From GCCOM matlab code interp_roms2gccom.m
'''
   This code interpolates ROMS data onto GCCOM grid 
   INPUTS:
          lonR := Vector of longitude values in ROI
          latR := Vector of latitude values in ROI
          Rdepth := Vector of depth values in ROMS data
          u := 3-D array of u-values
          lonBC,latBC,levBC := Vectors representing values at opening face
   OUTPUTS:
          u_interp,u_real
   u_interp,u_real = interp_roms2gccom(lonR,latR,depthR,u,lonBC,latBC,levBC)
'''
#---------------------------------------------------------------------
def interp_roms2gccom(lonR,latR,depthR,u,lonBC,latBC,levBC):
    #- Initialize empty return variables
    u_interp = np.zeros(lonBC.shape)
    u_real = np.zeros(lonBC.shape)
    for l in range(len(lonBC)):
        lon_reg = lonBC[l]
        lat_reg = latBC[l]
        lev_reg = levBC[l]
        #- indices of closest point  
        ind_lon, = np.where(lonR<=lon_reg)
        ind_lon = ind_lon[-1]
        ind_lat, = np.where(latR<=lat_reg)
        ind_lat = ind_lat[-1]
        ind_lev, = np.where(depthR<=lev_reg)
        ind_lev = ind_lev[0]

        x_01 =(lon_reg-lonR[ind_lon])/(lonR[ind_lon+1]-lonR[ind_lon])
        y_01 =(lat_reg-latR[ind_lat])/(latR[ind_lat+1]-latR[ind_lat])
        z_01 =(lev_reg-depthR[ind_lev])/(depthR[ind_lev-1]-depthR[ind_lev]); 

        #if (l == 0):
        #    print(str(ind_lon)+','+str(ind_lat)+','+str(ind_lev))
        #    print(str(x_01)+','+str(y_01)+','+str(z_01))

        #-  to test or troubleshoot this function, uncomment the following if..end 
        #if ((x_01>1)|(x_01<0)|(y_01>1)|(y_01<0)|(z_01>1)|(z_01<0)):
        #    print('Error in interp_roms2gccom: local x or y or z outside of range 0 to 1')
        u000 = u[ind_lon  ,ind_lat  ,ind_lev  ];
        u100 = u[ind_lon+1,ind_lat  ,ind_lev  ];
        u010 = u[ind_lon  ,ind_lat+1,ind_lev  ];
        u001 = u[ind_lon  ,ind_lat  ,ind_lev-1];
        u101 = u[ind_lon+1,ind_lat  ,ind_lev-1];
        u011 = u[ind_lon  ,ind_lat+1,ind_lev-1];
        u110 = u[ind_lon+1,ind_lat+1,ind_lev  ];
        u111 = u[ind_lon+1,ind_lat+1,ind_lev-1];

        u_real[l]=u[ind_lon,ind_lat,ind_lev];
        
        #------------------------------------------------------------------
        #- Remnant code commented out 
        #u_interp(l) = u(ind_lon,ind_lat,ind_lev) + \
        #               (lon_reg-lonR(ind_lon))/(lonR(ind_lon+1) - \
        #               lonR(ind_lon))*(u(ind_lon+1,ind_lat,ind_lev) - \
        #               u(ind_lon,ind_lat,ind_lev)) + \
        #               (lat_reg-latR(ind_lat))/(latR(ind_lat+1) - \
        #               latR(ind_lat))*(u(ind_lon,ind_lat+1,ind_lev) - \
        #               u(ind_lon,ind_lat,ind_lev)) + \
        #               (lev_reg-depthR(ind_lev))/(depthR(ind_lev-1) - \
        #               depthR(ind_lev))*(u(ind_lon,ind_lat,ind_lev-1) - \ 
        #               u(ind_lon,ind_lat,ind_lev))
        #------------------------------------------------------------------

        u_interp[l] = u000*(1-x_01)*(1-y_01)*(1-z_01) + \
                      u100*(  x_01)*(1-y_01)*(1-z_01) + \
                      u010*(1-x_01)*(  y_01)*(1-z_01) + \
                      u001*(1-x_01)*(1-y_01)*(  z_01) + \
                      u101*(  x_01)*(1-y_01)*(  z_01) + \
                      u011*(1-x_01)*(  y_01)*(  z_01) + \
                      u110*(  x_01)*(  y_01)*(1-z_01) + \
                      u111*(  x_01)*(  y_01)*(  z_01)

    #- Remnant plotting code
    # figure
    # scatter3(lonBC,latBC,levBC,5,u_interp,'filled');
    # hold on
    # scatter3(temp_lon,temp_lat,temp_lev,10,u_real,'filled');
    # colormap(jet);
    # colorbar;
    # legend('interpolated','ROMS')

    return u_interp,u_real




#----GCCOM_GHOST_BC.M -----------------------------------------------
#- From GCCOM matlab code gccom_ghost_bc.m
'''
   This code adds ghost points to GCCOM mesh
   INPUTS:
          variable := string denoting which face of cube
          u_in := 2-D array of interpolated values 
          IMax,JMax,KMax := Dimensions of 3D Mesh
   OUTPUTS:
          u_inall
   u_inall = gccom_ghost_BC('ubcw',u_in, IMax,JMax,KMax)
'''
#---------------------------------------------------------------------
def gccom_ghost_BC(variable,varin, IMax,JMax,KMax):
    if ((variable == 'ubcw') | (variable == 'ubce')):
        #- Initialize empty array
        #BC_var = zeros(3,JMax+3,KMax+3);
        BC_var = np.zeros([3,JMax+3,KMax+3])

        #- Fill in Ghost points closest to GCCOM with u-values of grid
        #BC_var[3,3:JMax+1,1:KMax-1] = varin(1:JMax-1,1:KMax-1); 
        BC_var[2,2:JMax,0:KMax-2] = varin[0:JMax-2,0:KMax-2]

        #- Fill in Ghost points closest + 1 to GCCOM with u-values of 
        #- previous ghost points
        #BC_var(2,3:JMax+1,3:KMax+1) = BC_var(3,3:JMax+1,3:KMax+1); 
        BC_var[1,2:JMax,2:KMax] = BC_var[2,2:JMax,2:KMax] 

        #- Fill in Ghost points closest + 2 to GCCOM with u-values of prev
        #BC_var(1,3:JMax+1,3:KMax+1) = BC_var(3,3:JMax+1,3:KMax+1);
        BC_var[0,2:JMax,2:KMax] = BC_var[2,2:JMax,2:KMax]

        #- Fill two ghost cells above and below using vertical BCs: 

        #- no-slip  on Y wall at upper bounds
        #BC_var(:,JMax+2,:)=-BC_var(:,JMax+1,:); 
        BC_var[:,JMax+1,:] = -BC_var[:,JMax,:]
        #BC_var(:,JMax+3,:)=-BC_var(:,JMax,:);
        BC_var[:,JMax+2,:] = -BC_var[:,JMax-1,:]

        #- no-slip on Y wall at lower bounds
        #BC_var(:,2,:)=-BC_var(:,3,:); 
        BC_var[:,1,:] = -BC_var[:,2,:]
        #BC_var(:,1,:)=-BC_var(:,4,:);
        BC_var[:,0,:] = -BC_var[:,3,:]

        #- no-slip at bottom Z wall
        #BC_var(:,:,2)=-BC_var(:,:,3); 
        BC_var[:,:,1] = -BC_var[:,:,2]
        #BC_var(:,:,1)=-BC_var(:,:,4);
        BC_var[:,:,0] = -BC_var[:,:,3]

        #-free-slip at top Z wall
        #BC_var(:,:,KMax+2)=BC_var(:,:,KMax+1); 
        BC_var[:,:,KMax+1] = BC_var[:,:,KMax]
        #BC_var(:,:,KMax+3)=BC_var(:,:,KMax  );
        BC_var[:,:,KMax+2] = BC_var[:,:,KMax-1]
        #---- END case UBCW -------------------------

    elif ((variable == 'vbcw') | (variable == 'vbce')):
        #- Fill in Ghost points closest to GCCOM with u-values of grid
        #BC_var = zeros(2, JMax+4, KMax+3);
        BC_var = np.zeros([2,JMax+4,KMax+3])

        #BC_var(2,3:JMax+2,3:KMax+1) = varin(1:JMax,1:KMax-1);
        BC_var[1,2:JMax+1,2:KMax] = varin[0:JMax-1,0:KMax-2]

        #BC_var(1,3:JMax+2,3:KMax+1) = BC_var(2,3:JMax+2,3:KMax+1);
        BC_var[0,2:JMax+1,2:KMax] = BC_var[1,2:JMax+1,2:KMax]
  
        #BC_var(:,2,:)=BC_var(:,3,:);
        BC_var[:,1,:] = BC_var[:,2,:]
        #BC_var(:,1,:)=BC_var(:,4,:);
        BC_var[:,0,:] = BC_var[:,3,:]
        #BC_var(:,JMax+3,:)=BC_var(:,JMax+2,:);
        BC_var[:,JMax+2,:] = BC_var[:,JMax+1,:]
        #BC_var(:,JMax+4,:)=BC_var(:,JMax+1,:  );  
        BC_var[:,JMax+3,:] = BC_var[:,JMax,:]  

        #BC_var(:,:,2)=-BC_var(:,:,3);
        BC_var[:,:,1] = -BC_var[:,:,2]
        #BC_var(:,:,1)=-BC_var(:,:,4);
        BC_var[:,:,0] = -BC_var[:,:,3]
        #BC_var(:,:,KMax+2)=BC_var(:,:,KMax+1);
        BC_var[:,:,KMax+1] = BC_var[:,:,KMax]
        #BC_var(:,:,KMax+3)=BC_var(:,:,KMax  );
        BC_var[:,:,KMax+2] = BC_var[:,:,KMax-1]
        #---- END case VBCW -------------------------

    return BC_var


