'''
bathy-project.py: Project bathy csv file from Lat/Lon data to UTM zone 10 meters
Usage: python bathy-project.py <input_bathy>
'''
import os,sys
#from gis_utils import *
from ll2utm import *

#- Open up input file 
input_args = sys.argv
if (len(sys.argv) < 2):
    input_file = './input/monterey.xyz'
else:
    input_file = input_args[1]       #- Grid_sdbay.dat

#- Make sure input file exists
if (os.path.exists(input_file) == False):
    print("Cannot find ascii grid file: ./"+input_file)
    sys.exit()

#- Open up file and read coordinates
print("Reading Geo File:"+input_file)
fid = open(input_file,'r')
lon,lat,z = [],[],[]
for line in fid:
    line = line.strip()
    cols = line.split(',')
    lon.append(float(cols[0]))
    lat.append(float(cols[1]))
    z.append(float(cols[2]))
fid.close()


print("Projecting Data to UTM Coordinates")
zone, easting, northing = LL_to_UTM(23, lat, lon)

#- Write output
out_file = "./input/mont_bathy.xyz"
ofid = open(out_file,'w')
for i in range(len(easting)):
    ofid.write('%9.2f %10.2f %8.2f\n' % (easting[i],northing[i],z[i]))
ofid.close()
print("Created UTM File:"+out_file)


