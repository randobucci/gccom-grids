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
    input_file = '../input/mont_coast.ll'
else:
    input_file = input_args[1]       

#- Make sure input file exists
if (os.path.exists(input_file) == False):
    print("Cannot find ascii coastline: ./"+input_file)
    sys.exit()

#- Open up file and read coordinates
print("Reading Geo File:"+input_file)
fid = open(input_file,'r')
lon,lat = [],[]
for line in fid:
    line = line.strip()
    cols = line.split()
    lon.append(float(cols[0]))
    lat.append(float(cols[1]))
fid.close()


print("Projecting Data to UTM Coordinates")
zone, easting, northing = LL_to_UTM(23, lat, lon)

#- Write output
out_file = "../input/mont_coast.xy"
ofid = open(out_file,'w')
for i in range(len(easting)):
    ofid.write('%9.2f %10.2f\n' % (easting[i],northing[i]))
ofid.close()
print("Created UTM File:"+out_file)


