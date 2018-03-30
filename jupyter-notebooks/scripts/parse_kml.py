#--------------------------------------------------------
#- File: parse_kml.py
#- Description: Extract coordinates from Google Earth KML
#-              and convert to UTM, saving results.
#- Author: Randy Bucciarelli
#- Date: June 15, 2017
#--------------------------------------------------------
import sys
from xml.dom.minidom import parseString

#- Define source files
input_args = sys.argv
if (len(sys.argv) < 2):
    src_dir = '../input/'
    kml_file = src_dir+'lj-outline.kml'
else:
    kml_file = input_args[1]

#- Open kml file and read into memory
file = open(kml_file)
theData = file.read()
file.close()

#-Load data string into DOM
theDom = parseString(theData)

lats,lons,verts = [],[],[]
#- Loop through dom and find coordinates
for d in theDom.getElementsByTagName('coordinates'):
    #- Get all the vertices in first polygon
    positions = d.firstChild.data.split()
    for p in positions:
        coords = p.split(',') 
        lats.append(float(coords[1]))
        lons.append(float(coords[0]))
        verts.append(0)

#- Write geographic coords to outfile
out_file = kml_file[0:len(kml_file)-4]+'_geo.xy'
thefile = open(out_file, 'w')
for i in range(len(lats)):
    thefile.write("%13.8f\t%13.8f\t%d\n" % (lons[i],lats[i],verts[i]))
thefile.close()
print("Created "+out_file)
