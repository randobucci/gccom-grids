#!/bin/bash
#------------------------------------------------------------------------- 
#- gridgen_gridbathy.sh: Will interpolate bathymetry values onto
#-                       grid center points (*CE).  
#- INPUT:  ../input/sd_bathy.dat (xyz bathymetry file)
#-         ../output/sdbay-grid.xy
#- 
#- OUTPUT: ../output/output_bathy_CE (depth values at center points)
#-         ../output/output_grid_CE (grid center points)
#-
#- REQUIRES: gridbathy (Run "gridbathy -h" for more information.)
#------------------------------------------------------------------------- 

dataDir='../input'
gridDir='../output'
binDir='../software/gridutils'
bathyFile=$dataDir'/mont_bathy.xyz'
gridFile=$gridDir'/mont_grid.xy'     #- Corner nodes of grid
gridFileCE=$gridDir'/mont_grid_CE.xy'   #- Center nodes of grid
nodeMode='CE'      #- Output eiher CO or CE coordinates
outFile=$gridDir'/output_bathy_CE'
outFile2=$gridDir'/output_grid_CE'

#- Interpolate bathy onto grid points using Sibsonian nearest neighbor
$binDir/gridbathy -b $bathyFile -g $gridFile -a 2 -i CO > $outFile
/bin/cp -r $gridFileCE $outFile2
