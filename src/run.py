import sys, os
import geomio

"""
To do:
    -fix closed volumes#done 
    -add variable stl#done
    -add inside meshlayer
    
 
"""
#name of the inputfile
inFile = "input/livedemo.svg"

#name of the stl file which will be generated
name = 'live.stl'

#wether it is a closed volume or an open surface
Volume = False

#how many layers are interpolated in total
numInterLayers = 8

#number of points that are computed per bezier segment
nPrec = 4

#calling the main function
geomio.geomioFront(inFile,numInterLayers, nPrec, name, Volume)
#plot the pointcloud. requiers open3d

#geomio.plotCloud3D(inFile,numInterLayers,nPrec)
#l,nl = geomio.getLayers(inFile)

#get the voundaries of the mesh/volume
#geomio.getBounds(inFile, numInterLayers, nPrec)



