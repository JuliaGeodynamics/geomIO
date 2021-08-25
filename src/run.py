import sys, os
import geomio 

#name of the inputfile
inFile = "input/demo.svg"

#name of the stl file which will be generated
name = 'deom.stl'

#wether it is a closed volume or an open surface
Volume = False

#how many layers are interpolated in total
numInterLayers = 5

#number of points that are computed per bezier segment
nPrec = 20

#calling the main function
#geomio.geomioFront(inFile,numInterLayers, nPrec, name, Volume)
#plot the pointcloud. requiers open3d

geomio.plotCloud3D(inFile,numInterLayers,nPrec)
#l,nl = geomio.getLayers(inFile)

#get the voundaries of the mesh/volume
#geomio.getBounds(inFile, numInterLayers, nPrec)



