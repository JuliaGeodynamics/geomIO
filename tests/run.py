import sys, os


x = os.getcwd()
#sys.path.append('/../geomIO/src/python')
os.chdir("..")
os.chdir("src/python")


import geomio
os.chdir(x)
#from geomIO.src.python import geomio
#name of the inputfile
inFile = "input/over.svg"

#name of the stl file which will be generated
name = 'folds2.stl'

#wether it is a closed volume or an open surface
Volume = False

#how many layers are interpolated in total
numInterLayers = 4

#number of points that are computed per bezier segment
nPrec = 40



#calling the main function
geomio.geomioFront(inFile,numInterLayers, nPrec, name, Volume)
#plot the pointcloud. requiers open3d

#geomio.plotCloud3D(inFile,numInterLayers,nPrec)
#l,nl = geomio.getLayers(inFile)

#get the voundaries of the mesh/volume
geomio.getBounds(inFile, numInterLayers, nPrec)



