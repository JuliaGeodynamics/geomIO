import sys, os

fpath = os.path.dirname(__file__)
geomIOpath = os.path.join(fpath, '../src/python')
sys.path.append(geomIOpath)
import geomio
import numpy as np

from stl import mesh
#name of the inputfile
inFile = "input/over.svg"

#name of the stl file which will be generated
name = ["2Layer.stl"]

#wether it is a closed volume or an open surface
Volume = False

#how many layers are interpolated in total
numInterLayers = 5

#number of points that are computed per bezier segment
nPrec = 12
#print(os.getcwd())
#grid = np.array([])
#grid = np.load("grid.npy")
#triangles = mesh.Mesh.from_file("output/fold.stl")


data = geomio.readSVG(inFile)
c,l =geomio.getPoints2D(inFile, nPrec, True)
#coors = geomio.getPoints2D(inFile,nPrec)
#mode bin or asc
#calling the main function
#geomio.geomioFront(inFile,numInterLayers, nPrec, name, Volume, xml=True)
#plot the pointcloud. requiers open3d
#Phase = geomio.rayTracing(inFile, numInterLayers, nPrec, grid)
#np.save("Phase.npy",Phase)
#np.save("phase.npy", Phase)
#geomio.plotCloud3D(inFile,numInterLayers,nPrec)
#l,nl = geomio.getLayers(inFile)

#get the voundaries of the mesh/volume
#geomio.getBounds(inFile, numInterLayers, nPrec)