import sys, os

fpath = os.path.dirname(__file__)
geomIOpath = os.path.join(fpath, '../src/python')
sys.path.append(geomIOpath)
import geomio
import numpy as np
import time as t
from stl import mesh
#name of the inputfile
inFile = "input/livedemo.svg"

#name of the stl file which will be generated
name = ["fold.stl"]#, "salt.stl"]

#wether it is a closed volume or an open surface
Volume = False

#how many layers are interpolated in total
numInterLayers = 3

#number of points that are computed per bezier segment
nPrec = 4
#print(os.getcwd())
#grid = np.array([])
#grid = np.load("grid.npy")
#triangles = mesh.Mesh.from_file("output/fold.stl")
import pickle
#os.chdir('/home/lucas/Desktop/pymarkers')
with open("gridL.pk1", 'rb') as inp:
    grid = pickle.load(inp)
x,y,z = grid[0], grid[1], grid[2]

#X,Y,Z = x[:,0,0], y[0,:,0], z[0,0,:]
X, Y = x[:,:,0], y[:,:,0]


#data = geomio.readSVG(inFile)
#c,l =geomio.getPoints2D(inFile, nPrec, True)
#coors = geomio.getPoints2D(inFile,nPrec)
#mode bin or asc
#calling the main function
#geomio.geomioFront(inFile,numInterLayers, nPrec, name, Volume, xml=True)
#os.chdir('/home/lucas/Desktop/geomIO/test')
#t1 = t.time()
Phase = geomio.rayTracing(inFile, numInterLayers, nPrec, grid)

#t2 = t.time()
#print(t2-t1)
#plot the pointcloud. requiers open3d
#os.chdir('/home/lucas/Desktop/pymarkers')
np.save("Phase.npy",Phase)
#np.save("phase.npy", Phase)
#geomio.plotCloud3D(inFile,numInterLayers,nPrec)
#l,nl = geomio.getLayers(inFile)

#get the voundaries of the mesh/volume
#geomio.getBounds(inFile, numInterLayers, nPrec)