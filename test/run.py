import sys, os

fpath = os.path.dirname(__file__)
geomIOpath = os.path.join(fpath, '../src/python')
sys.path.append(geomIOpath)
import geomio
import numpy as np
import time as t
from stl import mesh
#name of the inputfile
inFile = "input/dome.svg"

#name of the stl file which will be generated
name = ["dome.stl"]#, "salt.stl"]

#wether it is a closed volume or an open surface
Volume = False

#how many layers are interpolated in total
numInterLayers = 3

#number of points that are computed per bezier segment
nPrec = 4

# import pickle
# os.chdir('/home/lucas/Desktop/pymarkers')
# with open("grid.pk1", 'rb') as inp:
#     grid = pickle.load(inp)
# # x,y,z = grid[0], grid[1], grid[2]

# os.chdir('/home/lucas/Desktop/geomIO/test')
# Phase = geomio.rayTracing("output/dome.stl", grid)

#mode bin or asc
#calling the main function
geomio.geomioFront(inFile,numInterLayers, nPrec, name, Volume, xml=False)
# 
# t1 = t.time()
# Phase = geomio.rayTracing(inFile, numInterLayers, nPrec, grid)

# t2 = t.time()
# print(t2-t1)
#plot the pointcloud. requiers open3d
# os.chdir('/home/lucas/Desktop/pymarkers')
# with open('Phase.npy', 'wb') as f:
#     np.save(f,Phase)

#geomio.plotCloud3D(inFile,numInterLayers,nPrec)
#l,nl = geomio.getLayers(inFile)

#get the boundaries of the mesh/volume
#geomio.getBounds(inFile, numInterLayers, nPrec)