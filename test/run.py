import sys, os

fpath = os.path.dirname(__file__)
geomIOpath = os.path.join(fpath, '../src/python')
sys.path.append(geomIOpath)
import geomio
import numpy as np
import time as t
from meshGen import *
#name of the inputfile
inFile = "input/dome.svg"

#name of the stl file which will be generated
#name = ["dome.stl"]#, "salt.stl"]
name = ["dome1.stl", "dome2.stl", "dome3.stl"]
#wether it is a closed volume or an open surface
Volume = False

#how many layers are interpolated in total
numInterLayers = 3

#number of points that are computed per bezier segment
nPrec = 4
xBound =[0,250]
yBound =[100,200]
zBound =[-200,0]
nx,ny,nz = 40,60,80

x,y,z,X,Y,Z = testMesh(nx,ny,nz, xBound, yBound, zBound)
grid = list((X,Y,Z))
# t1 = t.time()
#Phase = geomio.rayTracing(name, grid)
# t2 = t.time()
# print(t2-t1)
#writeVTK(x,y,z,Phase, "dome")

#calling the main function
#geomio.geomioFront(inFile,numInterLayers, nPrec, name, Volume, xml=True)
# 

# Phase = geomio.rayTracing(inFile, numInterLayers, nPrec, grid)


#plot the pointcloud. requiers open3d
#geomio.plotCloud3D(inFile,numInterLayers,nPrec)

#get the boundaries of the mesh/volume
#geomio.getBounds(inFile, numInterLayers, nPrec)