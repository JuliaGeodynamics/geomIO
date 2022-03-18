import sys, os
import geomio
import numpy as np
import time as t
from meshGen import *


# Get geomIO context
ctx = geomio.context()

# Add SVG to geomIO context
ctx.add_svgInFile("input/dome.svg")

# Export STLs
ctx.export_svg2stl()


#name of the stl file which will be generated
#name = ["dome.stl"]#, "salt.stl"]
#for raytracing layers must be provided in a sequence from oldest to youngest structure
name = ["output/dome1.stl", "output/dome2.stl", "output/dome3.stl"]
#wether it is a closed volume or an open surface
Volume = False

#how many layers are interpolated in total
numInterLayers = 3

#number of points that are computed per bezier segment
nPrec = 4

xBound =[0,250]
yBound =[100,200]
zBound =[-200,0]
nx,ny,nz = 30,40,80

x,y,z,X,Y,Z = createMesh(nx,ny,nz, xBound, yBound, zBound)
grid = list((X,Y,Z))
t1 = t.time()
#Phase = geomio.rayTracing(name, grid)
t2 = t.time()
print(t2-t1)
#writeVTK(x,y,z,Phase, "domeNeW")

#calling the main function
geomio.geomioFront(ctx,numInterLayers, nPrec, name, Volume, xml=True)
# 

# Phase = geomio.rayTracing(inFile, numInterLayers, nPrec, grid)


#plot the pointcloud. requiers open3d
#geomio.plotCloud3D(inFile,numInterLayers,nPrec)

#get the boundaries of the mesh/volume
#geomio.getBounds(inFile, numInterLayers, nPrec)