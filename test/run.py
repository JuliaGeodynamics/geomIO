import sys,os


fpath = os.path.dirname(__file__)
geomIOpath = os.path.join(fpath, '../src/python')
sys.path.append(geomIOpath)



import geomio
import numpy as np
import time as t



import warnings
warnings.filterwarnings('ignore')

#### 1) Create geomIO context
ctx = geomio.context()

#### 2) Add SVG to geomIO context
ctx.addsvg("input/dome.svg")

#### 3) Export STLs
nInterPaths = 12     # Refine the grid in normal direction by adding N interpolated path instances
nBezCtrlPts = 8    # Refine the grid along the Bezier paths by adding additional control points
isVol       = False # Closed (volume) or or open STL mesh, i.e., do we have closed or open Bezier path instances
ctx.svg2stl(nInterPaths,nBezCtrlPts,isVol, "BIN")
#ctx.plot(nInterPaths,nBezCtrlPts,isVol)

#### 4) Transfer STL layered structure to rectilinear grid
#    4a) Add rectilinear grid to geomIO context
xBnds    = [0,250]
yBnds    = [110,190]
zBnds    = [-200,0]
nx,ny,nz = 400,40,400
ctx.addgrd(nx,ny,nz,xBnds,yBnds,zBnds)

#    4b) Use ALL STL structures to assign the structure to the rectilinear grid

#    4b) [optional] Only use specific STL structures to assign the structure to the rectilinear grid
t1 = t.time()
ctx.stl2grd(["dome1","dome2","dome3"])
#ctx.stl2grd(["dome1","dome3"])
t2 = t.time()
print(t2-t1)


#### 5) Export phase grid as vtr (VTK rectilinear grid)
ctx.grd2vtr()







#name of the stl file which will be generated
#name = ["dome.stl"]#, "salt.stl"]
#for raytracing layers must be provided in a sequence from oldest to youngest structure
#name = ["output/dome1.stl", "output/dome2.stl", "output/dome3.stl"]
#wether it is a closed volume or an open surface
#Volume = False

#how many layers are interpolated in total
#numInterLayers = 3

#number of points that are computed per bezier segment
#nPrec = 4



#x,y,z,X,Y,Z = createMesh(nx,ny,nz, xBnds, yBnds, zBnds)
#grid = list((X,Y,Z))
#t1 = t.time()
#Phase = geomio.rayTracing(name, grid)
#t2 = t.time()
#print(t2-t1)
#writeVTK(x,y,z,Phase, "domeNeW")

#calling the main function
#geomio.geomioFront(ctx,numInterLayers, nPrec, name, Volume, xml=True)
# 
#Phase = geomio.rayTracing(ctx, numInterLayers, nPrec, grid)


#plot the pointcloud. requiers open3d
#geomio.plotCloud3D(inFile,numInterLayers,nPrec)

#get the boundaries of the mesh/volume
#geomio.getBounds(inFile, numInterLayers, nPrec)