# -*- coding: utf-8 -*-

# We need the following python packages for this to work:
#
# svgpathtools, numpy, matplotlib, math, sys, os, scipy
# numpy-stl

from cmath import inf
from distutils.log import error
from svgpathtools import svg2paths, real, imag, Line, svg2paths2, Document
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as Pt
import matplotlib.patches as patches
import math
import sys,os,inspect

import scipy as sc
from scipy import interpolate


from pointcloud_delaunay import *

from read_svg import *
from curve_interpolations import *
from create_STL_surface import *

class rectLinGrid:
    xBnds       : tuple
    yBnds       : tuple
    zBnds       : tuple
    nx          : int
    ny          : int
    nz          : int
    X           : np.ndarray
    Y           : np.ndarray
    Z           : np.ndarray
    PHS         : np.ndarray
    grd         : list
    
    def __init__(self,nx : int ,ny: int, nz: int, xBnds:tuple, yBnds:tuple, zBnds:tuple):
        self.nx              = nx
        self.ny              = ny
        self.nz              = nx
        self.xBnds           = xBnds
        self.xBnds           = yBnds
        self.xBnds           = zBnds
        self.x               = np.linspace(xBnds[0],xBnds[1],nx)
        self.y               = np.linspace(yBnds[0],yBnds[1],ny)
        self.z               = np.linspace(zBnds[0],zBnds[1],nz)
        self.X,self.Y,self.Z = np.meshgrid(self.x,self.y,self.z, indexing ='ij') #matrix indexing{ij} array[i,j]
        self.PHS             = np.zeros_like(self.X)
        self.grd             = list((self.X,self.Y,self.Z))

class stlMesh(NamedTuple):
    nInterPaths: int 
    nBezCtrlPts: int
    isVol      : bool
class context:
    svg:   svgFileData
    rlgrd: rectLinGrid
    stl:   stlMesh

    def __init__(self):   
        self.fpath = os.path.dirname(__file__)
        self.svgInFile = ""
        self.cwd = os.path.split(inspect.getframeinfo(sys._getframe(1)).filename)[0]
        print('Current working directory: ', self.cwd)
        self.outDir=os.path.join(self.cwd,'output')
        if not os.path.isdir(self.outDir):
            os.makedirs(self.outDir)
            print(' -> No output directory. Create: ', self.outDir)
        else:
            print(' -> Existing ouput directory: ', self.outDir)

    def addsvg(self,fname):
        self.svgInFile = os.path.join(self.cwd,fname)
        print('Loading ... ', self.svgInFile)
        self.svg = readSVG(self.svgInFile)

    def svg2stl(self,nInterPaths: int, nBezCtrlPts: int,isVol:bool):
        self.stl = stlMesh(nInterPaths,nBezCtrlPts,isVol)
        obj = list(set(self.svg.CurveNames))
        for i in range(len(obj)):
            path    = splitPaths(self.svg, obj[i])
            fname   = self.outDir + '/' + obj[i] + '.stl'
            print(' -> Export: ',  fname)
            wSTL(self.svg, path,
             self.numInterLayers, self.nPrec, fname,self.volume)
    
    def addgrd(self,nx : int ,ny: int, nz: int, xBnds:tuple, yBnds:tuple, zBnds:tuple):
        print('Create rectilinear grid: ', str(nx)+'x'+str(ny)+'x'+str(nz)+','+str(xBnds)+'x'+str(yBnds)+'x'+str(zBnds))
        self.rlgrd = rectLinGrid(nx,ny,nz,xBnds,yBnds,zBnds)
        

    def stl2grd(self,inFile:list):
        assert hasattr(self, 'rlgrd'), "Initialize grid for geomIO context with creategrd(...)"
        assert hasattr(self, 'svg'),   "Load svg file into geomIO context with readsvg(...)"
        
        if inFile:
            # Only use specific STLs/SVG objects to set up the grid
            for i in range(len(inFile)):
                fname          = self.outDir+'/'+inFile[i]+'.stl'
                print(' -> Assign to grid: ',  fname)
                phs            = fastRayFile(fname, self.rlgrd.grd)
                self.rlgrd.PHS = self.rlgrd.PHS + phs
        else:
            # Use ALL svg objects to set up the grid
            svgobjs = np.unique(np.array(self.svg.CurveNames))
            for objName in np.nditer(svgobjs):
                fname          = self.outDir+'/'+str(objName)+'.stl'
                print(' -> Assign to grid: ',  fname)
                phs            = fastRayFile(fname, self.rlgrd.grd)
                self.rlgrd.PHS = self.rlgrd.PHS + phs



    
class opt_():
    pass



def getPoints2D(inFile, nPrec, xml:bool = False):
    numInter = 1
    data = readSVG(inFile)
    Coors = list()
    
    if xml :
        labels = list(data.CurveNames)
        #name = list(labels)
        for i in range(len(labels)):
            name = str(labels[i])+ ".stl"
            path = splitPaths(data, labels[i])
            lc = getPointCoords(data, path, numInter,nPrec)
            Coors.append(lc)
    
    


    
    # x = data.Curves
    # cPoints = list()
    # for i in range(len(x)):
    #     if data.Commented[i]:
    #         continue
    #     else:
            
    #         line = x[i]  
    #         points = controlPoints(line)
    #         cPoints.append(points)
    # t = np.linspace(0,1, nPrec)
    # Coordinates = np.array([])

    # for r in range(len(cPoints)):
    #     for p in range(len(cPoints[r])):
    #         seg = cPoints[r][p]
    #         for s in range(nPrec):
    #             B = deCastel(seg,t[s])
                
                 
                
    #             Coordinates = np.append([B],Coordinates)

        
    
    # newshape = int(len(Coordinates)/2)
    # Coordinates = np.reshape(Coordinates,(newshape,2))
    return Coors, labels


def scaling():
    
    """
    tbi
    """
    return
    


def checkPath(attributes, Layers):  # Obsolete?
    """
    

    Parameters
    ----------
    Attributes and Layers from svgpathtools

    Returns
    -------
    number of paths per layer and also checks wether 
    the number of paths per layer is the same for all
    layers
    """
    
    #Check for dicts
    if not isinstance(attributes, list):
        sys.exit("Attributes must be list")
    if not isinstance(Layers, dict):
        sys.exit("Layers must be dictonary")
    numPaths = 0
    Counter = list(Layers.values())
    for i in Layers.values():
        if i in Counter[0]:
            numPaths +=1
    checkPath = 0
    for p in Layers.values():
        checkPath = 0
        for r in range(len(Counter)):
            if p in Counter[r]:
                checkPath +=1
        if checkPath != numPaths:
            sys.exit("Invalid number of Paths in " + str(p))
        #if i in attributes.values():
        #    print("dasads")
    return numPaths

def line2coor(path, tol = 1e-6):    # Obsolete?
    """
    DEPRACATED

    Parameters
    ----------
    path: path(s) from svgpathtools
    tol : tolerance, depraceteed.

    Returns
    -------
    coord : coordinates of lines

    """
    coord = []
     # read coordinates
    for i in range(len(path)) :
        
        # if not isinstance(path[i], Line) :
            
        #     sys.exit("All paths should consist of lines segments only")
        
        line =  path[i]
        p    =  line[0]
        pe   =  line[1]
        x    =  int(real(p))
        y    =  -int(imag(p))

        coord.append([x, y])
    end = path.end
    x    =  int(real(end))
    y    =  -int(imag(end))
    coord.append([x, y])
    return coord
    



def sortLayers(inFile): # Obsolete
    """
    This function checks where are to little or to much segments 
    if there are not equal number of paths(or segments) per layer    

    Parameters
    ----------
    inFile : inputfile
    Returns
    -------
    None.

    """
    paths, attributes = svg2paths(inFile)
    Layers, numLayers = getLayers(inFile)
    numPaths = checkPath(attributes, Layers)
    lenghts = list() 
    for i in range(numPaths):
        lenghts.append(len(paths[i]._segments))
    c = 0
    for p in range(numLayers):

        for r in range(numPaths):

            if len(paths[c]._segments) != lenghts[r]:
                sys.exit("wrong number of segments in Layer " + str(p)+ " considering your lowermost Layer as Layer 0")
                #print(Layers[p])
            c +=1      
    return 

from raytest import *



def rayTracing(inFile, grid:list, numInter: int = 2, nPrec:int = 2):
    
    if isinstance(inFile, list) :
        Phase = np.zeros_like(grid[0])
        for i in range(len(inFile)):
            Ph = fastRayFile(inFile[i] , grid)
            Phase = Phase +Ph
    else:
        if ".stl" in inFile:
            Phase = fastRayFile(inFile, grid)
    
        elif ".svg" in inFile:
            Phase = OpenVolumeTest(inFile, numInter, nPrec, grid)
            #catch case multiple objects
        
        else:
            sys.exit("Input must be provided as .svg or .stl file ")
    return Phase



def getBounds(inFile, numInter, nPrec):
    """
    printing the min and ma boundaries of the mesh/volume
    to adapt grid of thermomechanical coed(eg LaMEM)
    """
    lc = getPointCoords(inFile, numInter,nPrec)
    xBound = np.array([np.amin(lc[:,0]),np.amax(lc[:,0])])
    yBound = np.array([np.amin(lc[:,1]),np.amax(lc[:,1])])
    zBound = np.array([np.amin(lc[:,2]),np.amax(lc[:,2])])
    print("---------------------------")
    print("x boundaries are:"+ str(xBound))
    print("y boundaries are:"+ str(yBound))
    print("z boundaries are:"+ str(zBound))
    print("---------------------------")
    return






#-----------Poisson rec------------

#poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8, width=0, scale=1.1, linear_fit=False)[0]
#inFile = "input/over.svg"
def plotCloud3D(inFile,nInterLayers, nPrec):
    """
    uses open3D to display the pointcloud
    """
    import open3d as o3d
    Layers, numLayers = getLayers(inFile)
    lc = getPointCoords(inFile,nInterLayers,nPrec)
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(lc[:,:3])
    #pcd.colors = o3d.utility.Vector3dVector(lc[:,3:6]/255)
    #pcd.normals = o3d.utility.Vector3dVector(lc[:,6:9])
    o3d.visualization.draw_geometries([pcd])
    return


def plotCas(curve, nPoints = 50):
    """
    Plots the Bezier segments with de Casteljaus algorithm
    rather useless
    """
    t = np.linspace(0,1,nPoints)

    plt.figure()
    for p in range(len(curve)):
        c = curve[p]
        for i in range(nPoints):
        
            B = deCastel(c,t[i])
            plt.plot(B[0],B[1],'.')
    plt.show()

    return


def plotBezier(cPoints):
    """
    pyplot equivalent for plotting
    """
    fig, ax = plt.subplots()
    
    for i in range(len(cPoints)):
        verts =[(cPoints[i].start),
                (cPoints[i].c1),
                (cPoints[i].c2),
                (cPoints[i].end)]
        codes = [
            Pt.MOVETO,
            Pt.CURVE4,
            Pt.CURVE4,
            Pt.CURVE4,
        ]
        
        path = Pt(verts, codes)
        
        
        patch = patches.PathPatch(path, facecolor='none', lw=2)
        ax.add_patch(patch)
        
        xs, ys = zip(*verts)
        ax.plot(xs, ys, 'x--', lw=2, color='black', ms=10)
    plt.show()
    return


import time

def geomioFront(ctx:context, numInterLayers: int, nPrec: int, name:list, volume:bool = False, mode:str = "ASCII", xml: bool = False):
    
    if xml :
        labels = list(set(ctx.svg.CurveNames))
        #name = list(labels)
        for i in range(len(labels)):
            #name = str(labels[i]) + ".stl"
            path = splitPaths(ctx.svg, labels[i])
            wSTL(ctx.svg, path, numInterLayers, nPrec, ctx.outDir + '/' +  labels[i] + '.stl' , volume, mode)
            
    else:

        ## what is the purpose of this part? 

        path = ctx.svg.Curves 
        name = name[0]
        #t1 = time.time()
        wSTL(ctx.svg, path, numInterLayers, nPrec, name, volume, mode)
        #t2 = time.time()
        #print(t2-t1)
    return
    

def geomioMulitLayer(inFile, numInterLayers, nPrec, name, volume = False, mode = "ASCII"):
    
    
    
    return




        
        

    