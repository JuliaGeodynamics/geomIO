# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""



# import np as np
# import sys, os

# from svglib.svglib import svg2rlg
# from reportlab.graphics import renderPDF, renderPM
# #import svgpathtools as pt
# import math 


# We need the following python packages for this to work:
#
# svgpathtools, numpy, matplotlib, math, sys, os, scipy
# numpy-stl

from svgpathtools import svg2paths, real, imag, Line, svg2paths2, Document
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as Pt
import matplotlib.patches as patches
import math
import sys,os
#import ipdb
import scipy as sc
from scipy import interpolate

#import stlWrite
from pointcloud_delaunay import *

from read_svg import *
from curve_interpolations import *
from create_STL_surface import *

#os.chdir("..")



def scaling():
    return
    
    
    
def plot_line(inp, tol):    # Obsolete?
    #deprecated
    
    line = np.array(inp)#*tol
     
    #plt.plot(line[:, 0], line[:, 1], color='b')
    
    plt.scatter(line[:, 0], line[:, 1], s=5)




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
def rayTracing(inFile, numInter, nPrec, grid):

    
    Phase = OpenVolumeTest(inFile, numInter, nPrec, grid)
    
    return Phase






#es ist halt wirklich die interpolation die von anfang bis end inter aber in der mitte ignoriert




def getBounds(inFile, numInter, nPrec):
    """
    printing the min and ma boundaries of the mesh/volume
    to adapt grid of thermomechanical coed(eg LaMEM)
    """
    lc = getCarthesian(inFile, numInter,nPrec)
    xBound = np.array([np.amin(lc[:,0]),np.amax(lc[:,0])])
    yBound = np.array([np.amin(lc[:,1]),np.amax(lc[:,1])])
    zBound = np.array([np.amin(lc[:,2]),np.amax(lc[:,2])])
    print("---------------------------")
    print("x boundaries are:"+ str(xBound))
    print("y boundaries are:"+ str(yBound))
    print("z boundaries are:"+ str(zBound))
    print("---------------------------")
    return


# dlnSC = Delaunay(lc)
# hull = ConvexHull(lc)




#-----------Poisson rec------------

#poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8, width=0, scale=1.1, linear_fit=False)[0]
#inFile = "input/over.svg"
def plotCloud3D(inFile,nInterLayers, nPrec):
    """
    uses open3D to display the pointcloud
    """
    import open3d as o3d
    Layers, numLayers = getLayers(inFile)
    lc = getCarthesian(inFile,nInterLayers,nPrec)
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(lc[:,:3])
    #pcd.colors = o3d.utility.Vector3dVector(lc[:,3:6]/255)
    #pcd.normals = o3d.utility.Vector3dVector(lc[:,6:9])
    o3d.visualization.draw_geometries([pcd])
    return

#plotCloud3D(inFile)

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

#plotCas(SEG)




import matplotlib.pyplot as plt


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

def geomioFront(inFile, numInterLayers, nPrec, name, volume = False, mode = "ASCII"):
    
    t1 = time.time()
    wSTL(inFile, numInterLayers, nPrec, name, volume, mode)
    t2 = time.time()
    print(t2-t1)
    return
    




        
        

    