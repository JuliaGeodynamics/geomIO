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

from svgpathtools import svg2paths, real, imag
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as Pt
import matplotlib.patches as patches
import sys

from pointcloud_delaunay import *
import meshwrite as wr
from read_svg import *
from curve_interpolations import *
from create_STL_surface import *
from typing import NamedTuple
from stl import mesh
#os.chdir("..")




def ComputeSurface(svgFileData,name, numInterLayers, prec, volume=False):
    """
        This computes a mesh and a triangulated surface from the curve "name", in svgFileData
        
    """

    CurveNames = svgFileData.CurveNames
    Commented  = svgFileData.Commented

    triangles = face=normals= Mesh = LayerCoords = None
    id = [i for i in range(len(CurveNames)) if ((CurveNames[i]==name) & (Commented[i]==False)) ]
    if len(id)==0:
        print("Curve: "+name+" is not present")

    elif len(id)==1:
        print("Curve: "+name+" is only present on a single layer; I will interprete this, but cannot create a triangular surface")
        
        # Compute intermediate curves:
        cPoints = controlPoints(svgFileData.Curves[id[0]])

        # Interpolate curves 
        LayerCoords, Mesh = CreateMesh([cPoints], [svgFileData.zCoord[id[0]]], prec)    

    else:
        # We can create a surface from this. The curves are already ordered from bot->top
        
        # Create list with curves & z-values on different layers
        cPoints = list()
        zCoor   = list()
        for i in id:
            points = controlPoints(svgFileData.Curves[i])
            cPoints.append(points)
            zCoor.append(svgFileData.zCoord[i])

        # Compute intermediate curves:
        cPoints, Zvals = InterpolateCurves(cPoints, zCoor, numInterLayers)

        # Compute 2D Array that has the surface
        LayerCoords, Mesh = CreateMesh(cPoints, Zvals, prec)
        
        # Triangulate surface & normals
        if volume==True:
            triangles, face, normals = triSurfClose_Compute(LayerCoords, numInterLayers, len(id))
            print("closed volume")
        else:
            triangles, face, normals = triSurfOpen_Compute(LayerCoords, numInterLayers, len(id))
            print("open surface")

        print("Found curve: " + name )


    return triangles, face, normals, Mesh, LayerCoords


def InterpolateCurves(cPoints, zCoor, numInterLayers):
        """
            Computes intermediate curves
        """

        # Compute z-values of intermediate layers:
        idx, Zvals, zAdded = compEmpty(zCoor, numInterLayers)
        Zvals = np.flip(Zvals)
        idx = np.flip(idx)
        
        print("Zvals="+str(Zvals))
        print("zCoor="+str(zCoor))

        # Interpolate curves
        interLayers = list()
        bezier      = list()
        for p in range(len(cPoints[0])):
            bezier = list()
            for r in range(len(cPoints)):
                bezier.append(cPoints[r][p])

            # obtain the intermediate curves (interpolated from fixed curves)
            seg = interp(bezier, zCoor, numInterLayers)
            interLayers.append(seg)

        reshape = list()
        for r in range(numInterLayers):
            layer = list()
            for p in range(len(cPoints[0])):    
                x = interLayers[p][r]
                layer.append(x)
            reshape.append(layer)

        print("idx="+str(idx))
        print("Zvals="+str(Zvals))
        print("zCoor="+str(zCoor))
        print("len(interLayers)="+str(len(interLayers)))
        print("len(cPoints)="+str(len(cPoints)))
        
        
        # Append the new curves @ the end of the already existing ones
        z_Total = []
        for i in range(len(zCoor)):
            z_Total.append(zCoor[i])  #

        for i in range(len(zAdded)):
            cPoints.append(reshape[i])
            z_Total.append(zAdded[i])   

        print("z_Total="+str(z_Total))

        # Now sort the curves according to z-value:
        li=[]
        for i in range(len(z_Total)):
            li.append([z_Total[i],i])
        li.sort()                           # sort from top->bottom to be consistent with rest of code
        sort_ind    =   [x[1] for x in li]

        cPoints_sort =   [cPoints[i]  for i in sort_ind]
        Zsort       =   [z_Total[i]  for i in sort_ind]

        print("Zsort="+str(Zsort))


        # This doesn'y "feel"
        #count = 0
        #for i in range(len(idx)):
        #    if idx[i] == 1:
        #        cPoints.insert(i, reshape[count])
        #        count +=1  

        print("len(cPoints)="+str(len(cPoints)))
        print("len(reshape)="+str(len(reshape)))

        return cPoints_sort, Zsort
        #return cPoints, Zvals


def CreateMesh(cPoints, Zvals, prec):
    """
        Create a mesh from various curves
    """
    t = np.linspace(0,1, prec)
    LayerCoors = np.array([])
    for r in range(len(Zvals)):
            
        for p in range(len(cPoints[0])):
            seg = cPoints[r][p]
            for s in range(prec):
                B = deCastel(seg,t[s])
                B = np.append(B,[Zvals[r]]) 
                    
                LayerCoors = np.append([B],LayerCoors)
        newshape = int(len(LayerCoors)/3)
        LayerCoors = np.reshape(LayerCoors,(newshape,3))
    
    Layer_X = np.reshape(LayerCoors[:,0], (prec*len(cPoints[0]), len(Zvals)))
    Layer_Y = np.reshape(LayerCoors[:,1], (prec*len(cPoints[0]), len(Zvals)))
    Layer_Z = np.reshape(LayerCoors[:,2], (prec*len(cPoints[0]), len(Zvals)))
    Mesh    = (Layer_X,Layer_Y,Layer_Z)

    return LayerCoors, Mesh

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




def geomioFront(inFile, numInterLayers, nPrec, name, volume = False, mode = "a"):
    
    wSTL(inFile, numInterLayers, nPrec, name, volume, mode)
    
    return
    



        
def wSTL(inFile, numInter, nPrec,name, volume = False, mode = "a"):
    """
     placeholder function for writing stl files
        uses the lib np stl
    """
    
    # = This can be replaced 

    import struct
    import stl
    from stl import mesh

    if 1==0:
        # Use the "old" method
        if volume:
            triangles, face = triSurfClose(inFile,numInter,nPrec)
        else:    
            triangles, face = triSurfOpen(inFile,numInter,nPrec)
        lc = getCarthesian(inFile, numInter,nPrec)
    else:
        # Employ the new approach, that can deal with multiple curves, scaling, Affinity Design etc.
        svgFileData = readSVG(inFile, True)
        CurveNames  = get_CurveNames(svgFileData) 
        triangles, face, normals, Mesh, lc = ComputeSurface(svgFileData,CurveNames[0], numInter, nPrec, volume)

    
    # Write *.stl file
    os.chdir("output")
    cube = mesh.Mesh(np.zeros(face.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(face):
        for j in range(3):
            cube.vectors[i][j] = lc[f[j],:]
    
    #Write the mesh to file "cube.stl"
    if mode == "a":
        cube.save(str(name),mode=stl.Mode.ASCII )
    elif mode == "b":
        cube.save(str(name),mode=stl.Mode.BINARY )
    os.chdir("..")
    return
    