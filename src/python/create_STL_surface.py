
from curve_interpolations import *
from svgpathtools import svg2paths, real, imag, Line, svg2paths2, Document
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as Pt
import matplotlib.patches as patches
import math
import sys,os
import scipy as sc
from scipy import interpolate
import struct
import stl
from stl import mesh
#from numba import jit

class stlMesh(NamedTuple):
    paths      : list
    nInterPaths: int 
    nBezCtrlPts: int
    isVol      : bool
    mode       : str


#@jit(nopython=True) 
def indexingOpen(numLayers, lc, nInter):
    nQuads = nInter+numLayers -1 # number of quads to be split in two triangles
    tri1 = np.array([])#first set of triangles
    tri2 = np.array([])#second set of second set of triangles
    for p in range(nQuads):
        numP = int(len(lc)/(nQuads+1))#number of points per layer
        nodeIdx = numP *p #index reffering to current layer
        #shape = (int(numP*2)-2,3) # shape of index array for triangles of current layer + next
        shapeS = (numP-1,3)# shape of indexing array of triangles 
        ind = np.zeros(shapeS, dtype = int)#indices to first ste of triangles
        ind2 = np.zeros(shapeS,dtype=int)# indices to second set of triangles
        #arrshift = numP-1
        #counter = 0 # number of triangles
        if p == 0:
            for i in range(numP-1):
                if i == 0:
                    continue
                ind[i-1,0] = i
                ind[i-1,1] = i+1
                ind[i-1,2] = numP+i
                
                ind2[i-1,0] = numP +i
                ind2[i-1,1] = numP + i+1
                ind2[i-1,2] = i+1
          #      counter +=1
            tri1 = ind[0:-1,:]
            tri2 = ind2[0:-1,:]
        else:
    
            for i in range(numP-1):
                if i == 0:
                    continue
                ind[i-1,0] = nodeIdx + i
                ind[i-1,1] = nodeIdx + i+1
                ind[i-1,2] = nodeIdx + numP+i
                
                ind2[i-1,0] = nodeIdx + numP +i
                ind2[i-1,1] = nodeIdx + numP + i+1
                ind2[i-1,2] = nodeIdx+ i+1
           #     counter +=1
                
            tri2 = np.concatenate((ind2[0:-1,:], tri2))    
            tri1 = np.concatenate((ind[0:-1,:], tri1))
    
    triangles = np.concatenate((tri1,tri2))
    normals = triNormals(triangles, lc)

    mesh1 = np.array([])
    #print(counter)
    
    
            
    for p in range(len(triangles)):
        #first ind = first tri
        mesh1 = np.concatenate((mesh1, normals[p]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,0]]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,1]]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,2]]))
        
    mesh1 = np.reshape(mesh1,(len(triangles)*4,3))
    
    return mesh1, triangles

      

#@jit(nopython=True) 
def triNormals(tri,lc):
    """
    Function to compute the normal of a triangle
    currently unused
    """
    normals = np.zeros((len(tri),3))
    for i in range(len(tri)):
        #print(i)
        p1 = lc[tri[i,0],:]
        p2 = lc[tri[i,1],:]
        p3 = lc[tri[i,2],:]
        norm = np.cross(p2-p1, p3-p1)
        normals[i,:] = norm
    return normals
def triSurfOpen(lc, numLayers, nInter, nPrec):
    """
    surface triangulation for open surfaces
    """
    #Layers, numLayers = getLayers(inFile)
    
    #lc = getPointCoords(inFile,nInter,nPrec) # get coordinates
    
    m1, t = indexingOpen(numLayers, lc, nInter)
    
    return m1, t
    
   

from scipy.spatial import ConvexHull, convex_hull_plot_2d,Delaunay 
    
#@jit(nopython=True) 
def layerTri(lc, numP, LLI, last = True):
    

    
    if (numP % 2) == 0:
        shape = (int((numP-2)/2),3)
        index1 = np.zeros(shape, dtype = int)
        index2 = np.zeros(shape, dtype = int)
        for i in range(len(index1)):
            index1[i,0] = i
            index1[i,1] = i+1
            index1[i,2] = numP-1-i
            
            index2[i,0] = i
            index2[i,1] = numP-1-i
            index2[i,2] = numP -2 -i
    
        triangles = np.concatenate((index1,index2))

    else:
        shape = (int((numP-3)/2),3)
        index1 = np.zeros(shape, dtype = int)
        index2 = np.zeros(shape, dtype = int)
        for i in range(len(index1)):
            index1[i,0] = i+1
            index1[i,1] = i+2
            index1[i,2] = numP-1-i
            
            index2[i,0] = i+1
            index2[i,1] = numP-i -1
            index2[i,2] = numP -2 -i
    
        first = np.array([[0,1,numP-1]])
        triangles = np.concatenate((index1,index2))
        triangles = np.concatenate((triangles, first))
        
    normals = triNormals(triangles, lc)

    mesh = np.array([])
    for p in range(len(triangles)):
        #first ind = first tri
        mesh = np.concatenate((mesh, normals[p]))
        mesh = np.concatenate((mesh, lc[triangles[p,0]]))
        mesh = np.concatenate((mesh, lc[triangles[p,1]]))
        mesh = np.concatenate((mesh, lc[triangles[p,2]]))
        
    mesh = np.reshape(mesh,(len(triangles)*4,3))
   
    if last and (numP % 2) == 0:
        trueIdx = np.concatenate((index1 + LLI, index2 + LLI))
        triangles = trueIdx
    elif last:
        trueIdx = np.concatenate((index1 + LLI, index2 + LLI))
        triangles = np.concatenate((trueIdx, first + LLI))

        
            
    return mesh, triangles
#@jit(nopython=True) 
def triangulateCover(lc):
    numP = len(lc)
    shape = len(lc)-2
    ind = np.zeros((shape,3))
    triangles = np.zeros((shape,3), dtype = int)
    for i in range(shape):
        triangles[i,0] = i 
        triangles[i,1] = i+1
        triangles[i,2] = numP-2-i
    normals = triNormals(triangles, lc)
    mesh = np.array([])
    for p in range(len(triangles)):
        #first ind = first tri
        mesh = np.concatenate((mesh, normals[p]))
        mesh = np.concatenate((mesh, lc[triangles[p,0]]))
        mesh = np.concatenate((mesh, lc[triangles[p,1]]))
        mesh = np.concatenate((mesh, lc[triangles[p,2]]))
        
    mesh = np.reshape(mesh,(len(triangles)*4,3))
    
    return  mesh, triangles

def triSurfClose(lc, numLayers, nInter, nPrec):
    """
    surface triangulation for closed volumes
    """
    #Layers, numLayers = getLayers(inFile)
    #lc = getPointCoords(inFile,nInter,nPrec)

    nQuads = nInter+numLayers -1
    tri1 = np.array([])
    tri2 = np.array([])
    
    #========Section where the "covers" are done : faulty
    #number of points per layer
    numP = int(len(lc)/(nQuads+1))
    #get x and y coordinates from first and last layer
    #cov1 = lc[0:numP,0:2]
    cov1Z = lc[0:numP,:]
    #last layer index, start of last layer
    LLI = len(lc)-numP
    #cov2 = lc[LLI:-1,0:2]
    cov2Z =lc[LLI:-1,:]
    dummy = np.array([lc[len(lc)-1]])
    cov2Z = np.concatenate((cov2Z,dummy))
    # startCov = Delaunay(cov1,furthest_site=True)
    # endCov = Delaunay(cov2, furthest_site=True)
    
    #cCF, cTF = layerTri(cov1Z,numP,LLI, False)
    #cCL, cTL = layerTri(cov2Z, numP, LLI, True)
    cCF, cTF = triangulateCover(cov1Z)
    cCL, cTL = triangulateCover(cov2Z)
    cTL += LLI
    #coordinatses of both covers
    cover = np.concatenate((cCF,cCL))
    #index of the triangles in lc
    covTri = np.concatenate((cTF,cTL))
    nQuads = nInter+numLayers -1 # number of quads to be split in two triangles
    tri1 = np.array([])#first set of triangles
    tri2 = np.array([])#second set of second set of triangles
    for p in range(nQuads):
        numP = int(len(lc)/(nQuads+1))#number of points per layer
        nodeIdx = numP *p #index reffering to current layer
        #shape = (int(numP*2)-2,3) # shape of index array for triangles of current layer + next
        shapeS = (numP-1,3)# shape of indexing array of triangles 
        ind = np.zeros(shapeS, dtype = int)#indices to first ste of triangles
        ind2 = np.zeros(shapeS,dtype=int)# indices to second set of triangles
        #arrshift = numP-1
        #counter = 0 # number of triangles
        if p == 0:
            for i in range(numP-1):
                if i == 0:
                    continue
                ind[i-1,0] = i
                ind[i-1,1] = i+1
                ind[i-1,2] = numP+i
                
                ind2[i-1,0] = numP +i
                ind2[i-1,1] = numP + i+1
                ind2[i-1,2] = i+1
          #      counter +=1
            tri1 = ind[0:-1,:]
            tri2 = ind2[0:-1,:]
        else:
    
            for i in range(numP-1):
                if i == 0:
                    continue
                ind[i-1,0] = nodeIdx + i
                ind[i-1,1] = nodeIdx + i+1
                ind[i-1,2] = nodeIdx + numP+i
                
                ind2[i-1,0] = nodeIdx + numP +i
                ind2[i-1,1] = nodeIdx + numP + i+1
                ind2[i-1,2] = nodeIdx+ i+1
           #     counter +=1
                
            tri2 = np.concatenate((ind2[0:-1,:], tri2))    
            tri1 = np.concatenate((ind[0:-1,:], tri1))
    
    triangles = np.concatenate((tri1,tri2))
    normals = triNormals(triangles, lc)

    mesh1 = np.array([])
    #print(counter)
    
    
            
    for p in range(len(triangles)):
        #first ind = first tri
        mesh1 = np.concatenate((mesh1, normals[p]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,0]]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,1]]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,2]]))
        
    mesh1 = np.reshape(mesh1,(len(triangles)*4,3))
    # CoorR = np.zeros_like(mesh1)

   
    #mesh1 = np.reshape(mesh1,(len(triangles)*4,3))
    mesh1 = np.concatenate((mesh1,cover))
    triangles = np.concatenate((triangles,covTri))

    return mesh1, triangles

from pointcloud_delaunay import *
import vtk
def plotInteractive(svg:svgFileData, stlctx:stlMesh):
    lc = getPointCoords(svg, stlctx.paths, stlctx.nInterPaths,stlctx.nBezCtrlPts)


    if stlctx.isVol:
        triangles, face = triSurfClose(lc,len(stlctx.paths), stlctx.nInterPaths,stlctx.nBezCtrlPts)
    else:           
        triangles, face = triSurfOpen(lc, len(stlctx.paths), stlctx.nInterPaths,stlctx.nBezCtrlPts)
    import plotly.figure_factory as ff
    fig = ff.create_trisurf(x=lc[:,0], y=lc[:,1], z=lc[:,2],
                         colormap="Portland",
                         simplices=face,
                         title="Mesh")
    fig.show('browser')
    return
    

def writeSTL(fname, svg:svgFileData,stlctx:stlMesh, mode = "ASCII"):
    """
    placeholder function for writing stl files
    uses the lib np stl
    """
    lc = getPointCoords(svg, stlctx.paths, stlctx.nInterPaths,stlctx.nBezCtrlPts)


    if stlctx.isVol:
        triangles, face = triSurfClose(lc,len(stlctx.paths), stlctx.nInterPaths,stlctx.nBezCtrlPts)
    else:           
        triangles, face = triSurfOpen(lc, len(stlctx.paths), stlctx.nInterPaths,stlctx.nBezCtrlPts)
  

    cube = mesh.Mesh(np.zeros(face.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(face):
        for j in range(3):
            cube.vectors[i][j] = lc[f[j],:]
    
    if mode == "ASCII":       
        cube.save(str(fname),mode=stl.Mode.ASCII )
    elif mode == "BIN":
        cube.save(str(fname),mode=stl.Mode.BINARY )
    return

