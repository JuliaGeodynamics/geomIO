#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 17:03:02 2021

@author: lucas
"""
from time import time
import numpy as np
import matplotlib.path as mpltPath
from rectLinGrid import *
import matplotlib.tri as mtri

from numba import jit
from numba.np.extensions import cross2d

@jit(nopython=True)
def insideTriangle2D(point:np.array, vertices:np.array):
    """
    determine if point is inside 2D triangle using barycentric weights
    https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution    
    https://math.stackexchange.com/questions/51326/determining-if-an-arbitrary-point-lies-inside-a-triangle-defined-by-three-points
    Parameters
    ----------
    point : Point to test for
    vertices : vertices of 2D triangle
        

    Returns
    -------
    Bool

    """
    #point = point.astype(float)
    #vertices = vertices.astype(float)
    a = vertices[0]
    b = vertices[1]
    c = vertices[2]
    
    xd = cross2d(a,b) + cross2d(b,c) + cross2d(c,a)
    #xd = np.cross(a,b) + np.cross(b,c) + np.cross(c,a)
    if xd ==0:
        xd+=0.00001

    xa = cross2d(b,c) + cross2d(point, b-c)
    xb = cross2d(c,a) + cross2d(point, c-a)
    xc = cross2d(a,b) + cross2d(point, a-b)    
    #xa = np.cross(b,c) + np.cross(point, b-c)
    #xb = np.cross(c,a) + np.cross(point, c-a)
    #xc = np.cross(a,b) + np.cross(point, a-b)
    
    wa = xa/xd
    wb = xb/xd
    wc = xc/xd
    
    
    if 0<wa<1 and 0<wb<1 and 0<wc<1:
        return True
    else:
        return False


@jit(nopython=True)    
def createPlane(p:np.array,q:np.array,r:np.array):
    """
    Creates a plane in Coordinate form
    ax + by +cz = d
    where a,b,c are the components of the normal array and d is d
    """    
    PQ = q-p
    PR = r-p
    OP = p
    normal = np.cross(PQ,PR)
    d = np.dot(normal, OP)
    return normal, d

@jit(nopython=True)
def upperTest(v1:np.array,v2:np.array,v3:np.array, zVal:np.array):
    """
    https://abiturma.de/mathe-lernen/geometrie/lagebeziehungen-und-schnitt/schnitt-gerade-ebene
    """
    Coor, d = createPlane(v1,v2,v3)
    r = ((d - zVal[0]*Coor[0]-zVal[1]*Coor[1])/Coor[2])-zVal[2]
    InterPoint = zVal+r*np.array([0,0,1])
    return InterPoint[2]






def rayTracingRectilinear(triangles, lc, grid):
    
    
    x,y,z = grid[0], grid[1], grid[2]
    X,Y = x[:,:,0], y[:,:,0]
    Z = z[0,0,:]
    Phase = np.zeros_like(x) # Phase array(comes notflattend)

    nx, ny, nz = Phase.shape[0],Phase.shape[1],Phase.shape[2]

    numberInter = np.zeros_like(X)
    #interTri = np.array([])
    InterPointZ = np.array([])
    #numberInter = np.array([])
    for j in range(ny):
        for i in range(nx):
            print("computing cell " + "x " + str(i) +", y " + str(j))
            
            for p in range(len(triangles)):
                             
                v1 = lc[triangles[p,0]]
                v2 = lc[triangles[p,1]]
                v3 = lc[triangles[p,2]]
                vertex = np.array([[v1[0], v1[1]],[v2[0], v2[1]],[v3[0], v3[1]]]) # transforming the triangle vertices into 2D
                if insideTriangle2D(np.array([X[i,j], Y[i,j]]), vertex):                         # Use barycentric weights to check intersection
                    
                    numberInter[i,j] +=1                      # count the triangles a point intersects
                    #interTri = np.append([p], interTri)     # save the index to triangle that intersects
                       # Save the Z coordinate of the intersection Point on the triangle
                       #test 
                    SecPoint =upperTest(v1,v2,v3,np.array([X[i,j], Y[i,j], Z[0]]))
                    InterPointZ = np.append([upperTest(v1,v2,v3,np.array([X[i,j], Y[i,j], Z[0]]))], InterPointZ)
                    #print(SecPoint)
                else:
                    continue


                
    InterPointZ = np.flip(InterPointZ)            
    print("intersection 2D completed")
     # used to count indeces for InterTri and InterPointZ



    for k in range(nz):
        triIdx = 0
        for j in range(ny):
            
                for i in range(nx):
                    #print("computing cell " + "x " + str(i) +", y " + str(j))
                    #print(triIdx)
                    
                    
                    if numberInter[i,j] == 0:
                        Phase[i,j,k]= 0
                    else:
                        intersections = 0  
                        for t in range(int(numberInter[i,j])):
                            if Z[k] <= InterPointZ[int(triIdx)]: # get the Z elevation of the intersection point and see 
                                intersections += 1
                            triIdx += 1
                        if intersections % 2 == 0: 
                            Phase[i, j,k] = 0
                        else:
                            Phase[i,j,k] = 1

    
    return Phase



from stl import mesh


def fastRayFile(inFile:str, grid):
    
    
    triangles = mesh.Mesh.from_file(inFile)
    vertices = triangles.points
    x,y,z = grid[0], grid[1], grid[2]
    X,Y = x[:,:,0], y[:,:,0]
    Z = z[0,0,:]
    Phase = np.zeros_like(x) # Phase array(comes notflattend)

    nx, ny, nz = Phase.shape[0],Phase.shape[1],Phase.shape[2]

    numberInter = np.zeros_like(X) # count the intersections on the x,y plane

    InterPointZ = np.array([]) # Z -  coordinate of intersectionpoints
    
    for j in range(ny):
        for i in range(nx):
            print("computing cell " + "x " + str(i) +", y " + str(j))
            for p in range(len(vertices)):
                             
                
                vertex = np.array([[vertices[p,0], vertices[p,1]],[vertices[p,3], vertices[p,4]],
                                   [vertices[p,6], vertices[p,7]]]) # transforming the triangle vertices into 2D
                if insideTriangle2D(np.array([X[i,j], Y[i,j]]), vertex):  # Use barycentric weights to check intersection
                    
                    numberInter[i,j] +=1                      # count the triangles a point intersects
                    #interTri = np.append([p], interTri)     # save the index to triangle that intersects
                       # Save the Z coordinate of the intersection Point on the triangle
                       #test 
                    SecPoint = upperTest(vertices[p,0:3],vertices[p,3:6],vertices[p,6:9],np.array([X[i,j], Y[i,j], Z[0]]))
                    InterPointZ = np.append([SecPoint], InterPointZ)
                    # InterPointZ = np.insert(InterPointZ, upperTest(vertices[p,0:3],vertices[p,3:6],vertices[p,6:9],np.array([X[i,j], Y[i,j], Z[0]])))
                    #InterPointZ = np.append([upperTest(vertices[p,0:3],vertices[p,3:6],vertices[p,6:9],np.array([X[i,j], Y[i,j], Z[0]]))], InterPointZ)
                    #print(SecPoint)
                else:
                    continue
    InterPointZ = np.flip(InterPointZ)            
    print("intersection 2D completed")
     # used to count indeces for InterTri and InterPointZ
    for k in range(nz):
        triIdx = 0
        for j in range(ny):
            
                for i in range(nx):
                    #print("computing cell " + "x " + str(i) +", y " + str(j))
                    #print(triIdx)
                    
                    
                    if numberInter[i,j] == 0:
                        Phase[i,j,k]= 0
                    else:
                        intersections = 0  
                        for t in range(int(numberInter[i,j])):
                            if Z[k] <= InterPointZ[int(triIdx)]: # get the Z elevation of the intersection point and see 
                                intersections += 1
                            triIdx += 1
                        if intersections % 2 == 0: 
                            Phase[i, j,k] = 0
                        else:
                            Phase[i,j,k] = 1
     
     
    """
    triIdx = 0
    for j in range(ny):
        
            for i in range(nx):
                #print("computing cell " + "x " + str(i) +", y " + str(j))
                #print(triIdx)
                if numberInter[i,j] == 0:
                    Phase[i,j,:]= 0
                    continue
                else:
                    for k in range(nz):
                        intersections = 0  
                        for t in range(int(numberInter[i,j])):
                            if Z[k] <= InterPointZ[int(triIdx)]: # get the Z elevation of the intersection point and see 
                                intersections += 1
                            #triIdx += 1
                        if intersections % 2 == 0: 
                            Phase[i, j,k] = 0
                        else:
                            Phase[i,j,k] = 1
                    triIdx += numberInter[i,j]
    """
                        

    
    return Phase




@jit(nopython=True)
def rayTrc_core(Phase,X,Y,Z,nx,ny,nz,nHit,zHit,vertices):

    for j in range(ny):
        for i in range(nx):
            print("computing cell " + "x " + str(i) +", y " + str(j))
            for p in range(len(vertices)):
                             
                
                vertex = np.array([[vertices[p,0], vertices[p,1]],[vertices[p,3], vertices[p,4]],
                                   [vertices[p,6], vertices[p,7]]]) # transforming the triangle vertices into 2D
                if insideTriangle2D(np.array([X[i,j], Y[i,j]]), vertex):  # Use barycentric weights to check intersection
                    
                    nHit[i,j] +=1                      # count the triangles a point intersects
                    #interTri = np.append([p], interTri)     # save the index to triangle that intersects
                       # Save the Z coordinate of the intersection Point on the triangle
                       #test 
                    SecPoint = upperTest(vertices[p,0:3],vertices[p,3:6],vertices[p,6:9],np.array([X[i,j], Y[i,j], Z[0]]))
                    zHit = np.append([SecPoint], zHit)

                else:
                    continue
    
    zHit = np.flip(zHit)            
    print("intersection 2D completed")
     # used to count indeces for InterTri and InterPointZ
    for k in range(nz):
        triIdx = 0
        for j in range(ny):
                for i in range(nx):
                    #print("computing cell " + "x " + str(i) +", y " + str(j))
                    #print(triIdx)
                    
                    
                    if nHit[i,j] == 0:
                        Phase[i,j,k]= 0
                    else:
                        intersections = 0  
                        for t in range(int(nHit[i,j])):
                            if Z[k] <= zHit[int(triIdx)]: # get the Z elevation of the intersection point and see 
                                intersections += 1
                            triIdx += 1
                        if intersections % 2 == 0: 
                            Phase[i, j,k] = 0
                        else:
                            Phase[i,j,k] = 1
                            
    return Phase


def rayTrc_rlgrd(inFile:str, rlgrd:rectLinGrid):
    
    # STL *surface* mesh
    triangles = mesh.Mesh.from_file(inFile)
    vertices = triangles.points
    
    # 2D slices of rectilinear grid
    X,Y,Z = rlgrd.X[:,:,0],rlgrd.Y[:,:,0],rlgrd.Z[0,0,:]
    
    #
    Phase = rlgrd.PHS
    nHit = np.zeros_like(X) # count the intersections on the x,y plane
    zHit = np.array([]) # Z -  coordinate of intersectionpoints   
    Phase = rayTrc_core(Phase,X,Y,Z,rlgrd.nx,rlgrd.ny,rlgrd.nz,nHit,zHit,vertices)
    
    return Phase


















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




def OpenVolumeTest(inFile, nInter, nPrec, grid):
    """
    surface triangulation for open surfaces
    """
    data = readSVG(inFile)
    path = data.Curves
    Layers, numLayers = getLayers(inFile)
    lc = getPointCoords(data, path, nInter,nPrec)

    
    nQuads = nInter+numLayers -1
    tri1 = np.array([])
    tri2 = np.array([])
    for p in range(nQuads):
        numP = int(len(lc)/(nQuads+1))
        nodeIdx = numP *p
        shape = (int(numP*2)-2,3)
        shapeS = (numP,3)
        ind = np.zeros(shapeS, dtype = int)
        ind2 = np.zeros(shapeS,dtype=int)
        arrshift = numP-1
        counter = 0
        if p == 0:
            for i in range(numP-1):
                if i == 0:
                    continue
                ind[i,0] = i
                ind[i,1] = i+1
                ind[i,2] = numP+i
                
                ind2[i,0] = numP +i
                ind2[i,1] = numP + i+1
                ind2[i,2] = i+1
                counter +=1
            tri1 = ind
            tri2 = ind2
        else:
    
            for i in range(numP-1):
                if i == 0:
                    continue
                ind[i,0] = nodeIdx + i
                ind[i,1] = nodeIdx + i+1
                ind[i,2] = nodeIdx + arrshift+i
                
                ind2[i,0] = nodeIdx + arrshift +i
                ind2[i,1] = nodeIdx + arrshift + i+1
                ind2[i,2] = nodeIdx+ i+1
                counter +=1
                
            tri2 = np.concatenate((ind2, tri2))    
            tri1 = np.concatenate((ind, tri1))
    
    triangles = np.concatenate((tri1,tri2))

 
    Phase = rayTracingRectilinear(triangles,lc,grid)
    return Phase

    



    
        
    