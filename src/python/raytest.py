#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 17:03:02 2021

@author: lucas
"""
from time import time
import numpy as np
import matplotlib.path as mpltPath

#tri = np.load("triS.npy")
#grid = np.load("gridVals.npy")*1e3
#surf = np.delete(surface,np.arange(surface.size,3))
#surf = np.reshape(surf, (1728,3))



def line2Pl(p1,p2):
    line = np.array([[p1[0],p1[1],p1[2],1],
                     [p2[0],p2[1],p2[2],1]])
    pluecker = np.array([
        p1[0]*p2[1]-p2[0]*p1[1],
        p1[0]*p2[2]-p2[0]*p1[2],
        p1[0]-p2[0],
        p1[1]*p2[2]- p2[1]*p1[2],
        p1[2]-p2[2],
        p2[1]-p1[1]])
    return pluecker

def sideOP(a, b):
    operator = a[0]*b[4]+a[1]*b[5]+a[2]*b[3]+a[3]*b[2]+a[4]*b[0]+a[5]*b[1]
    return operator

def intersec(p1,p2,v1,v2,v3, eps = 0.1):
    """
    Uses pluecker Coordinate system to do determine line triangle intersection
    currently unused
    Documentation: https://members.loria.fr/SLazard/ARC-Visi3D/Pant-project/files/Line_Triangle.html
    """
    
    line = line2Pl(p1,p2)
    l1 = line2Pl(v1,v2)
    l2 = line2Pl(v2,v3)
    l3 = line2Pl(v3,v1)
    s1= sideOP(line,l1)
    s2 = sideOP(line,l2)
    s3 = sideOP(line,l3)
    if s1 and s2 and s3 != 0:
        
        if s1 <0 and s2 <0 and s3 <0:
            return True
        elif s1 >0 and s2 >0 and s3 >0:
            return True
        #elif (-eps)< s1< eps and (-eps)< s2< eps  or (-eps)< s2< eps  and (-eps)< s3< eps  or (-eps)< s13 eps  and (-eps)< s1< eps :
         #   return True
        else: 
            return False
    else:
        
        #here comes the coplanar stuff
        return

def insideTriangle2D(point, vertices):
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
    point = point.astype(float)
    vertices = vertices.astype(float)
    a = vertices[0]
    b = vertices[1]
    c = vertices[2]
    xd = np.cross(a,b) + np.cross(b,c) + np.cross(c,a)
    if xd ==0:
        xd+=0.00001
    
    xa = np.cross(b,c) + np.cross(point, b-c)
    xb = np.cross(c,a) + np.cross(point, c-a)
    xc = np.cross(a,b) + np.cross(point, a-b)
    
    wa = xa/xd
    wb = xb/xd
    wc = xc/xd
    
    
    if 0<wa<1 and 0<wb<1 and 0<wc<1:
        return True
    else:
        return False
    
    
    #this is a point
    # v0 = vertices[0]
    # #these are vectors
    # v1 = vertices[1]- v0
    # v2 = vertices[2] - v0
    
    # a = (np.cross(point,v2)-np.cross(v0,v2))/np.cross(v1,v2)
    
    # b = (np.cross(point,v1)-np.cross(v0,v1))/np.cross(v1,v2)
    
    # c = a+b
    # if c<1 and a >0 and b>0:
    #     return True
    # else:
    #     return False
#---------test for inside 2d triangle(passed)--------    
# ver = np.array([[1,1],[4,2],[2,7]])
# point = np.array([2,3])
# x = insideTriangle2D(point,ver)

def interTRiangleFast(v1,v2,v3, zVal):
    
    """
    https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution    
    
    """
    avg = (v1+v2+v3)/3
    if zVal < avg :
        return True
    else:
        return False
    
def createPlane(p,q,r):
    """
    Creates a plane in Coordinate form
    ax + by +cz = d
    where a,b,c are the components of the normal array and d is d
    """    
    PQ = q-p
    PR = r-p
    OP = p
    a = np.asanyarray([PQ,PR,OP])
    b = np.array([2,2,2])
    normal = np.cross(PQ,PR)
    d = np.dot(normal, OP)
    
    
    
    return normal, d

def upperTest(v1,v2,v3, zVal):
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
     #test with more intersection(flip current file by 90°)
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

    numberInter = np.zeros_like(X)
    #interTri = np.array([])
    InterPointZ = np.array([])
    
    for j in range(ny):
        for i in range(nx):
            print("computing cell " + "x " + str(i) +", y " + str(j))
            for p in range(len(vertices)):
                             
                
                vertex = np.array([[vertices[p,0], vertices[p,1]],[vertices[p,3], vertices[p,4]],
                                   [vertices[p,6], vertices[p,7]]]) # transforming the triangle vertices into 2D
                if insideTriangle2D(np.array([X[i,j], Y[i,j]]), vertex):                         # Use barycentric weights to check intersection
                    
                    numberInter[i,j] +=1                      # count the triangles a point intersects
                    #interTri = np.append([p], interTri)     # save the index to triangle that intersects
                       # Save the Z coordinate of the intersection Point on the triangle
                       #test 
                    #SecPoint =upperTest(v1,v2,v3,np.array([X[i,j], Y[i,j], Z[0]]))
                    InterPointZ = np.append([upperTest(vertices[p,0:3],vertices[p,3:6],vertices[p,6:9],np.array([X[i,j], Y[i,j], Z[0]]))], InterPointZ)
                    #print(SecPoint)
                else:
                    continue
    InterPointZ = np.flip(InterPointZ)            
    print("intersection 2D completed")
     # used to count indeces for InterTri and InterPointZ
     #test with more intersection(flip current file by 90°)
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

# def fastRayFile(inFile:str, grid):
    
#     zVals, lowerBound, Layers, numPoints, points = sortGrid(grid)
#     triangles = mesh.Mesh.from_file(inFile)


#     numPoints = len(points)                               # number of cells per Layer
#     numberInter = np.zeros(numPoints)                     # number of intersections with triangles
#     interTri = np.array([])                               # index wich traingles are inntersected
#     InterPointZ = np.array([])                            # Z coordinate of 3D intersection Point     
    
#     for i in range(numPoints):
#         print("Computing 2D for column" + str(i))         # check all x  and y coordinates for intersection
#         for p in range(len(triangles)):                 
#                 v1 = triangles.v0[p,:]
#                 v2 = triangles.v1[p,:]
#                 v3 = triangles.v2[p,:]
#                 vertex = np.array([[v1[0], v1[1]],[v2[0], v2[1]],[v3[0], v3[1]]]) # transforming the triangle vertices into 2D
#                 if insideTriangle2D(points[i,:], vertex):                         # Use barycentric weights to check intersection
                    
#                     numberInter[i] +=1                      # count the triangles a point intersects
#                     interTri = np.append([p], interTri)     # save the index to triangle that intersects
#                     # Save the Z coordinate of the intersection Point on the triangle
#                     InterPointZ = np.append([upperTest(v1,v2,v3,np.array([points[i,0],points[i,1], lowerBound]))], InterPointZ)
#                 else:
#                     continue
#     print("intersection 2D completed")
#     triCounter = 0 # used to count indeces for InterTri and InterPointZ
#     Phase = np.array([]) # Phase array(comes flattend)
    

                            
#     for s in range(numPoints):  # test for all coordinates
               
#             print("computing for column" + str(s))
#             for p in range(Layers):
                
#                 if numberInter[s] ==0:                        # if no intersection whole column will be ignored
            
#                     Phase = np.append([0],Phase)              # Add zeros for every cell with no intersection
#                 else:
#                     intersections = 0                         # Count number of total intersections
#                     for t in range(int(numberInter[s])):
#                         # triIndex = int(interTri[triCounter])  # Call the triangle that intersects
#                         # v1 = lc[triangles[triIndex,0]]        # Get its coordinates
#                         # v2 = lc[triangles[triIndex,1]]        # this is completely redundant
#                         # v3 = lc[triangles[triIndex,2]]
#                         # v1 = v1[2]
#                         # v2 = v2[2]
#                         # v3 = v3[2]
#                         if zVals[p] <= InterPointZ[triCounter]: # get the Z elevation of the intersection point and see 
#                             intersections += 1                  # wether is underneath the surface or not
#                         else:
#                             continue
#                         # if interTRiangleFast(v1,v2,v3, zVals[p]) == True:
#                         #     intersections += 1
#                         # else:
#                         #     continue hier war btw der fehler weil hier kein counter lief
#                     if intersections % 2 == 0:                  # if even number of intersections point is above surface 
#                         Phase = np.append([0], Phase)
#                     else:
#                         Phase = np.append([1],Phase)            # else it is underneath
#             triCounter +=int(numberInter[s])                    # update indices
                    

                    
#     return Phase
    
    


#normalvektoren weg



        #even uneven!!!
def gridWrap(surf, grid):
    upBound = np.amax(grid[2,:])
    cells = len(grid[0,:])
    Phase = np.zeros([cells])
    for i in range(cells):
        top = grid[:,i].copy()
        top[2] = upBound
        cS = 1
        cE = 4
        print(i)
        for p in range(int(len(surf)/4)):
            triangle = surf[cS:cE,:]
            v1=triangle[0]
            v2=triangle[1]
            v3=triangle[2]
            
            #dreiecke initialisieren
            nInter = 0
            if intersec(grid[:,i],top,v1,v2,v3):
                nInter +=1
                if nInter == 0:
                    Phase[i]= 0
                elif nInter %2 !=0:
                    
                    Phase[i]=1
            cS +=4
            cE +=4
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
#from numba import jit

class Tri(object):
    
    def __init__(self, v1,v2,v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        
    @classmethod 
    def ray(cls, coor1, coor2):
        if intersec(coor1, coor2, v1,v2,v3):
            return True
        else:
            return False
        
        
#@jit(nopython=True) 
def mainTest(triangles, lc, grid):
    upBound = np.amax(grid[2,:])
    cells = len(grid[0,:])
    Phase = np.zeros([cells])
    
    
    for p in range(len(triangles)):
        #loop over all triangles
        print(p)
        #currentTri = np.array([lc[triangles[p,0]],lc[triangles[p,1]],lc[triangles[p,2]]])
        v1 = lc[triangles[p,0]]
        v2 = lc[triangles[p,1]]
        v3 = lc[triangles[p,2]]
        for r in range(len(Phase)):
            #loop over coordinates
            top = grid[:,r].copy()
            top[2] = upBound
            nInter = 0
            
            if intersec(grid[:,r],top,v1,v2,v3):
                nInter +=1
                if nInter == 0:
                    Phase[r]= 0
                elif nInter %2 !=0:
                    
                    Phase[r]=1
    return Phase

def fileTest(inFile, grid):
    
    
    return

def OpenVolumeTest(inFile, nInter, nPrec, grid):
    """
    surface triangulation for open surfaces
    """
    data = readSVG(inFile)
    path = data.Curves
    Layers, numLayers = getLayers(inFile)
    lc = getCarthesian(data, path, nInter,nPrec)

    
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

    



    
        
    