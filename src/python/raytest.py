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


def sortGrid(grid):
    
    
    
    return

def fastRay(triangles, lc, grid):
    
    
    
    Phase = np.array([])
    zVals = np.unique(grid[2,:]) 
    Layers = len(np.unique(grid[2,:]))                     # number of different z coordinates
    numPoints = np.unique(grid[2,:], return_counts = True)[1][0]  # number of points per layer
    # I need every 56th(layers) x and y coordinate better have function that determines this
    x = np.array([])
    y = np.array([])
    #count = 0
    for s in range(int(len(grid[0,:])/Layers)):
        x = np.append([x],grid[0,s*Layers])
        y = np.append([y],grid[0,s*Layers])

    points = np.concatenate((x,y))
    #points = np.array([[grid[0,0:Layers], grid[1,0:Layers]]])     # x and y coordinates of these cells
    numPoints = points.size                                 # number of cells per z coordinate
    #inter = np.zeros_like(points)                          # number of intersections with trianglew
    numberInter = np.zeros_like(points)                     # number of intersections with triangles
    interTri = np.array([])                                 # index wich traingles are inntersected
    
    for i in range(len(points)):                            # check all x  and y coordinates for intersection
        for p in range(len(triangles)):                 
                v1 = lc[triangles[p,0]]
                v2 = lc[triangles[p,1]]
                v3 = lc[triangles[p,2]]
                vertex = np.array([[v1[0], v1[1]],[v2[0], v2[1]],[v3[0], v3[1]]])
                if insideTriangle2D:
                    #inter[i]= 1
                    numberInter[i] +=1                      # count the intersections
                    interTri = np.append([p], interTri)     # save the index to triangle that intersects
                else:
                    continue
    triCounter = 0

    ph = np.zeros(numPoints)                             # test for all coordinates
    for s in range(numPoints):
         
       
            for p in range(numZ):
                if numberInter[s] ==0:                          # if no intersection whole column will be ignored
            
                    Phase = np.append([0],Phase)              # das wird ned funktioneren, eher concat
                else:
                    intersections = 0
                    for t in range(numberInter[s]):
                        triIndex = int(interTri[t+ triCounter])
                        v1 = lc[triangles[triIndex,0]]
                        v2 = lc[triangles[triIndex,1]]
                        v3 = lc[triangles[triIndex,2]]
                        v1 = v1[2]
                        v2 = v2[2]
                        v3 = v3[2]
                        triCounter +=1
                        if interTRiangleFast(v1,v2,v3, zVals[p]) == True:
                            intersections += 1
                        else:
                            continue
                    if intersections % 2 == 0:
                        Phase = np.append([0], Phase)
                    else:
                        Phase = np.append([1],Phase)
                        
                        print("ich stehe nur hier damit der code kompiliert ;)" )
    return Phase



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

#Ph = gridWrap(tri,grid)
#Ph = np.reshape(Ph,(16,16,16),order='C')


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
from numba import jit

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


def OpenVolumeTest(inFile, nInter, nPrec, grid):
    """
    surface triangulation for open surfaces
    """
    Layers, numLayers = getLayers(inFile)
    lc = getCarthesian(inFile,nInter,nPrec)

    
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

    #Phase = np.zeros((grid.size/3))

    #redo in parallel
    print(len(triangles))
    
            
            
        #first ind = first tri
        # mesh1 = np.concatenate((mesh1, normals[p]))
        # mesh1 = np.concatenate((mesh1, lc[triangles[p,0]]))
        # mesh1 = np.concatenate((mesh1, lc[triangles[p,1]]))
        # mesh1 = np.concatenate((mesh1, lc[triangles[p,2]]))
        
    #mesh1 = np.reshape(mesh1,(len(triangles)*4,3))
    #np.save("surf",mesh1)
    #Phase = mainTest(triangles, lc, grid)
    Phase = fastRay(triangles, lc, grid)
    
    return Phase

    



    
        
    