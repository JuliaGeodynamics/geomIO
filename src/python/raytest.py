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

from stl import mesh



# A utility function to calculate area
# of triangle formed by (x1, y1),
# (x2, y2) and (x3, y3)
@jit(nopython=True) 
def area(x1, y1, x2, y2, x3, y3):
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1)
                + x3 * (y1 - y2)) / 2.0)
 
 
# A function to check whether point P(x, y)
# lies inside the triangle formed by
# A(x1, y1), B(x2, y2) and C(x3, y3)
@jit(nopython=True) 
def isInside(x1, y1, x2, y2, x3, y3, x, y):

    if x1==x2 and x1==x3 and y1==y2 and y1==y3:
        return False
 
    # Calculate area of triangle ABC
    A = area (x1, y1, x2, y2, x3, y3)
    # Calculate area of triangle PBC
    A1 = area (x, y, x2, y2, x3, y3)
    # Calculate area of triangle PAC
    A2 = area (x1, y1, x, y, x3, y3)
    # Calculate area of triangle PAB
    A3 = area (x1, y1, x2, y2, x, y)
    # Check if sum of A1, A2 and A3
    # is same as A
    
    if(abs(A - (A1 + A2 + A3)) < 1e-4):
        return True
    else:
        return False
 
# Driver program to test above function
# Let us check whether the point P(10, 15)
# lies inside the triangle formed by
# A(0, 0), B(20, 0) and C(10, 30)
#if (isInside(0, 0, 20, 0, 10, 30, 10, 15)):
#    print('Inside')
#else:
#    print('Not Inside')





@jit(nopython=True)
def createPlane(p:np.array,q:np.array,r:np.array):
    """
    Creates a plane in Coordinate form
    ax + by +cz = d
    where a,b,c are the components of the normal array and d is d
    """    
    normal = np.cross(q-p,r-p)
    d = np.vdot(normal, p) # This is still causes a NUMBA performance warning
    return normal, d

@jit(nopython=True)
def rayXplane(v1:np.array,v2:np.array,v3:np.array, zVal:np.array):
    Coor, d = createPlane(v1,v2,v3)
    r = ((d - zVal[0]*Coor[0]-zVal[1]*Coor[1])/Coor[2])-zVal[2]
    InterPoint = zVal+r*np.array([0,0,1])
    return InterPoint[2]


@jit(nopython=True)
def rayTrc_core(x,y,z,vertices,zHit,verbose):
    #x,y,z = grid[0], grid[1], grid[2]
    X,Y = x[:,:,0], y[:,:,0]
    Z = z[0,0,:]
    Phase = np.zeros_like(x) # Phase array(comes notflattend)
    nx, ny, nz = Phase.shape[0],Phase.shape[1],Phase.shape[2]
    nHit = np.zeros_like(X,dtype=np.int32) # count the intersections on the x,y plane
    #zHit = np.array([]) # Z -  coordinate of intersectionpoints
    
    for j in range(ny):
        for i in range(nx):
            if verbose:
                print("computing cell " + "x " + str(i) +", y " + str(j))
            for p in range(len(vertices)):
                
                if isInside(vertices[p,0],vertices[p,1], vertices[p,3], vertices[p,4], vertices[p,6], vertices[p,7],  X[i,j], Y[i,j]):                    
                    nHit[i,j] +=1                      # count the triangles a point intersects
                    SecPoint = rayXplane(vertices[p,0:3],vertices[p,3:6],vertices[p,6:9],np.array([X[i,j], Y[i,j], Z[0]]))
                    zHit = np.append([SecPoint], zHit)
                else:
                    continue
    zHit = np.flip(zHit)            
    print("intersection 2D completed")
     # used to count indices for InterTri and zHit
    for k in range(nz):
        triIdx = 0
        for j in range(ny):
                for i in range(nx):
                    #print("computing cell " + "x " + str(i) +", y " + str(j))

                    if nHit[i,j] == 0:
                        Phase[i,j,k]= 0
                    else:
                        intersections = 0  
                        for t in range(nHit[i,j]):
                            if Z[k] <= zHit[triIdx]: # get the Z elevation of the intersection point and see 
                                intersections += 1
                            triIdx += 1
                        if intersections % 2 == 0: 
                            Phase[i,j,k] = 0
                        else:
                            Phase[i,j,k] = 1

    return Phase


## ALternative 
# sort zHit[i,j,:] vertically
#def rayTrc_core(x,y,z,vertices):
#    #x,y,z = grid[0], grid[1], grid[2]
#    X,Y = x[:,:,0], y[:,:,0]
#    Z = z[0,0,:]
#    Phase = np.zeros_like(x,dtype=np.float32) # Phase array(comes notflattend)
#    nx, ny, nz = Phase.shape[0],Phase.shape[1],Phase.shape[2]
#    nHit = np.zeros_like(X,dtype=int) # count the intersections on the x,y plane
#    zHit = np.ones_like(x,dtype=np.float64) * Z[nz-1] # Z -  coordinate of intersectionpoints
#
#    for j in range(ny):
#        for i in range(nx):
#            print("computing cell " + "x " + str(i) +", y " + str(j))
#            for p in range(len(vertices)):
#                
#                if isInside(vertices[p,0],vertices[p,1], vertices[p,3], vertices[p,4], vertices[p,6], vertices[p,7],  X[i,j], Y[i,j]):                    
#                    nHit[i,j] +=1                      # count the triangles a point intersects
#                    zHit[i,j,nHit[i,j]-1] = rayXplane(vertices[p,0:3],vertices[p,3:6],vertices[p,6:9],np.array([X[i,j], Y[i,j], Z[0]]))                   
#                else:
#                    continue   
#    print("intersection 2D completed")
#     # used to count indices for InterTri and zHit
#    for j in range(ny):
#        for i in range(nx):
#            if nHit[i,j] != 0:
#                hitCount = 0
#                zTest = np.sort(zHit[i,j,:])   
#                for k in range(nz):
#                    #if Z[k] > zHit[i,j,hitCount]:
#                    #if Z[k] > zTest[hitCount]: 
#                        hitCount += 1
#                    if hitCount % 2 == 0:
#                        Phase[i,j,k] = 1
#                    else:
#                        Phase[i,j,k] = 0
#    
#    return Phase

def rayTrc_rlgrd(inFile:str, grid):
       
    triangles = mesh.Mesh.from_file(inFile)
    vertices = triangles.points.astype(np.float64)
    x,y,z = grid[0], grid[1], grid[2]
    # X,Y = x[:,:,0], y[:,:,0]
    # Z = z[0,0,:]
    # Phase = np.zeros_like(x) # Phase array(comes notflattend)
    # nx, ny, nz = Phase.shape[0],Phase.shape[1],Phase.shape[2]
    # nHit = np.zeros_like(X) # count the intersections on the x,y plane
    zHit = np.array([]) # Z -  coordinate of intersectionpoints

    Phase = rayTrc_core(x,y,z, vertices,zHit,verbose=True)    
    #Phase = rayTrc_core(x,y,z, vertices)    # Alternative
   
    return Phase
 
 
 






