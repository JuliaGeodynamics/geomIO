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
#import ipdb
import scipy as sc
from scipy import interpolate

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
    upBound = np.amax(grid[2,:])
    cells = len(grid[0,:])
    Phase = np.zeros([cells])
    
            
    for p in range(len(triangles)):
        print(p)
        #currentTri = np.array([lc[triangles[p,0]],lc[triangles[p,1]],lc[triangles[p,2]]])
        v1 = lc[triangles[p,0]]
        v2 = lc[triangles[p,1]]
        v3 = lc[triangles[p,2]]
        for r in range(len(Phase)):
            top = grid[:,r].copy()
            top[2] = upBound
            nInter = 0
            if intersec(grid[:,r],top,v1,v2,v3):
                nInter +=1
                if nInter == 0:
                    Phase[r]= 0
                elif nInter %2 !=0:
                    
                    Phase[r]=1
            
            
        #first ind = first tri
        # mesh1 = np.concatenate((mesh1, normals[p]))
        # mesh1 = np.concatenate((mesh1, lc[triangles[p,0]]))
        # mesh1 = np.concatenate((mesh1, lc[triangles[p,1]]))
        # mesh1 = np.concatenate((mesh1, lc[triangles[p,2]]))
        
    #mesh1 = np.reshape(mesh1,(len(triangles)*4,3))
    #np.save("surf",mesh1)

    
    
    return Phase

# from scipy.spatial import ConvexHull, convex_hull_plot_2d,Delaunay 
    

# def layerTri(lc, numP, LLI, last = True):
    

    
#     if (numP % 2) == 0:
#         shape = (int((numP-2)/2),3)
#         index1 = np.zeros(shape, dtype = int)
#         index2 = np.zeros(shape, dtype = int)
#         for i in range(len(index1)):
#             index1[i,0] = i
#             index1[i,1] = i+1
#             index1[i,2] = numP-1-i
            
#             index2[i,0] = i
#             index2[i,1] = numP-1-i
#             index2[i,2] = numP -2 -i
    
#         triangles = np.concatenate((index1,index2))

#     else:
#         shape = (int((numP-3)/2),3)
#         index1 = np.zeros(shape, dtype = int)
#         index2 = np.zeros(shape, dtype = int)
#         for i in range(len(index1)):
#             index1[i,0] = i+1
#             index1[i,1] = i+2
#             index1[i,2] = numP-1-i
            
#             index2[i,0] = i+1
#             index2[i,1] = numP-i -1
#             index2[i,2] = numP -2 -i
    
#         first = np.array([[0,1,numP-1]])
#         triangles = np.concatenate((index1,index2))
#         triangles = np.concatenate((triangles, first))
        
#     normals = triNormals(triangles, lc)

#     mesh = np.array([])
#     for p in range(len(triangles)):
#         #first ind = first tri
#         mesh = np.concatenate((mesh, normals[p]))
#         mesh = np.concatenate((mesh, lc[triangles[p,0]]))
#         mesh = np.concatenate((mesh, lc[triangles[p,1]]))
#         mesh = np.concatenate((mesh, lc[triangles[p,2]]))
        
#     mesh = np.reshape(mesh,(len(triangles)*4,3))
   
#     if last and (numP % 2) == 0:
#         trueIdx = np.concatenate((index1 + LLI, index2 + LLI))
#         triangles = trueIdx
#     elif last:
#         trueIdx = np.concatenate((index1 + LLI, index2 + LLI))
#         triangles = np.concatenate((trueIdx, first + LLI))

        
            
#     return mesh, triangles

# def triangulateCover(lc):
#     numP = len(lc)
#     shape = len(lc)-2
#     ind = np.zeros((shape,3))
#     triangles = np.zeros((shape,3), dtype = int)
#     for i in range(shape):
#         triangles[i,0] = i 
#         triangles[i,1] = i+1
#         triangles[i,2] = numP-2-i
#     normals = triNormals(triangles, lc)
#     mesh = np.array([])
#     for p in range(len(triangles)):
#         #first ind = first tri
#         mesh = np.concatenate((mesh, normals[p]))
#         mesh = np.concatenate((mesh, lc[triangles[p,0]]))
#         mesh = np.concatenate((mesh, lc[triangles[p,1]]))
#         mesh = np.concatenate((mesh, lc[triangles[p,2]]))
        
#     mesh = np.reshape(mesh,(len(triangles)*4,3))
    
#     return  mesh, triangles

# def triSurfClose(inFile, nInter, nPrec):
#     """
#     surface triangulation for closed volumes
#     """
#     Layers, numLayers = getLayers(inFile)
#     lc = getCarthesian(inFile,nInter,nPrec)

#     nQuads = nInter+numLayers -1
#     tri1 = np.array([])
#     tri2 = np.array([])
    
#     #========Section where the "covers" are done : faulty
#     #number of points per layer
#     numP = int(len(lc)/(nQuads+1))
#     #get x and y coordinates from first and last layer
#     #cov1 = lc[0:numP,0:2]
#     cov1Z = lc[0:numP,:]
#     #last layer index, start of last layer
#     LLI = len(lc)-numP
#     #cov2 = lc[LLI:-1,0:2]
#     cov2Z =lc[LLI:-1,:]
#     dummy = np.array([lc[len(lc)-1]])
#     cov2Z = np.concatenate((cov2Z,dummy))
#     # startCov = Delaunay(cov1,furthest_site=True)
#     # endCov = Delaunay(cov2, furthest_site=True)
    
#     #cCF, cTF = layerTri(cov1Z,numP,LLI, False)
#     #cCL, cTL = layerTri(cov2Z, numP, LLI, True)
#     cCF, cTF = triangulateCover(cov1Z)
#     cCL, cTL = triangulateCover(cov2Z)
#     cTL += LLI
#     #coordinatses of both covers
#     cover = np.concatenate((cCF,cCL))
#     #index of the triangles in lc
#     covTri = np.concatenate((cTF,cTL))

    

#     #=======section for inside: working    

#     for p in range(nQuads):

#         nodeIdx = numP *p
#         shape = (int(numP*2)-2,3)
#         shapeS = (numP,3)
#         ind = np.zeros(shapeS, dtype = int)
#         ind2 = np.zeros(shapeS,dtype=int)
#         arrshift = numP-1
#         counter = 0
#         if p == 0:
#             for i in range(numP):

#                 ind[i,0] = i
#                 ind[i,1] = i+1
#                 ind[i,2] = numP+i
                
#                 ind2[i,0] = numP +i
#                 ind2[i,1] = numP + i+1
#                 ind2[i,2] = i+1
#                 counter +=1
#             tri1 = ind
#             tri2 = ind2
#         else:
    
#             for i in range(numP):

#                 ind[i,0] = nodeIdx + i
#                 ind[i,1] = nodeIdx + i+1
#                 ind[i,2] = nodeIdx + arrshift+i
                
#                 ind2[i,0] = nodeIdx + arrshift +i
#                 ind2[i,1] = nodeIdx + arrshift + i+1
#                 ind2[i,2] = nodeIdx+ i+1
#                 counter +=1
                
#             tri2 = np.concatenate((ind2, tri2))    
#             tri1 = np.concatenate((ind, tri1))
    
#     triangles = np.concatenate((tri1,tri2))
#     normals = triNormals(triangles, lc)

#     mesh1 = np.array([])
#     #print(counter)
    
    
            
#     for p in range(len(triangles)):
#         #first ind = first tri
#         mesh1 = np.concatenate((mesh1, normals[p]))
#         mesh1 = np.concatenate((mesh1, lc[triangles[p,0]]))
#         mesh1 = np.concatenate((mesh1, lc[triangles[p,1]]))
#         mesh1 = np.concatenate((mesh1, lc[triangles[p,2]]))
#         #vertex = Tri(lc[triangles[p,0]],lc[triangles[p,1]],lc[triangles[p,1]])
        
#     mesh1 = np.reshape(mesh1,(len(triangles)*4,3))
#     mesh1 = np.concatenate((mesh1,cover))
#     triangles = np.concatenate((triangles,covTri))
#     #for d in range(2):
#     # mesh1 = np.concatenate((mesh1,mesh2))
#     # mesh1 = np.concatenate((mesh1,mesh3))

#     # triangles = np.concatenate((triangles,triCov1))
#     # faceLast = triCov2+LLI
#     # triangles = np.concatenate((triangles,faceLast))        
    
#     #return cover, covTri
#     return mesh1, triangles

# # inFile = 'input/slab.svg'

# # mesh1, tri = triSurfClose(inFile, 5,20)

    
        
    