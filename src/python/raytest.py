#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 17:03:02 2021

@author: lucas
"""
from time import time
import numpy as np
import matplotlib.path as mpltPath

tri = np.load("triS.npy")
grid = np.load("gridVals.npy")*1e3
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
            if intersec(grid[:,i],top,v1,v2,v3):
                Phase[i]=1
            cS +=4
            cE +=4
    return Phase

Ph = gridWrap(tri,grid)
Ph = np.reshape(Ph,(16,16,16),order='C')

    
        
    