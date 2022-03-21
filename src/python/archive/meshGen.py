# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 16:11:52 2022

@author: lumoser
"""

import numpy as np

import sys,os
#archive

def createMesh(nx : int ,ny: int, nz: int, xBound:tuple, yBound:tuple, zBound:tuple):
    
    """
    function to generate a 3D meshgrid for testing the raytraycing
    nx,ny,nz = number of nodes in each direction
    xBound, yBound, zBound = lower and upper boundary in each direction
    returns: x,y,z = list of coordinates along each axis
    X,Y,Z = the mesh
    
    """
    
    x = np.linspace(xBound[0],xBound[1],nx)
    y = np.linspace(yBound[0],yBound[1],ny)
    z = np.linspace(zBound[0],zBound[1],nz)
    
    X,Y,Z = np.meshgrid(x,y,z, indexing ='ij') #matrix indexing{ij} array[i,j]
    # X = np.transpose(X, (1,0,2)) 
    # Y = np.transpose(Y, (1,0,2))
    # Z = np.transpose(Z, (1,0,2))
    
    return x,y,z, X,Y,Z

#x,y,z,X,Y,Z = testMesh(4,6,8, [0,250],[100,200],[-200,0])
    
def writeVTK(x,y,z, Phase, name:str):
    """
    function to generate vtk-files from rectilinear grids
    input : x,y,z = the list of coordinates along each axis
    Phase: array of shape(nx,ny,nz) with the data
    """
    from uvw import RectilinearGrid, DataArray
    
    print("------------------------")
    print("Writing vtk file to output/"+ str(name)+".vtr")
    print("------------------------")
    grid = RectilinearGrid(str(name)+".vtr", (x, y, z), compression=True)
    grid.addPointData(DataArray(Phase, range(3), 'Phase'))
    grid.write()
    
    return