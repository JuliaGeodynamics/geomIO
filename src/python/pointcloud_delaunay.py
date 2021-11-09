import ipdb
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import vtk
from vtk.util import numpy_support as VN


def findAlpha(PointCloudArray):
    alpha = 0    
    return alpha


def delny3D(PointCloudArray):
    ## generate VTK pointcloud
    PC_points = vtk.vtkPoints()
    PC_vertices = vtk.vtkCellArray()
    
    for i in range(0,len(PointCloudArray)):
    
        Point = PointCloudArray[i,:]
    
        id1 = PC_points.InsertNextPoint(Point)
        PC_vertices.InsertNextCell(1)
        PC_vertices.InsertCellPoint(id1)
    
    PC_polydata = vtk.vtkPolyData()
    PC_polydata.SetPoints(PC_points)
    PC_polydata.SetVerts(PC_vertices)
    
    
    ## generate 3D delaunay
    delny = vtk.vtkDelaunay3D()
    delny.SetInputData(PC_polydata)
    delny.SetOffset(4.0)
    delny.SetAlpha(20)
    delny.Update()
    
    delnyPolyData = delny.GetOutput()
    return delnyPolyData
