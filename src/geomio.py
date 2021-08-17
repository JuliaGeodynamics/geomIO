# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# import numpy as np
# import sys, os

# from svglib.svglib import svg2rlg
# from reportlab.graphics import renderPDF, renderPM
# #import svgpathtools as pt
# import math 





from svgpathtools import svg2paths, real, imag, Line

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as Pt
import matplotlib.patches as patches
import math
import sys,os
import ipdb
import scipy as sc
from scipy import interpolate
#import stlWrite
from pointcloud_delaunay import *
import meshwrite as wr
os.chdir("..")



#redo for naming layers

def getLayers(inFile):
    """
    Parameters
    ----------
    inFile : Input svg File
        The .svg file used. Warning: must be Inkscape SVG

    Returns
    -------
    Layers : dict
        Dictonary containing the information which path belongs to which layer.
        If Layers represent geological ages then the first entry of the dict is 
        the first layer created in Inkscape(and usually the oldest)
    """
    f = open(inFile)
    text = f.readlines()
    index = []
    numLayers = 0
    for i in range(len(text)):
        if "<g" in text[i]:
            numLayers += 1
            index.append(i)
            
    index.append(int(len(text)))
    #print(index )
    #Layers = list()
    Layers = dict()        
    #print(index)
    initNum = 0
    compNum = 0
    
    #-------------BUG!!!------------------
    for i in range(len(index)-1):
        #print(i)#
        if i > 0:
            if compNum != initNum:
                #print(initNum,compNum)
                sys.exit("Number of Paths must be equal per Layer. Invalid number in:" + str(layer))
        compNum = 0
        for p in range(index[i],index[i+1]):
            if "inkscape:label"  in text[p]:
                layerSTR = text[p]
                layerSTR = layerSTR.split("\"")
                #print(layerSTR)
                layer = layerSTR[1]
                #print(layer)
            else:
                layerSTR = "Layer" + str(p)
            if "id=\"path"  in text[p]:
                
                if i == 0:
                    initNum += 1
                    compNum +=1
                    pathSTR = text[p]
                    pathSTR = pathSTR.split("\"")
                    pathSTR = pathSTR[1]
                    Layers[pathSTR] = layer
                else:                    
                    compNum +=1
                    #print (compNum)
                    pathSTR = text[p]
                    pathSTR = pathSTR.split("\"")
                    pathSTR = pathSTR[1]
                    Layers[pathSTR] = layer

        #Layers.append(Paths)
    return Layers, numLayers

def plot_line(inp, tol):
    #depracated
    
    line = np.array(inp)#*tol
     
    #plt.plot(line[:, 0], line[:, 1], color='b')
    
    plt.scatter(line[:, 0], line[:, 1], s=5)




def checkPath(attributes, Layers):
    """
    

    Parameters
    ----------
    Attributes and Layers from svgpathtools

    Returns
    -------
    number of paths per layer and also checks wether 
    the number of paths per layer is the same for all
    layers
    """
    
    #Check for dicts
    if not isinstance(attributes, list):
        sys.exit("Attributes must be list")
    if not isinstance(Layers, dict):
        sys.exit("Layers must be dictonary")
    numPaths = 0
    Counter = list(Layers.values())
    for i in Layers.values():
        if i in Counter[0]:
            numPaths +=1
    checkPath = 0
    for p in Layers.values():
        checkPath = 0
        for r in range(len(Counter)):
            if p in Counter[r]:
                checkPath +=1
        if checkPath != numPaths:
            sys.exit("Invalid number of Paths in " + str(p))
        #if i in attributes.values():
        #    print("dasads")
    return numPaths

def line2coor(path, tol = 1e-6):
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
    


class Bezier(object):
    """
    basic class for quadratic bezier curves
    
    """

    def __init__(self, start,c1,c2, end):
        self.start = start
        self.c1 = c1
        self.c2 = c2
        self.end = end
        
        return   
    
def controlPoints(line):
    """
    

    Parameters
    ----------
    line : line of svgpathtools module

    Returns
    -------
    cPoints: Class storing the controlpoints for quadratic bezier curves

    """
    cPoints = list()

    
    for i in range(len(line)):
        
        points = line._segments[i]
        start = points.start#.copy()
        c1 = points.control1#.copy()
        c2 = points.control2#.copy()
        end = points.end#.copy()
        start = np.array([real(start),-imag(start)])
        c1 = np.array([real(c1),-imag(c1)])
        c2 = np.array([real(c2),-imag(c2)])
        end = np.array([real(end),-imag(end)])
        coors = Bezier(start,c1,c2, end)
        #ipdb.set_trace()
        cPoints.append(coors)
        del coors
    return cPoints


def sortLayers(inFile):
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

def getCpoints(inFile):
    """
    This function generates a list of quadratic bezier Controlpoints
    stored as bezier class for each layer and puts it toghter in a list

    Parameters
    ----------
    inFile : inputfile

    Returns
    -------
    cPoints : List of bezier segments per layer

    """
    paths, attributes = svg2paths(inFile)
    cPoints = list()
    for i in range(len(paths)):
        line = paths[i]
        points = controlPoints(line)
        cPoints.append(points)
    
    return cPoints
    
def makeBezier(cPoints):
    """
    
    """
    #if len(cPoints) != 4:
     #   sys.exit("bezier Curves must be cubic")
    #bezierList = list()
    #for i in range(numInter):
    start = cPoints[:,0]
    c1 = cPoints[:,1]
    c2 = cPoints[:,2]
    end = cPoints[:,3]
    Points = Bezier(start,c1,c2,end)
        
    return Points

def getZvalues(inFile):
    Layers, numLayers = getLayers(inFile)
    names = list()
    for z in Layers.values():
        names.append(z)
    zCoors = np.array([])
    zCoor = np.array([], dtype = float)    
    for i in range(len(names)):
        dig = []
        if "p" in names[i]:
            dig = names[i].split("p")

            if "m" in dig[0]:
                dig[0] = dig[0].replace("m", "")
                dig[0] = -int(dig[0])
                #print(dig[0])
                power = int(len(dig[1]))
                num = int(dig[0])- (int(dig[1])*(10**-(int(len(dig[1])))))
            else:
                power = int(len(dig[1]))
                num = int(dig[0])+ (int(dig[1])*(10**-(int(len(dig[1])))))
                
            #print(num)
        else:
            if "m" in names[i]:
                num = names[i].replace("m","")
                num = -int(num)
            else:
                num = int(names[i])
                    

        zCoor = np.append([num], zCoor) 

    return zCoor
#Interp1D not necessary

def interZlayers(zCoor,numLayers):
    #print(zCoor)
    Inter = zCoor[-1]- zCoor[0]
    #print(Inter)
    Inter = Inter/numLayers
    zValues = np.linspace(zCoor[0] +Inter, zCoor[-1]- Inter, numLayers)
    #print(zValues)
    return zValues

#es ist halt wirklich die interpolation die von anfang bis end inter aber in der mitte ignoriert

def compEmpty(zCoor, numLayers):
    zValues = interZlayers(zCoor, numLayers)
    zCoor = np.append([zValues], zCoor)
    vals = np.sort(zCoor)
    #print(vals)
    idx = np.zeros(len(vals))
    for i in range(len(idx)):
        if vals[i] in zValues and vals[i] != vals[i-1]:
            idx[i] = 1
            
    return idx, vals
#write getcoords
#wrpng starting points, need to loop over starting points

#this function ignores all medium layers somehow, gotta interp from layer to layer
#another bug, doesnt work with only two layers? :O ----> fixed, was @interwrap
def interp(cPoints, zCoor, numPoints = 5):
    """
    This is the 1D interpolation in y direction

    Parameters
    ----------
    cPoints : List of controlpoints for all zusammengehörigen 
    bezierabschnitten
    numPoints : Number of interpolationpoints
         The default is 5.

    Returns
    -------
    Volume : list of array, contains all the controlpoint tuples for start, c1,c2, end

    """

    if not isinstance(cPoints, list):
        sys.exit("Input for interpolation must be list of controllpoint class")
    
    start= np.zeros([len(cPoints),2])
    c1 = np.zeros([len(cPoints),2])
    c2 = np.zeros([len(cPoints),2])
    end = np.zeros([len(cPoints),2])
    
    #shift the interpolation layers away from first layer
    Inter = zCoor[-1]- zCoor[0]
    Inter = Inter/numPoints
    zValues = np.linspace(zCoor[0] +Inter, zCoor[-1]- Inter, numPoints)
    
    for d in range(len(zValues)):
        if d ==0 :
            continue
        if zValues[d]  != zValues[d-1]:
            zValues[d] = zValues[d]+2
    #print(zValues)
        
    #print(zValues)
    #zValues = 
    #print(zValues) 
    
    for i in range(len(cPoints)):
        start[i,0] = cPoints[i].start[0]
        start[i,1] = cPoints[i].start[1]
        c1[i,0] = cPoints[i].c1[0]
        c1[i,1] = cPoints[i].c1[1]
        c2[i,0] = cPoints[i].c2[0]
        c2[i,1] = cPoints[i].c2[1]
        end[i,0] = cPoints[i].end[0]
        end[i,1] = cPoints[i].end[1]
    #print(start[:,0])


    # sX = np.linspace(start[0,0], start[-1,0], numPoints)
    # c1X = np.linspace(c1[0,0], c1[-1,0], numPoints)
    # c2X = np.linspace(c2[0,0], c2[-1,0], numPoints)
    # endX =  np.linspace(end[0,0], end[-1,0], numPoints)

    #interpolation for X
    startX = interpolate.interp1d(zCoor, start[:,0])
    c1X = interpolate.interp1d(zCoor, c1[:,0])    
    c2X = interpolate.interp1d(zCoor, c2[:,0])
    endX = interpolate.interp1d(zCoor, end[:,0])
    
        #interpolation for Y
    startY = interpolate.interp1d(zCoor, start[:,1])
    c1Y = interpolate.interp1d(zCoor, c1[:,1])    
    c2Y = interpolate.interp1d(zCoor, c2[:,1])
    endY = interpolate.interp1d(zCoor, end[:,1])
    
    #print(c1X(zValues))
    VolumeX = np.array([startX(zValues),c1X(zValues),
                        c2X(zValues), endX(zValues) ])
    VolumeY = np.array([startY(zValues),c1Y(zValues),
                        c2Y(zValues), endY(zValues) ])
    # print("-----------------------------------")
    # print(VolumeY[0])
    # print(start[:,1])
    
    Segment =[]
    for r in range(numPoints):
        CP = np.array([])
        CP = np.array([VolumeX[:,r],VolumeY[:,r]])
        BEZ = makeBezier(CP)
        Segment.append(BEZ)
    #Volume = list([startI,c1I, c2I, endI])
    #Volume = startX(zValues)
    #print(Volume)
    #Volume = list([np.array([sX,startI(sX)]),np.array([c1X,c1I(c1X)]), 
                   #np.array([c2X,c2I(c2X)]), np.array([endX,endI(endX)])])

    return Segment 
    
def interInter(cPoints, interZ, numPoints):
    
    #Inter = interZ[1]- interZ[0]
    #Inter = Inter/numPoints
    zValues = np.linspace(interZ[0], interZ[-1], numPoints)
    
    start= np.zeros([2,2])
    c1 = np.zeros([2,2])
    c2 = np.zeros([2,2])
    end = np.zeros([2,2])
    
    for i in range(len(cPoints)):
        start[i,0] = cPoints[i].start[0]
        start[i,1] = cPoints[i].start[1]
        c1[i,0] = cPoints[i].c1[0]
        c1[i,1] = cPoints[i].c1[1]
        c2[i,0] = cPoints[i].c2[0]
        c2[i,1] = cPoints[i].c2[1]
        end[i,0] = cPoints[i].end[0]
        end[i,1] = cPoints[i].end[1]
        
    startX = interpolate.interp1d(interZ, start[:,0])
    c1X = interpolate.interp1d(interZ, c1[:,0])    
    c2X = interpolate.interp1d(interZ, c2[:,0])
    endX = interpolate.interp1d(interZ, end[:,0])
    
        #interpolation for Y
    startY = interpolate.interp1d(interZ, start[:,1])
    c1Y = interpolate.interp1d(interZ, c1[:,1])    
    c2Y = interpolate.interp1d(interZ, c2[:,1])
    endY = interpolate.interp1d(interZ, end[:,1])
    
    VolumeX = np.array([startX(zValues),c1X(zValues),
                        c2X(zValues), endX(zValues) ])
    VolumeY = np.array([startY(zValues),c1Y(zValues),
                        c2Y(zValues), endY(zValues) ])
    
    Segment =[]
    for r in range(numPoints):
        CP = np.array([])
        CP = np.array([VolumeX[:,r],VolumeY[:,r]])
        BEZ = makeBezier(CP)
        Segment.append(BEZ)
    
    
    return Segment, zValues
    

    
def interWrap(inFile, numInterLayers):
    zCoor = getZvalues(inFile)
    
    idx, vals = compEmpty(zCoor, numInterLayers)

    cPoints = getCpoints(inFile)
    interLayers= list()
    iS = list()
    zC = np.array([])
    #terrible hardcoding bug
    bezier = list()
    for p in range(len(cPoints[0])):
        bezier = list()
        for r in range(len(cPoints)):
            bezier.append(cPoints[r][p])
            if r== len(cPoints)-1:
                continue
            
            i1 = cPoints[r][p]
            i2 = cPoints[r+1][p]#
            intervall = list((i1,i2))
            #print(intervall)
            interZ = np.array([zCoor[r],zCoor[r+1]])
            #print(intervallZ)
            ss, z = interInter(intervall, interZ, 50)
            iS.append(ss)
            
            
                
        zC = np.append([zC],z)
            #bezier = list([cPoints[0][p],cPoints[1][p],cPoints[2][p],cPoints[3][p]])
        seg = interp(bezier, zCoor, numInterLayers)
        #print(seg)
        interLayers.append(seg)
    print(iS)
    reshape = list()
    for r in range(numInterLayers):
        layer = list()
        for p in range(len(cPoints[0])):    
    
            x = interLayers[p][r]
            layer.append(x)
        reshape.append(layer)
    #learn enumerate
    count = 0
    # print(idx)
    # print(len(idx))
    # print(len(reshape))
    for i in range(len(idx)):
        if idx[i] == 1:
            #print("c")
            cPoints.insert(i, reshape[count])
            count +=1
            
        
    return cPoints, vals               

        



def deCastel(cPoints, t = 0.5):
    """
    De Casteljau algorthm wich determines the absolute coordinates
    on a bezier curve for a given segment

    Parameters
    ----------
    cPoints : 4 controlpoints start, c1,c2,end.
    t :  where the coordinate is computed.

    Returns
    -------
    B : coordinates

    """
    Control = np.asarray([cPoints.start, cPoints.c1, cPoints.c2,cPoints.end])
    N,d = Control.shape
    order = np.arange(N)
    coeff = [math.factorial(N-1)
             // (math.factorial(i)* math.factorial(N-1-i))
             for i in range(N)]

    Px = (Control.T *coeff).T

    B = (np.power.outer(1-t, order[::-1])
         *np.power.outer(t,order)) @Px

    return B
#todo : correct z axis position in pointcloud
#numLayers is number of interpolation steps
#bug with only 2 inter
#
def getCarthesian(inFile, numInterLayers, prec):
    
    cPoints, vals = interWrap(inFile, numInterLayers)
    #print(vals)
    #print(cPoints)
    t = np.linspace(0,1, prec)
    carthCoors = list()
    LayerCoors = np.array([])
    print(len(vals))
    for r in range(len(vals)):
        
        for p in range(len(cPoints[0])):
            seg = cPoints[r][p]
            for s in range(prec):
                B = deCastel(seg,t[s])
                
                B = np.append(B,[vals[r]]) 
                
                LayerCoors = np.append([B],LayerCoors)
    newshape = int(len(LayerCoors)/3)
    LayerCoors = np.reshape(LayerCoors,(newshape,3))
    return LayerCoors




def triNormals(tri,lc):
    normals = np.zeros((len(tri),3))
    for i in range(len(tri)):
        p1 = lc[tri[i,0],:]
        p2 = lc[tri[i,1],:]
        p3 = lc[tri[i,2],:]
        norm = np.cross(p2-p1, p3-p1)
        normals[i,:] = norm
    return normals

def triSurfOpen(inFile, nInter, nPrec):
    
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

def triSurfClose(inFile, nInter, nPrec):
    
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
            for i in range(numP):

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
    
            for i in range(numP):

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


    


inFile = "input/slab.svg"
paths, attributes = svg2paths(inFile)
#sortLayers(inFile)

# f = open(inFile)#
# text = f.readlines()
Layers, numLayers = getLayers(inFile)
#checkPath(attributes, Layers)
interLayers= list()
cPoints = getCpoints(inFile)
zCoor = getZvalues(inFile)
# for p in range(len(cPoints[0])):
#         bezier = list([cPoints[0][p],cPoints[1][p],cPoints[2][p],cPoints[3][p]])
#         seg = interp(bezier, zCoor, 5)
#         interLayers.append(seg)
idx, val = compEmpty(zCoor, 5)
lc = getCarthesian(inFile,5,30)

#triangles, face = triSurfOpen(inFile,5,30)

    
#bezier = list([cPoints[0][0],cPoints[1][0],cPoints[2][0],cPoints[3][0]])

#cPointsnew, vals = interWrap(inFile, numInterLayers=6)


def getBounds():
    return





    
    # mesh[p] = lc[triangles[p,0]]
    # mesh[p+1] = lc[triangles[p,1]]
    # mesh[p+2] = lc[triangles[p,2]]
    # mesh[p+3] = normals[p,:]
     

#from here it gets messy
from scipy.spatial import ConvexHull, convex_hull_plot_2d,Delaunay 
# dlnSC = Delaunay(lc)
# hull = ConvexHull(lc)
import numpy
import struct
import stl
from stl import mesh

# cube = mesh.Mesh(np.zeros(face.shape[0], dtype=mesh.Mesh.dtype))
# for i, f in enumerate(face):
#     for j in range(3):
#         cube.vectors[i][j] = lc[f[j],:]

#Write the mesh to file "cube.stl"
#cube.save('debug.stl')


# def sw(ver):
#     class stl(object):
        
#         def __init__(self, n, v1, v2, v3):
#             self.n = n
#             self.v1 = v1
#             self.v2 = v2
#             self.v3 = v3
#             return
        
#     n = ver[:,0]
#     v1 = ver[:,1]
#     v2 = ver[:,2]
#     v3 = ver[:,3]
    
#     mesh = stl(n,v1,v2,v3)
#     with open("slab.stl", mode="wb")as fh:
#         fh.write(69)
#         fh.write(struct.pack('@i'))

        

class Stl(object):
    dtype = numpy.dtype([
        ('normals', numpy.float32, (3, )),
        ('v0', numpy.float32, (3, )),
        ('v1', numpy.float32, (3, )),
        ('v2', numpy.float32, (3, )),
        ('attr', 'u2', (1, )),
    ])

    def __init__(self, header, data):
        self.header = header
        self.data = data

    @classmethod
    def from_file(cls, filename, mode='rb'):
        with open(filename, mode) as fh:
            header = fh.read(80)
            size, = struct.unpack('@i', fh.read(4))
            data = numpy.fromfile(fh, dtype=cls.dtype, count=size)
            return Stl(header, data)

    def to_file(self, filename, mode='wb'):
        with open(filename, mode) as fh:
            fh.write(self.header)
            fh.write(struct.pack('@i', self.data.size))
            self.data.tofile(fh)

# head = "head"
# head = bytes(head, 'utf-8')
# file = Stl(head,mesh)
# file.to_file('mesh.stl')

    
def write(ver):
    with open("slab.stl", mode="wb")as fh:
        writer = wr.Binary_STL_Writer(fh)
        writer.add_faces(ver)
        writer.close()
    return
# sim = dlnSC.simplices
# vertex = dlnSC.vertices
# #vertex = lc[ind]
# from stl import mesh
# volume = mesh.Mesh(np.zeros(sim.shape[0], dtype=mesh.Mesh.dtype))
# for i, f in enumerate(sim):
#     for j in range(3):
#         volume.vectors[i][j] = vertex[f[j],:]
# volume.save('slab3.stl')

# dln = delny3D(lc)

# head = "head"
# head = bytes(head, 'utf-8')
# file = Stl(head,dlnSC.points)
# file.to_file('rick.stl')

# test = Stl.from_file('rick.stl')
# from vtk.util.numpy_support import vtk_to_numpy
# numpy_coordinates = vtk_to_numpy(dln)

# reader = vtk.vtkXMLUnstructuredGridReader()
# reader.SetFileName( "slabvtk.vtu" )
# reader.Update()


# Point_cordinates = reader.GetOutput().GetPoints().GetData()

# #write(vertex)
#     #file = stlWrite(dln.vertices)
# # fig = plt.figure()
# # ax = plt.axes(projection='3d')
# # ax.scatter3D(lc[:,0], lc[:,1], lc[:,2])
# # plt.show()

# # #dln.write()
# writer = vtk.vtkXMLUnstructuredGridWriter()
# #writer.SetFileType("stl")
# writer.SetFileName('slabvtk.vtu')
# writer.SetInputData(dln)
# writer.Write()

#pcd = o3d.geometry.PointCloud()
#pcd.points = o3d.utility.Vector3dVector(lc[:,:3])
#pcd.colors = o3d.utility.Vector3dVector(lc[:,3:6]/255)
#pcd.normals = o3d.utility.Vector3dVector(lc[:,6:9])
#o3d.visualization.draw_geometries([pcd])
#hull = pcd.compute_convex_hull()
#----------pivoting method
# distances = pcd.compute_nearest_neighbor_distance()
# avg_dist = np.mean(distances)
# radius = 3 * avg_dist

# #bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd,o3d.utility.DoubleVector([radius, radius * 2]))

#-----------Poisson rec------------

#poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8, width=0, scale=1.1, linear_fit=False)[0]
inFile = "input/slab.svg"
def plotCloud3D(inFile):
    
    import open3d as o3d
    Layers, numLayers = getLayers(inFile)
    lc = getCarthesian(inFile,50,25)
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

# for i in range(len(cPoints)):

#     plotBezier(cPoints[i])
def groupSegments(cPoints):
    """
    turns the layer list into a list of equal segments
    """
    segments= list()
    for i in range(len(cPoints[0])):
        segs = list()
        for p in range(len(cPoints)):
            c = cPoints[p][i]
            segs.append(c)
        segments.append(segs)
    return segments



        
        



def Dim(cPoints, nLayers = 15):
    segments = groupSegments(cPoints)
    adLayers = list()
    for i in range(len(segments)):
        segs = list()
        iPoints = interp(segments[i], nLayers)
        print(len(iPoints))
        newbez = makeBezier(iPoints)
        

        adLayers.append(newbez)
    return  adLayers

#l = Dim(cPoints)
#segments = groupSegments(cPoints)
    

# def getCoors(cPoints, nPoints=30):
#     for i in range(len(cPoints)):

    
    
    




        
        

    