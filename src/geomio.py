# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# import np as np
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
            if compNum != initNum and "$" not in str(text[index[i]:index[i+1]])  :
                print(initNum,compNum)
                #sys.exit("Number of Paths must be equal per Layer. Invalid number in:" + str(layer))
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

def scaling():
    return
    
    
    
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

def ignore(inFile):
    
    
    Layers, numLayers = getLayers(inFile)
    ignore = np.zeros(len(Layers))
    for q,p in enumerate(Layers.values()):
        if "$" in p:
            ignore[q] = 1
        elif "§" in p:
            ignore[q] =1
    return ignore

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
    
    ###-----------hier rein basteln---------
    cPoints = list()
    ig  = ignore(inFile)
    pathsR = list()
    for r in range(len(ig)):
        if ig[r]  == 0:
            pathsR.append(paths[r])
    for i in range(len(pathsR)):
        line = pathsR[i]
        points = controlPoints(line)
        cPoints.append(points)
    
    return cPoints
    
def makeBezier(cPoints):
    """
    This functions converts the svgpathtools bezier to a simple bezier class
    with start, end controlpoint 1 and 2
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
    """
    this function reads the Z-coordinate values passed with the individual 
    layernames

    Parameters
    ----------
    inFile : name of the Inputfile

    Returns
    -------
    zCoor : list of Z coordinates

    """
    Layers, numLayers = getLayers(inFile)
    names = list()
    for z in Layers.values():
        names.append(z)
    zCoors = np.array([])
    zCoor = np.array([], dtype = float)    
    for i in range(len(names)):
        dig = []
        if "$"  in names[i]:
            continue
        if "§" in names[i]:
            continue
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
    """
    This function generates the Z coordinates where the Layers are interpolated
    numLayers refers to the number of Layers to interpolate(numInterLayers)
    
    """
    #print(zCoor)
    Inter = zCoor[-1]- zCoor[0]
    #print(Inter)
    Inter = Inter/numLayers
    zValues = np.linspace(zCoor[0] +Inter, zCoor[-1]- Inter, numLayers, False)
    #print(zValues)
    return zValues

#es ist halt wirklich die interpolation die von anfang bis end inter aber in der mitte ignoriert

def compEmpty(zCoor, numLayers):
    """
    this function returns a complete array of all z Coordinates
    as well as an index array, where 1= given layer and 0 = interpolated layer
    """
    zValues = interZlayers(zCoor, numLayers)
    zCoor = np.append([zValues], zCoor)
    vals = np.sort(zCoor)
    vals = np.flip(vals)
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
def interp(cPoints, zCoor, numPoints = 5, meth = 'linear'):
    """
    This is the 1D interpolation in y direction

    Parameters
    ----------
    cPoints : List of controlpoints for all zusammengehörigen 
    bezierabschnitten
    numPoints : Number of interpolationpoints
         The default is 5.
    meth = Crystal meth. Jk interpolation method, currently only linear is allowed

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
    # Inter = zCoor[-1]- zCoor[0]
    # Inter = Inter/numPoints
    # zValues = np.linspace(zCoor[0] +Inter, zCoor[-1]- Inter, numPoints)
    zValues = interZlayers(zCoor,numPoints)

            
    #print(zValues)
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
    startX = interpolate.interp1d(zCoor, start[:,0], meth)
    c1X = interpolate.interp1d(zCoor, c1[:,0],meth)    
    c2X = interpolate.interp1d(zCoor, c2[:,0],meth)
    endX = interpolate.interp1d(zCoor, end[:,0],meth)
    
        #interpolation for Y
    startY = interpolate.interp1d(zCoor, start[:,1],meth)
    c1Y = interpolate.interp1d(zCoor, c1[:,1],meth)    
    c2Y = interpolate.interp1d(zCoor, c2[:,1],meth)
    endY = interpolate.interp1d(zCoor, end[:,1],meth)
    
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
    

    

    
def interWrap(inFile, numInterLayers):
    """
    Function to call interpolation,
    returns a list of all new controlpoints properly sorted as well as all 
    z coordinates
    
    """
    zCoor = getZvalues(inFile)
    
    idx, vals = compEmpty(zCoor, numInterLayers)

    cPoints = getCpoints(inFile)
    interLayers= list()
#terrible hardcoding bug
    bezier = list()
    for p in range(len(cPoints[0])):
        bezier = list()
        for r in range(len(cPoints)):
            bezier.append(cPoints[r][p])

        seg = interp(bezier, zCoor, numInterLayers)

        interLayers.append(seg)

    reshape = list()
    for r in range(numInterLayers):
        layer = list()
        for p in range(len(cPoints[0])):    
    
            x = interLayers[p][r]
            layer.append(x)
        reshape.append(layer)
    #learn enumerate
    count = 0

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
    """
    This function computes the coordinates for all
    the points in the pointcloud, prec being the number 
    of points to compute per bezier segment 
    """
    cPoints, vals = interWrap(inFile, numInterLayers)

    #print(cPoints)
    t = np.linspace(0,1, prec)
    carthCoors = list()
    LayerCoors = np.array([])
    #print(len(vals))
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


 

# inFile = "input/over.svg"
# paths, attributes = svg2paths(inFile)
# #sortLayers(inFile)
# #seg = interDisc(inFile, 5)
# # f = open(inFile)#
# # text = f.readlines()
# Layers, numLayers = getLayers(inFile)
# #checkPath(attributes, Layers)
# #interLayers= list()
# cPoints = getCpoints(inFile)
#zCoor = getZvalues(inFile)
# for p in range(len(cPoints[0])):
#         bezier = list([cPoints[0][p],cPoints[1][p],cPoints[2][p],cPoints[3][p]])
#         seg = interp(bezier, zCoor, 5)
#         interLayers.append(seg)
# idx, val = compEmpty(zCoor, 5)
# lc = getCarthesian(inFile,5,30)
# a,b = callIntervall(inFile,5)


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

def triSurfOpen(inFile, nInter, nPrec):
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
    np.save("surf",mesh1)

    
    
    return mesh1, triangles

from scipy.spatial import ConvexHull, convex_hull_plot_2d,Delaunay 

#def fakeNormals(triangles, last = True):
    

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
    # if last:
        
    #     normals = np.zeros((len(triangles),3))
    #     normals[:,2] = 1
    # else:
    #     normals = np.zeros((len(triangles),3))
    #     normals[:,2] = -1
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

def triSurfClose(inFile, nInter, nPrec):
    """
    surface triangulation for closed volumes
    """
    Layers, numLayers = getLayers(inFile)
    lc = getCarthesian(inFile,nInter,nPrec)

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

    

    #=======section for inside: working    

    for p in range(nQuads):

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
    mesh1 = np.concatenate((mesh1,cover))
    triangles = np.concatenate((triangles,covTri))
    #for d in range(2):
    # mesh1 = np.concatenate((mesh1,mesh2))
    # mesh1 = np.concatenate((mesh1,mesh3))

    # triangles = np.concatenate((triangles,triCov1))
    # faceLast = triCov2+LLI
    # triangles = np.concatenate((triangles,faceLast))        
    
    #return cover, covTri
    return mesh1, triangles

# inFile = 'input/slab.svg'

# mesh1, tri = triSurfClose(inFile, 5,20)





def getBounds(inFile, numInter, nPrec):
    """
    printing the min and ma boundaries of the mesh/volume
    to adapt grid of thermomechanical coed(eg LaMEM)
    """
    lc = getCarthesian(inFile, numInter,nPrec)
    xBound = np.array([np.amin(lc[:,0]),np.amax(lc[:,0])])
    yBound = np.array([np.amin(lc[:,1]),np.amax(lc[:,1])])
    zBound = np.array([np.amin(lc[:,2]),np.amax(lc[:,2])])
    print("---------------------------")
    print("x boundaries are:"+ str(xBound))
    print("y boundaries are:"+ str(yBound))
    print("z boundaries are:"+ str(zBound))
    print("---------------------------")
    return


# dlnSC = Delaunay(lc)
# hull = ConvexHull(lc)



def wSTL(inFile, numInter, nPrec,name, volume = False, mode = "a"):
    """
    placeholder function for writing stl files
    uses the lib np stl
    """
    
    import struct
    import stl
    from stl import mesh
    if volume:
        triangles, face = triSurfClose(inFile,numInter,nPrec)
    else:
            
        triangles, face = triSurfOpen(inFile,numInter,nPrec)
    lc = getCarthesian(inFile, numInter,nPrec)
    os.chdir("output")
    cube = mesh.Mesh(np.zeros(face.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(face):
        for j in range(3):
            cube.vectors[i][j] = lc[f[j],:]
    
    #Write the mesh to file "cube.stl"
    if mode == "a":
        
        cube.save(str(name),mode=stl.Mode.ASCII )
    elif mode == "b":
        cube.save(str(name),mode=stl.Mode.BINARY )
    os.chdir("..")
    return


def wMulti(inFile, numInter, nPrec, names, volume):
    import struct
    import stl
    from stl import mesh
    
    
# lc = getCarthesian(inFile,5,25)
# getBounds(inFile,5,25)

#wSTL(inFile,25,50,'foldlay.stl')

        

class Stl(object):
    dtype = np.dtype([
        ('normals', np.float32, (3, )),
        ('v0', np.float32, (3, )),
        ('v1', np.float32, (3, )),
        ('v2', np.float32, (3, )),
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
            data = np.fromfile(fh, dtype=cls.dtype, count=size)
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
    header = 0
    with open("slab.stl", mode="wb")as fh:
#        writer = wr.Binary_STL_Writer(fh)
 #       writer.add_faces(ver)
  #      writer.close()
         fh.write(header)
         fh.write(struct.pack)
  
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
# from vtk.util.np_support import vtk_to_np
# np_coordinates = vtk_to_np(dln)

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
#inFile = "input/over.svg"
def plotCloud3D(inFile,nInterLayers, nPrec):
    """
    uses open3D to display the pointcloud
    """
    import open3d as o3d
    Layers, numLayers = getLayers(inFile)
    lc = getCarthesian(inFile,nInterLayers,nPrec)
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
    rather useless
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


def geomioFront(inFile, numInterLayers, nPrec, name, volume = False, mode = "a"):
    
    wSTL(inFile, numInterLayers, nPrec, name, volume, mode)
    
    return

#l = Dim(cPoints)
#segments = groupSegments(cPoints)
    

# def getCoors(cPoints, nPoints=30):
#     for i in range(len(cPoints)):

    
    
    




        
        

    