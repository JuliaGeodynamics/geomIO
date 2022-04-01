from svgpathtools import svg2paths
#from read_svg import getLayers
from read_svg import *
import numpy as np
from svgpathtools import svg2paths, real, imag, Line, svg2paths2, Document
from scipy import interpolate



import matplotlib.pyplot as plt
from matplotlib.path import Path as Pt
import matplotlib.patches as patches
import math
import sys,os

import scipy as sc

def splitPaths(data, key:str):
    paths = list()
    for i in range(len(data.Curves)):
        if data.CurveNames[i] == key:
            paths.append(data.Curves[i])
            
    return paths


def getZvalues(inFile):  # obsolete
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
        if "$" in names[i]:
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

def ignore(inFile): # Obsolete
    
    Layers, numLayers = getLayers(inFile)
    ignore = np.zeros(len(Layers))
    for q,p in enumerate(Layers.values()):
        if "$" in p:
            ignore[q] = 1
        elif "§" in p:
            ignore[q] =1
    return ignore


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


def getCpoints(data, path):
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

    
    
    # paths, attributes = svg2paths(inFile)
    # data = readSVG(inFile)
    # c = data.Curves
    # ###-----------hier rein basteln---------
    cPoints = list()
    # ig  = ignore(inFile)
    # pathsR = list()
    # for r in range(len(ig)):
    #     if ig[r]  == 0:
    #         pathsR.append(paths[r])
    # for i in range(len(pathsR)):
    #     line = pathsR[i]
    #     points = controlPoints(line)
    #     cPoints.append(points)
    
    for i in range(len(path)):
        line = path[i]
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

def interZlayers(zCoor,numLayers):
    """
    This function generates the Z coordinates where the Layers are interpolated
    numLayers refers to the number of Layers to interpolate(numInterLayers)
    
    """
    #print(zCoor)
    Inter = zCoor[-1]- zCoor[0]
    #print(Inter)
    Inter = Inter/numLayers
    #replace with interp1D
    x = zCoor.copy()
    y = np.arange(0, len(zCoor))
    f = interpolate.interp1d(y,x)
    ynew = np.linspace(0.1,y[-1], numLayers, endpoint= False)
    #zValues = np.linspace(zCoor[0] +Inter, zCoor[-1]- Inter, numLayers, False)
    zValues = f(ynew)
    #print(zValues)
    return zValues


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
    
    
def interWrap(data, path, numInterLayers):
    """
    Function to call interpolation,
    returns a list of all new controlpoints properly sorted as well as all 
    z coordinates
    
    """
    #zCoor = getZvalues(inFile)

    
    zCoor = np.array(list(set(data.zCoord)))
            
    
    idx, vals = compEmpty(zCoor, numInterLayers)

    cPoints = getCpoints(data, path)
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
def getPointCoords(data, path, numInterLayers, prec):
    """
    This function computes the coordinates for all
    the points in the pointcloud, prec being the number 
    of points to compute per bezier segment 
    """
    cPoints, vals = interWrap(data, path, numInterLayers)

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
    CoorR = np.zeros_like(LayerCoors)
    CoorR[:,0] = LayerCoors[:,0]
    CoorR[:,1] = LayerCoors[:,2]
    CoorR[:,2] = LayerCoors[:,1]
    return CoorR#LayerCoors
    #return LayerCoors
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

