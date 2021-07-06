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
os.chdir("..")





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
    
    #Layers = list()
    Layers = dict()        
    #print(index)
    initNum = 0
    compNum = 0
    for i in range(len(index)-1):
        #print(i)#
        if i > 0:
            if compNum != initNum:

                sys.exit("Number of Paths must be equal per Layer. Invalid number in:" + str(layer))
        compNum = 0
        for p in range(index[i],index[i+1]):
            if "inkscape:label"  in text[p]:
                layerSTR = text[p]
                layerSTR = layerSTR.split("\"")
                #print(layerSTR)
                layer = layerSTR[1]
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
    
#write getcoords
#wrpng starting points, need to loop over starting points


def interp(cPoints, numPoints = 40):
    """
    This is the 1D interpolation in y direction

    Parameters
    ----------
    cPoints : List of controlpoints for all zusammengehörigen 
    bezierabschnitten
    numPoints : Number of interpolationpoints
         The default is 40.

    Returns
    -------
    Volume : list of array, contains all the controlpoint tuples for start, c1,c2, end

    """

    if not isinstance(cPoints, list):
        sys.exit("Input for interpolation must be list of controllpoint class")
    startTest = np.zeros([len(cPoints),2])
    startTest[:,0]= 1
    start= np.zeros([len(cPoints),2])
    c1 = np.zeros([len(cPoints),2])
    c2 = np.zeros([len(cPoints),2])
    end = np.zeros([len(cPoints),2])
    
    
    for i in range(len(cPoints)):
        start[i,0] = cPoints[i].start[0]
        start[i,1] = cPoints[i].start[1]
        c1[i,0] = cPoints[i].c1[0]
        c1[i,1] = cPoints[i].c1[1]
        c2[i,0] = cPoints[i].c2[0]
        c2[i,1] = cPoints[i].c2[1]
        end[i,0] = cPoints[i].end[0]
        end[i,1] = cPoints[i].end[1]
    # print(start)
    # print(start[:,0])
    # print(start[1,:])
    sX = np.linspace(start[0,0], start[-1,0], numPoints)
    c1X = np.linspace(c1[0,0], c1[-1,0], numPoints)
    c2X = np.linspace(c2[0,0], c2[-1,0], numPoints)
    endX =  np.linspace(end[0,0], end[-1,0], numPoints)
    #print(sX)    
    startI = interpolate.interp1d(start[:,0], start[:,1])
    c1I = interpolate.interp1d(c1[:,0], c1[:,1])    
    c2I = interpolate.interp1d(c2[:,0], c2[:,1])
    endI = interpolate.interp1d(end[:,0], end[:,1])    
    #print(c1)
    #print(c1X)
    #print(startI(sX))
    #print(c1I(c1X))
    #Volume = list([startI,c1I, c2I, endI])
    Volume = list([np.array([sX,startI(sX)]),np.array([c1X,c1I(c1X)]), 
                   np.array([c2X,c2I(c2X)]), np.array([endX,endI(endX)])])

    return Volume 
                       

        
inFile = "input/volume.svg"
paths, attributes = svg2paths(inFile)
sortLayers(inFile)

#f = open(inFile)
#text = f.readlines()
Layers, numLayers = getLayers(inFile)
checkPath(attributes, Layers)


cPoints = getCpoints(inFile)

    
bezier = list([cPoints[0][0],cPoints[1][1],cPoints[2][2],cPoints[3][3]])

vol = interp(bezier)

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


# curve = cPoints[0]
# plotCas(curve,20)



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


def makeBezier(cPoints):
    """
    
    """
    if len(cPoints) != 4:
        sys.exit("bezier Curves must be cubic")
    bezierList = list()
    for i in range(len(cPoints[0])):
        start = cPoints[0][i]
        c1 = cPoints[1][i]
        c2 = cPoints[2][i]
        end = cPoints[3][i]
        Points = Bezier(start,c1,c2,end)
        bezierList.append(Points)
    return bezierList
        
        



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

l = Dim(cPoints)
segments = groupSegments(cPoints)
    

# def getCoors(cPoints, nPoints=30):
#     for i in range(len(cPoints)):

    
    
    




        
        

    