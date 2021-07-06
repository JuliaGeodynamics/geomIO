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
    


    
    
def controlPoints(line):
    cPoints = list()
    class Bezier(object):
    
        def __init__(self, start,c1,c2, end):
            self.start = start
            self.c1 = c1
            self.c2 = c2
            self.end = end
            
            return
    
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



import matplotlib.pyplot as plt
from matplotlib.path import Path as Pt
import matplotlib.patches as patches

# verts = [
#    (0., 0.),   # P0
#    (0.2, 1.),  # P1
#    (1., 0.8),  # P2
#    (0.8, 0.),  # P3
# ]

# codes = [
#     Pt.MOVETO,
#     Pt.CURVE4,
#     Pt.CURVE4,
#     Pt.CURVE4,
# ]

# path = Pt(verts, codes)

# fig, ax = plt.subplots()
# patch = patches.PathPatch(path, facecolor='none', lw=2)
# ax.add_patch(patch)

# xs, ys = zip(*verts)
# ax.plot(xs, ys, 'x--', lw=2, color='black', ms=10)

# ax.text(-0.05, -0.05, 'P0')
# ax.text(0.15, 1.05, 'P1')
# ax.text(1.05, 0.85, 'P2')
# ax.text(0.85, -0.05, 'P3')

# ax.set_xlim(-0.1, 1.1)
# ax.set_ylim(-0.1, 1.1)
# plt.show()



def plotBezier(cPoints):
    
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
        

        




















# #print(paths)
# p = line[0]
# pol = p.poly()
# print(pol(1))
# pe = line[1]
#
# import matplotlib.patches as mpatches

# #plt.figure()
# for i in range(len(paths)):
    
#     p = []
#     line = paths[i]
#     cor = line2coor(line)
#     # seg =line._segments
    
#     array = np.array(cor)#*tol
#     #print(array)
    
     
#     #plt.plot(array[:, 0], array[:, 1], color='b')
    
#     #plt.scatter(array[:, 0], array[:, 1], s=5)
    
#     #plt.show()
#arr = np.array([pe])

# 

# # bezier = line._segments
# # c = bezier[0].control1
# # cc = real(c)
# c = controlPoints(line)
# ccc = c[0]
# Path = mpath.Path

# fig, ax = plt.subplots()
# pp1 = mpatches.PathPatch(
#     Path([ccc[0].start, ccc[0].c1, ccc[0].c2, ccc[0].end],
#          [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.CLOSEPOLY]),
#     fc="none", transform=ax.transData)
# plt.show()


#drawing = svg2rlg("input/vector.svg")




