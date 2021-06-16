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
        start = np.array([real(start),imag(start)])
        c1 = np.array([real(c1),imag(c1)])
        c2 = np.array([real(c2),imag(c2)])
        end = np.array([real(end),imag(end)])
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
            
        
        
        
    
    # pathList = list()
    # for i in layerInfo.keys():
    #     ly = list()
    #     for p in range(len(attributes)):
    #         if attributes[p]['id'] in i: 
    #             ly.append(paths[p])
    #     pathList.append(ly)
        
    return 



            
                       
    
    
        
        
        
        
        
        
inFile = "input/vector1.svg"
paths, attributes = svg2paths(inFile)
sortLayers(inFile)

f = open(inFile)
text = f.readlines()
Layers, numLayers = getLayers(inFile)
checkPath(attributes, Layers)

x = list(Layers.values())
pathList =sortLayers(inFile)




# #print(paths)
# p = line[0]
# pol = p.poly()
# print(pol(1))
# pe = line[1]

#arr = np.array([pe])

import matplotlib.path as mpath
import matplotlib.patches as mpatches

#plt.figure()
for i in range(len(paths)):
    
    p = []
    line = paths[i]
    cor = line2coor(line)
    # seg =line._segments
    
    array = np.array(cor)#*tol
    #print(array)
    
     
    #plt.plot(array[:, 0], array[:, 1], color='b')
    
    #plt.scatter(array[:, 0], array[:, 1], s=5)
    
    #plt.show()

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




