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
    Paths = dict()        
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
                    Paths[pathSTR] = layer
                else:                    
                    compNum +=1
                    #print (compNum)
                    pathSTR = text[p]
                    pathSTR = pathSTR.split("\"")
                    pathSTR = pathSTR[1]
                    Paths[pathSTR] = layer

        #Layers.append(Paths)
    return Paths, numLayers

def plot_line(inp, tol):
    
    line = np.array(inp)#*tol
     
    #plt.plot(line[:, 0], line[:, 1], color='b')
    
    plt.scatter(line[:, 0], line[:, 1], s=5)


def path2numpy(inFile):
    return

def layerSeq(attributes, Layers):
    
    #Check for dicts
    if not isinstance(attributes, list):
        sys.exit("Attributes must be dictonary")
    if not isinstance(Layers, dict):
        sys.exit("Layers must be dictonary")
        
        
    for i in Layers.keys():
        print(i)
        #if i in attributes.values():
        #    print("dasads")
    return

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
    layerInfo = getLayers(inFile)
    Layers = list()
    for i in layerInfo.keys():
        ly = list()
        for p in range(len(attributes)):
            if attributes[p]['id'] in i: 
                ly.append(paths[p])
        Layers.append(ly)
        
    return Layers
            
                       
    
    
        
        
        
        
        
        
inFile = "input/vector1.svg"
paths, attributes = svg2paths(inFile)
#asd = sortLayers(inFile)

f = open(inFile)
text = f.readlines()
Layers, numLayers = getLayers(inFile)
layerSeq(attributes, Layers)





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

# bezier = line._segments
# c = bezier[0].control1
# cc = real(c)
c = controlPoints(line)
ccc = c[0]
# Path = mpath.Path

# fig, ax = plt.subplots()
# pp1 = mpatches.PathPatch(
#     Path([ccc[0].start, ccc[0].c1, ccc[0].c2, ccc[0].end],
#          [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.CLOSEPOLY]),
#     fc="none", transform=ax.transData)
# plt.show()


#drawing = svg2rlg("input/vector.svg")
#x = drawing.getBounds()

# x,y,z  = pt.svg2paths2("input/dom.svg")
# x_1 = x[0]
# b = x_1.point(0)

# Bez= p[56]



#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

# def poly_plot(pdict, tol, isOpen = False) :

#     for poly in pdict.values() :

#         for i in range(len(poly)):
        
#             coord = np.array(poly[i])*tol
            
#             plt.scatter(coord[:, 0], coord[:, 1], s=5)
        
#             if not isOpen:
            
#                 coord = np.vstack((coord, coord[0, :]))

#             plt.plot(coord[:, 0], coord[:, 1], color='k')

#     plt.axis('equal')
       
# #------------------------------------------------------------------------------

# def poly_clip(clip, subj, solution) :

#     pc = pyclipper.Pyclipper()
    
#     pc.AddPath(clip, pyclipper.PT_CLIP,    True)
#     pc.AddPath(subj, pyclipper.PT_SUBJECT, True)

#     res = pc.Execute2(pyclipper.CT_DIFFERENCE, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD)

#     sol   = [i.Contour for i in res.Childs]
#     holes = [i.IsHole  for i in res.Childs]

#     for i in range(len(sol)) :
        
#         if not holes[i] :

#             solution.append(sol[i])

# #------------------------------------------------------------------------------

# def poly_get_box(inp) :
    
#     poly = np.array(inp)
    
#     box = np.vstack((poly.min(axis=0), poly.max(axis=0)))
    
#     return box

# #------------------------------------------------------------------------------

# def points_in_box(points, box) :
    
#     pt   = np.array(points)
#     xmin = np.array(pt[:, 0] > box[0, 0])
#     xmax = np.array(pt[:, 0] < box[1, 0])
#     ymin = np.array(pt[:, 1] > box[0, 1])
#     ymax = np.array(pt[:, 1] < box[1, 1])

#     I = xmin & xmax & ymin & ymax
    
#     return pt[I, :]

# #------------------------------------------------------------------------------

# def add_path(coord, ID, pdict) :

#     if ID in pdict :
        
#         pdict[ID].append(coord)
        
#     else :
        
#         pdict[ID] = [coord]

# #------------------------------------------------------------------------------

# def read_path(path, att, solids, holes, lines, tol) :
    
#     # check attribute
#     A = att[0]
          
#     if A not in ['S', 'H', 'L'] :
        
#         sys.exit("Attributes should define solids (S), holes (H) or lines (L) " + att)

#     if A in ['S', 'H'] and not path.isclosed() :
        
#         sys.exit("Solid and hole paths must be closed " + att)

#     if A in ['L'] and path.isclosed() :
        
#         sys.exit("Line paths must be opened " + att)

#     coord = []
    
#     # read coordinates
#     for i in range(len(path)) :
        
#         if not isinstance(path[i], Line) :
            
#             sys.exit("All paths should consist of lines segments only")
        
#         line =  path[i]
#         p    =  line[0]
#         pe   =  line[1]
#         x    =  int(real(p)/tol)
#         y    = -int(imag(p)/tol)

#         coord.append([x, y])
        
#     if A in ['L'] :
        
#         x    =  int(real(pe)/tol)
#         y    = -int(imag(pe)/tol)
        
#         coord.append([x, y])

    
#     # get attribute identifier (key)
#     ID = int(att[1:])
    
#     # store path
#     if A in ['S'] :
#         add_path(coord, ID, solids)
    
#     if A in ['H'] :
#         add_path(coord, ID, holes)
    
#     if A in ['L'] :
#         add_path(coord, ID, lines)
      
#------------------------------------------------------------------------------

#plt.close('all')

#tol = 1e-6

# Rread SVG into a list of path objects and list of dictionaries of attributes 
# paths, attributes = svg2paths('test.svg')
# paths, attributes = svg2paths('input/cavern.svg')

# solids = {}
# holes  = {}
# lines  = {}

# for i in range(len(paths)) :
    
#     if 'A' in attributes[i] :
   
#         att = attributes[i]['A']
        
#         read_path(paths[i], att, solids, holes, lines, tol)

#     else :
        
#         sys.exit("All paths should define attributes")

# # scoord = [coord[i] for i in sorted(coord.keys())]

# plt.figure()

# poly_plot(solids, tol)
# poly_plot(holes,  tol)
# poly_plot(lines,  tol, isOpen = True)

#------------------------------------------------------------------------------

# clip = solids[1][0]
# subj = [lines[1][0], lines[2][0], lines[3][0], lines[4][0]]

# pc = pyclipper.Pyclipper()

# pc.AddPath(clip, pyclipper.PT_CLIP, True)
# pc.AddPaths(subj, pyclipper.PT_SUBJECT, False)

# solution = pc.Execute2(pyclipper.CT_INTERSECTION, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD)

# plt.figure()
# poly_plot(solids, tol)

# sol = [i.Contour for i in solution.Childs]

# plot_line(sol[0], tol)
# plot_line(sol[1], tol)
# plot_line(sol[2], tol)
# plot_line(sol[3], tol)

#------------------------------------------------------------------------------

# scoord = scoord[0:5]

#http://www.angusj.com


# sol = []

# for i in range(0, len(scoord)) :
   
#     tmp = []
    
#     for j in range(len(sol)) :

#         poly_clip(scoord[i], sol[j], tmp)
    
#     sol = []
    
#     if len(tmp) :
#         sol.extend(tmp) 
    
#     sol.append(scoord[i])

# # for i in range(len(sol)) :
# #     plt.figure()
# #     pplot([sol[i]])

# plt.figure()
# poly_plot(sol)


# box = poly_get_box(sol[0])

# inbox = points_in_box(sol[1], box)


# si = [scoord[0]]

# # for i in range(1, len(scoord)) :
# for i in range(1,2) :

#     sj = []
    
#     for j in range(len(si)) :

#         pclip(scoord[i], si[j], sj)
   


# plt.figure()
# pplot(solution)


# subj = scoord[0:2]
# clip = scoord[4]

# plt.figure()
# pplot(subj)


# plt.figure()
# pplot(clip)

# # subj = coord['S1']
# # clip = coord['S4']

# pc = pyclipper.Pyclipper()
# pc.AddPath(clip, pyclipper.PT_CLIP, True)
# pc.AddPaths(subj, pyclipper.PT_SUBJECT, True)

# solution = pc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD)


# # hole   = [i.IsHole for i in solution.Childs]

# # parents = [i.Parent for i in solution.Childs]

# # solution = [i.Contour for i in solution.Childs]


# plt.figure()

# pplot(solution)







# coord = []

# for path in paths:
#     for line in path:
#         for point in line:
#             x = real(point)
#             y = imag(point)
#             # print(x, y)
#             coord.append((x,y))
            
# coord = np.array(coord)


# # # plot grid
# # plt.figure(1)
# # plt.plot(coord[:, 0], coord[:, 1])
# # plt.scatter(coord[:, 0], coord[:, 1], s=5)
# # plt.title('Geometry')
# # plt.xlabel('x')
# # plt.ylabel('y')
# # plt.axis('equal')


# # plot grid
# plt.figure(1)

# for path in paths:
    
# #    if not path.isclosed() :
# #        continue
    
#     # plt.figure()

#     for line in path:
        
#         p1 = line[0]
#         p2 = line[1]
        
#         p1x =  real(p1)
#         p1y = -imag(p1)
#         p2x =  real(p2)
#         p2y = -imag(p2)
        
#         plt.plot((p1x, p2x), (p1y, p2y), color='k')
#         # plt.scatter((p1x, p2x), (p1y, p2y), color='g', s=5)

# plt.title('Geometry')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.axis('equal')
# plt.show() 
# plt.draw()



#a = np.array([[1,2,3],[4,5,6],[0,0,1]])


#coord = coord[coord[:,1].argsort()] # First sort doesn't need to be stable.
#coord = coord[coord[:,0].argsort(kind='mergesort')]
#coord = coord[coord[:,0].argsort(kind='mergesort')]

#B = np.lexsort(np.transpose(coord))
#B = coord[ind] 



# tol = 1e-3

# ind     = np.lexsort((coord[:, 1], coord[:, 0]))
# uind    = np.zeros(ind.shape, dtype=int)
# idx     = ind[0]
# uind[0] = idx
# val     = coord[idx]
# nval    = np.zeros(val.shape)
# num     = 1


# for i in range(1,len(ind)) :
    
#     nidx    = ind[i]
#     nval[:] = coord[nidx]
    
#     if np.linalg.norm(nval - val) > tol :
        
#         val[:] = nval[:]
#         idx    = nidx
#         num   += 1
        
#     uind[i] = idx

        
        
        
        
        
    
# print(num)
    


