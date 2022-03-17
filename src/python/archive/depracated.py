#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 16:55:51 2021

@author: lucas
"""

def interTRiangleFast(v1,v2,v3, zVal):
    
    """
    https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution    
    
    """
    avg = (v1+v2+v3)/3
    if zVal < avg :
        return True
    else:
        return False
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
            nInter = 0
            if intersec(grid[:,i],top,v1,v2,v3):
                nInter +=1
                if nInter == 0:
                    Phase[i]= 0
                elif nInter %2 !=0:
                    
                    Phase[i]=1
            cS +=4
            cE +=4
    return Phase

class Tri(object):
    
    def __init__(self, v1,v2,v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        
    @classmethod 
    def ray(cls, coor1, coor2):
        if intersec(coor1, coor2, v1,v2,v3):
            return True
        else:
            return False
        
        
#@jit(nopython=True) 
def mainTest(triangles, lc, grid):
    upBound = np.amax(grid[2,:])
    cells = len(grid[0,:])
    Phase = np.zeros([cells])
    
    
    for p in range(len(triangles)):
        #loop over all triangles
        print(p)
        #currentTri = np.array([lc[triangles[p,0]],lc[triangles[p,1]],lc[triangles[p,2]]])
        v1 = lc[triangles[p,0]]
        v2 = lc[triangles[p,1]]
        v3 = lc[triangles[p,2]]
        for r in range(len(Phase)):
            #loop over coordinates
            top = grid[:,r].copy()
            top[2] = upBound
            nInter = 0
            
            if intersec(grid[:,r],top,v1,v2,v3):
                nInter +=1
                if nInter == 0:
                    Phase[r]= 0
                elif nInter %2 !=0:
                    
                    Phase[r]=1
    return Phase

def fileTest(inFile, grid):
    
    
    return
# def fastRayFile(inFile:str, grid):
    
#     zVals, lowerBound, Layers, numPoints, points = sortGrid(grid)
#     triangles = mesh.Mesh.from_file(inFile)


#     numPoints = len(points)                               # number of cells per Layer
#     numberInter = np.zeros(numPoints)                     # number of intersections with triangles
#     interTri = np.array([])                               # index wich traingles are inntersected
#     InterPointZ = np.array([])                            # Z coordinate of 3D intersection Point     
    
#     for i in range(numPoints):
#         print("Computing 2D for column" + str(i))         # check all x  and y coordinates for intersection
#         for p in range(len(triangles)):                 
#                 v1 = triangles.v0[p,:]
#                 v2 = triangles.v1[p,:]
#                 v3 = triangles.v2[p,:]
#                 vertex = np.array([[v1[0], v1[1]],[v2[0], v2[1]],[v3[0], v3[1]]]) # transforming the triangle vertices into 2D
#                 if insideTriangle2D(points[i,:], vertex):                         # Use barycentric weights to check intersection
                    
#                     numberInter[i] +=1                      # count the triangles a point intersects
#                     interTri = np.append([p], interTri)     # save the index to triangle that intersects
#                     # Save the Z coordinate of the intersection Point on the triangle
#                     InterPointZ = np.append([upperTest(v1,v2,v3,np.array([points[i,0],points[i,1], lowerBound]))], InterPointZ)
#                 else:
#                     continue
#     print("intersection 2D completed")
#     triCounter = 0 # used to count indeces for InterTri and InterPointZ
#     Phase = np.array([]) # Phase array(comes flattend)
    

                            
#     for s in range(numPoints):  # test for all coordinates
               
#             print("computing for column" + str(s))
#             for p in range(Layers):
                
#                 if numberInter[s] ==0:                        # if no intersection whole column will be ignored
            
#                     Phase = np.append([0],Phase)              # Add zeros for every cell with no intersection
#                 else:
#                     intersections = 0                         # Count number of total intersections
#                     for t in range(int(numberInter[s])):
#                         # triIndex = int(interTri[triCounter])  # Call the triangle that intersects
#                         # v1 = lc[triangles[triIndex,0]]        # Get its coordinates
#                         # v2 = lc[triangles[triIndex,1]]        # this is completely redundant
#                         # v3 = lc[triangles[triIndex,2]]
#                         # v1 = v1[2]
#                         # v2 = v2[2]
#                         # v3 = v3[2]
#                         if zVals[p] <= InterPointZ[triCounter]: # get the Z elevation of the intersection point and see 
#                             intersections += 1                  # wether is underneath the surface or not
#                         else:
#                             continue
#                         # if interTRiangleFast(v1,v2,v3, zVals[p]) == True:
#                         #     intersections += 1
#                         # else:
#                         #     continue hier war btw der fehler weil hier kein counter lief
#                     if intersections % 2 == 0:                  # if even number of intersections point is above surface 
#                         Phase = np.append([0], Phase)
#                     else:
#                         Phase = np.append([1],Phase)            # else it is underneath
#             triCounter +=int(numberInter[s])                    # update indices
                    

                    
#     return Phase
    
    


#normalvektoren weg


def multiLayerRay(inFile: list, grid):
    
    triangles = mesh.Mesh.from_file(inFile)
    vertices = triangles.points
    x,y,z = grid[0], grid[1], grid[2]
    X,Y = x[:,:,0], y[:,:,0]
    Z = z[0,0,:]
    Phase = np.zeros_like(x) # Phase array(comes notflattend)

    nx, ny, nz = Phase.shape[0],Phase.shape[1],Phase.shape[2]

    numberInter = np.zeros_like(X)
    #interTri = np.array([])
    InterPointZ = np.array([])
    
    for j in range(ny):
        for i in range(nx):
            print("computing cell " + "x " + str(i) +", y " + str(j))
            for p in range(len(vertices)):
                             
                
                vertex = np.array([[vertices[p,0], vertices[p,1]],[vertices[p,3], vertices[p,4]],
                                   [vertices[p,6], vertices[p,7]]]) # transforming the triangle vertices into 2D
                if insideTriangle2D(np.array([X[i,j], Y[i,j]]), vertex):                         # Use barycentric weights to check intersection
                    
                    numberInter[i,j] +=1                      # count the triangles a point intersects
                    #interTri = np.append([p], interTri)     # save the index to triangle that intersects
                       # Save the Z coordinate of the intersection Point on the triangle
                       #test 
                    #SecPoint =upperTest(v1,v2,v3,np.array([X[i,j], Y[i,j], Z[0]]))
                    InterPointZ = np.append([upperTest(vertices[p,0:3],vertices[p,3:6],vertices[p,6:9],np.array([X[i,j], Y[i,j], Z[0]]))], InterPointZ)
                    #print(SecPoint)
                else:
                    continue
    InterPointZ = np.flip(InterPointZ)            
    print("intersection 2D completed")
     # used to count indeces for InterTri and InterPointZ
     #test with more intersection(flip current file by 90°)
    for k in range(nz):
        triIdx = 0
        for j in range(ny):
            
                for i in range(nx):
                    #print("computing cell " + "x " + str(i) +", y " + str(j))
                    #print(triIdx)
                    
                    
                    if numberInter[i,j] == 0:
                        Phase[i,j,k]= 0
                    else:
                        intersections = 0  
                        for t in range(int(numberInter[i,j])):
                            if Z[k] <= InterPointZ[int(triIdx)]: # get the Z elevation of the intersection point and see 
                                intersections += 1
                            triIdx += 1
                        if intersections % 2 == 0: 
                            Phase[i, j,k] = 0
                        else:
                            Phase[i,j,k] = 1
    
    
    return Phase


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
    """
    Uses pluecker Coordinate system to do determine line triangle intersection
    currently unused
    Documentation: https://members.loria.fr/SLazard/ARC-Visi3D/Pant-project/files/Line_Triangle.html
    """
    
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

def callIntervall(inFile, numInterLayers):
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
            #bezier.append(cPoints[r][p])
            if r== len(cPoints)-1:
                continue
            
            i1 = cPoints[r][p]
            i2 = cPoints[r+1][p]#
            intervall = list((i1,i2))
            #print(intervall)
            interZ = np.array([zCoor[r],zCoor[r+1]])
            #print(intervallZ)
            ss, z = interInter(intervall, interZ, 5)
            iS.append(ss)
            
            
            if p<0:
                continue
            zC = np.append([zC],z)
            #bezier = list([cPoints[0][p],cPoints[1][p],cPoints[2][p],cPoints[3][p]])
        #seg = interp(bezier, zCoor, numInterLayers)
        #print(seg)
        #interLayers.append(seg)
    #print(iS)
    reshape = list()
    for r in range(numInterLayers):
        layer = list()
        for p in range(len(cPoints[0])):    
    
            x = 0
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
            
        
    return zC, iS                         

def Inter(segment, numInter,inFile, x):
    zCoor = getZvalues(inFile)
    Layers, numLayers = getLayers(inFile)
    
    
    idx, zVal = compEmpty(zCoor, numInter)

    
    stab = zCoor[1]- zCoor[0]
    stab = stab/numLayers
    zValues = np.linspace(zCoor[0]+ stab, zCoor[-1], numInter, False)
    
    
    xvals = segment[:,0]
    yvals= segment[:,1]

    xVals = np.reshape(xvals,(numLayers, len(x)))
    yVals = np.reshape(yvals,(numLayers,len(x)))
    volume = np.array([])
    volX = np.array([])
    volY = np.array([])
    for i in range(len(x)):
        interX = interpolate.interp1d(zCoor, xVals[:,i])
        interY = interpolate.interp1d(zCoor, yVals[:,i])
        #coors = np.array([[interX(zValues)],[interY(zValues)]])
        #volume = np.append([coors], volume)
        xC = np.array([interX(zValues)])
        yC = np.array([interY(zValues)])
        volX = np.append([volX], xC)
        volY = np.append([volY], yC)
    volX = np.reshape(volX,(numInter, len(x)),'F')
    volY = np.reshape(volY,(numInter, len(x)),'F')

    
    #volume = np.reshape(volume,(len(x)*numInter,2))        
    return volX, volY


#remove hardcoding

def interDisc(inFile, numInterLayers):
    zCoor = getZvalues(inFile)
    Layers, numLayers = getLayers(inFile)
    idx, vals = compEmpty(zCoor, numInterLayers)
    x = np.linspace(0,1,10)
    cPoints = getCpoints(inFile)
    volume = np.array([])
    for i in range(len(cPoints[0])):
        segment = np.array([])

        for r in range(len(cPoints)):
            B = deCastel(cPoints[r][i],x)
            #print(B)
            if r == 0: 
                segment = B
            else:
                segment = np.concatenate((segment,B))
            #segment = np.append([B],segment)
        #segment = np.reshape(segment,((numLayers-1)*numInterLayers,2))
        
        xCoor, yCoor = Inter(segment, numInterLayers, inFile, x)
        volume = np.append([volume], re)
        #volume = np.concatenate((volume, re))
    #volume = np.reshape(volume,(numLayers+numInterLayers, (i+1)*len(x)))
        #Volume = Inter(segment, numInterLayers)
    return volume

def interInter(cPoints, interZ, numPoints):
    
    Inter = interZ[1]- interZ[0]
    Inter = Inter/numPoints
    zValues = np.linspace(interZ[0]+ Inter, interZ[-1], numPoints, False)
    
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