#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 16:55:51 2021

@author: lucas
"""




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