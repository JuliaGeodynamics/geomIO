#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 14:32:22 2021

@author: lucas
"""
import sys, os
import geomio as gm

#name of the inputfile
inFile = "input/over.svg"

#name of the stl file which will be generated
name = 'mesh.stl'

#wether it is a closed volume or an open surface
Volume = False

#how many layers are interpolated in total
numInterLayers = 25

#number of points that are computed per bezier segment
nPrec = 50

#calling the main function
gm.geomioFront(inFile,numInterLayers, nPrec, name, Volume)

#plot the pointcloud. requiers open3d
#gm.plotCloud3D(inFile,numInterLayers,nPrec)

#get the voundaries of the mesh/volume
gm.getBounds(inFile, numInterLayers, nPrec)