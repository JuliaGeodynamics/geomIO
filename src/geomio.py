# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import sys, os

from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM
import svgpathtools as pt


os.chdir("..")
#drawing = svg2rlg("input/vector.svg")
#x = drawing.getBounds()
f = open("input/vector.svg")
p = f.readlines()
x,y,z  = pt.svg2paths2("input/vector.svg")
x_1 = x[0]
b = x_1.point(0)

Bez= p[56]