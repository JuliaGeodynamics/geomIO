
#=
using Test, geomIO, PythonCall



# 
inFile          = "./input/over.svg"
outFile         = "folds2.stl"
Volume          = false
numInterLayers  = 4
nPrec           = 40;

# read & write
geomioFront(inFile,outFile, numInterLayers, nPrec, Volume)
@test isfile("./output/$outFile")

# number of layers
l,nl = getLayers(inFile)
@test nl==3

# Construct *.stl mesh from *.svg file & read it into julia
tri = triSurfOpen(inFile, numInterLayers, nPrec)
@test  tri[200][1][1] ≈ 127.50775f0

# Example in which we don't name the layer "m300/p300", but "m100/300" 
inFile          = "./input/fold.svg"
tri = triSurfOpen(inFile, numInterLayers, nPrec)
@test  tri[200][1][3] ≈ -100.0f0


# Closed surface
inFile          = "./input/volume.svg"
tri     = triSurfClose(inFile, numInterLayers, nPrec)
@test  tri[100][1][2] ≈ -156.21323f0

# Slab example
inFile          = "./input/slab.svg"
tri     = triSurfOpen(inFile, numInterLayers, nPrec)
@test  tri[100][1][2] ≈ -94.11989f0

# Example
inFile          = "./input/over.svg"
tri     = triSurfOpen(inFile, numInterLayers, nPrec)
@test  tri[100][1][2] ≈ -62.098507f0

# read Inkscape SVG file with multiple lines
inFile      = "./input/over2.svg"       # inkscape file
svgFileData = pygeomio.readSVG(inFile, false)
@test pyconvert(String,svgFileData[0][1]) == "NewCurve"
@test pyconvert(Complex,svgFileData[1][1][1][1]) ≈ 143.15726999999998 + 116.18955im
@test pyconvert(String,svgFileData[2][4]) == "Reference"
@test pyconvert(Int64,svgFileData[3]) == 5
@test pyconvert(Vector,svgFileData[5])==[0.0,0.0,0.0,100.0,nothing,nothing]
@test pyconvert(Vector,svgFileData[6])[3] ≈ 10.523572

# read Affinity Design SVG file
inFile      = "./input/example_AffinityDesigner.svg"       # inkscape file
svgFileData = pygeomio.readSVG(inFile, false)
@test pyconvert(String,svgFileData[0][1]) == "Curve_example"
@test pyconvert(Complex,svgFileData[1][1][1][1]) ≈ 647.662 + 1576.46im
@test pyconvert(String,svgFileData[2][4]) == "Reference"
@test pyconvert(Int64,svgFileData[3]) == 3
@test pyconvert(Vector,svgFileData[5])==[0.0,100.0,100.0,100.0,nothing]
@test pyconvert(Vector,svgFileData[6])[3] ≈ 267.841
=#