using Test, geomIO



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


# Closed surface
inFile          = "./input/volume.svg"
tri     = triSurfClose(inFile, numInterLayers, nPrec)
@show  tri[100][1][2] ≈ -156.21323f0