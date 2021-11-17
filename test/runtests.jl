using Test, geomIO





inFile          = "./input/over.svg"
outFile         = "folds2.stl"
Volume          = false
numInterLayers  = 4
nPrec           = 40;
geomioFront(inFile,outFile, numInterLayers, nPrec, Volume)
@test isfile("./output/$outFile")

# number of layers
l,nl = getLayers(inFile)
@test nl==3

# construct triangular stl mesh from *.,svg file & read it into julia
tri = triSurfOpen(inFile, numInterLayers, nPrec)
@test  tri[200][1][1] ≈ 127.50775f0