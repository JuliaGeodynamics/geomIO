using Test, geomIO


inFile          = "./input/over.svg"
outFile         = "folds2.stl"
Volume          = false
numInterLayers  = 4
nPrec           = 40;
geomioFront(inFile,outFile, numInterLayers, nPrec, Volume)
@test isfile("./output/$outFile")

