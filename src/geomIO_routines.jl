# This provides interfaces to the various python routines

export geomioFront

"""
    geomioFront(inFile::String, outFile::String, numInterLayers::Int64, nPrec::Int64, Volume::Bool)

Main routine that reads an *.svg file and creates an `*.stl` triangulates surface out of it,

- `inFile`  : Name of the input `*.svg` file (including directory)
- `outFile` : Name of the `*.stl` file which will be generated
- `numInterLayers` : The number of internal layers in the file that will be interpreted
- `nPrec` : Number of points that are computed per bezier segment
- `Volume` : Boolean that indicates whether we are reading a closed volume or nothing

The resulting file will be written to a file in the directory `./output` (relative to the current directory)

"""
function geomioFront(inFile::String, outFile::String, numInterLayers::Int64, nPrec::Int64, Volume::Bool)

    # Check whether an output directory and create one if it does not exist
    if !isdir("output")
        mkdir("output")
    end

    # Call python code
    pygeomio.geomioFront(inFile,numInterLayers, nPrec, outFile, Volume)

    # output
    println("Wrote ./output/$(outFile)")
    
    return nothing 
end