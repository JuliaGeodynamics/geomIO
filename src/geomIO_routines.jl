# This provides interfaces to the various python routines

using MeshIO, GeometryBasics    # STL in julia

export save, load # load and save STL meshes

export geomioFront, getBounds, getLayers, triSurfOpen, triSurfClose, svg2paths

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


"""
    getBounds(inFile::String, numInterLayers::Int64, nPrec::Int64)

Prints the boundaries of the mesh or volume
"""
function getBounds(inFile::String, numInterLayers::Int64, nPrec::Int64)
    
    return pygeomio.getBounds(inFile, numInterLayers, nPrec)
end

"""
    l,nl = getLayers(inFile::String)

Returns the number of layers `nl` in the file as well as the path's per layer, `l`    
"""
function getLayers(inFile::String)
    l,nl = pygeomio.getLayers(inFile)

    return l, pyconvert(Int64,nl)
end


"""
    tri = triSurfOpen(inFile::String, numInterLayers::Int64, nPrec::Int64)

Reads the SVG file `inFile` and reconstructs a triangular surface from it. The surface is *not* closed.  
`tri` is a triangular `stl` mesh that is consistent with `MeshIO`, and `FileIO`

# Example
=========
```julia
julia> using geomIO
julia> inFile="./test/input/slab.svg"
julia> numInterLayers, nPrec = 4,40
julia> tri = triSurfOpen(inFile, numInterLayers, nPrec)
Mesh{3, Float32, Triangle}:
 Triangle(Float32[241.55135, -79.856834, -100.0], Float32[241.55135, -79.856834, -100.0], Float32[241.55135, -79.856834, -100.0])
 Triangle(Float32[241.71376, -79.992004, 162.5], Float32[241.86136, -80.17122, 162.5], Float32[241.55135, -79.856834, 250.0])
 Triangle(Float32[241.86136, -80.17122, 162.5], Float32[241.99454, -80.39218, 162.5], Float32[241.71376, -79.992004, 250.0])
 Triangle(Float32[241.99454, -80.39218, 162.5], Float32[242.11372, -80.65258, 162.5], Float32[241.86136, -80.17122, 250.0])
 Triangle(Float32[242.11372, -80.65258, 162.5], Float32[242.2193, -80.95012, 162.5], Float32[241.99454, -80.39218, 250.0])
 ⋮
 Triangle(Float32[237.53671, -78.19536, 31.25], Float32[239.08516, -78.67469, 31.25], Float32[239.08516, -78.67469, -100.0])
 Triangle(Float32[239.08516, -78.67469, 31.25], Float32[240.4272, -79.22747, 31.25], Float32[240.4272, -79.22747, -100.0])
 Triangle(Float32[240.4272, -79.22747, 31.25], Float32[241.55135, -79.856834, 31.25], Float32[241.55135, -79.856834, -100.0])
 Triangle(Float32[241.55135, -79.856834, -100.0], Float32[241.55135, -79.856834, -100.0], Float32[241.55135, -79.856834, -100.0])
```
You can save this with:
```julia
julia> save("surf.stl",tri)
```

"""
function triSurfOpen(inFile::String, numInterLayers::Int64, nPrec::Int64)

    mesh, triangles = pygeomio.triSurfOpen(inFile,numInterLayers,nPrec)
    
    # reconstruct the normals to the triangles
    lc          =   pygeomio.getCarthesian(inFile,numInterLayers,nPrec)
    normals     =   pygeomio.triNormals(triangles, lc)
    
    # Load points & vertexes
    points      = connect(pyconvert(Array,lc),Point{3})

    # Create stl surface     
    tri = CreateSTL(pyconvert(Array,triangles), points, pyconvert(Array,normals))

    return tri
end

"""
    tri =  triSurfClose(inFile::String, numInterLayers::Int64, nPrec::Int64)

Reads the SVG file `inFile` and reconstructs a triangular surface from it. The surface is closed.  
The output `tri` is a triangular `stl` mesh, consistent with `MeshIO`. 

# Example
```julia
julia> 
```

"""
function  triSurfClose(inFile::String, numInterLayers::Int64, nPrec::Int64)

    mesh, triangles = pygeomio.triSurfClose(inFile,numInterLayers,nPrec)
    
    # reconstruct the normals to the triangles
    lc          = pygeomio.getCarthesian(inFile,numInterLayers,nPrec)
    normals     = pygeomio.triNormals(triangles, lc)
    
    # Load points & vertexes
    points      = connect(pyconvert(Array,lc),Point{3})

    # Create stl surface     
    tri = CreateSTL(pyconvert(Array,triangles), points, pyconvert(Array,normals))

    return tri
end

"""
    paths, attributes = svg2paths(inFile::String) 

Retrieves the `paths` and `atrributes` of the svg file.       
"""
function svg2paths(inFile::String)
    paths, attributes = pygeomio.svg2paths(inFile)

    return paths, attributes 
end


# Internal routine to create an STL surface from the python info
function CreateSTL(triangles, points, normalspy)

    # Reinterpolate the mesh 
    triangle_count = size(triangles,1)
    faces           = Array{GLTriangleFace}(undef, triangle_count)
    vertices        = Array{Point3f}(undef, triangle_count * 3)
    normals         = Array{Vec3f}(undef, triangle_count * 3)
    for i=1:triangle_count
        j=i-1;
        faces[j+1] = GLTriangleFace(j * 3 + 1, j * 3 + 2, j * 3 + 3)

        normals[j*3+1] = Vec3f(normalspy[i,:])
        normals[j*3+2] = normals[j*3+1] # hurts, but we need per vertex normals
        normals[j*3+3] = normals[j*3+1]

        vertices[j*3+1] = points[triangles[i,1] .+ 1]
        vertices[j*3+2] = points[triangles[i,2] .+ 1]
        vertices[j*3+3] = points[triangles[i,3] .+ 1]
    end

    # Create triangular mesh (STL format)
    return tri  = Mesh(meta(vertices; normals=normals), faces)

end
