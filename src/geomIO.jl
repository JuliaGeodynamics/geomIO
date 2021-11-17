module geomIO

# Try forcing using conda python installation
using  Conda, PyCall
const pygeomio = PyNULL()

function __init__()

  # Required python packages (installs them in the local conda installation if possible)
  module_name   = ["svgpathtools","numpy","stl","matplotlib","scipy","ipdb","vtk"]
  packages      = ["svgpathtools","numpy","numpy-stl","matplotlib","scipy","ipdb","vtk"]

  if PyCall.conda
    println("Employing the build-in julia conda python distribution for geomIO")
    
    # Import packages
    chn           =   "conda-forge"
    for i=1:length(module_name)
      pyimport_conda(module_name[i], packages[i],  chn)
    end
   
  else  
    println("Employing an existing python distribution for geomIO, located in $(PyCall.libpython).")
    println("Ensure that the following packages are installed there: $(packages)")
    println("Alternatively, you can install manually instruct PyCall to use the build-in conda distribution with:")
    println("julia> using Pkg")
    println("julia> ENV[\"PYTHON\"] =\"\"")
    println("julia> Pkg.add(\"PyCall\")")
    println("Restart julia after this")
    
    # Import packages
    for i=1:length(module_name)
      pyimport(module_name[i], packages[i])
    end

  end

  pushfirst!(PyVector(pyimport("sys")."path"), (@__DIR__)*"/python")  # relative to /src of geomIO
  copy!(pygeomio, pyimport("geomio"))

end
export pygeomio

include("./geomIO_routines.jl")


end # module
