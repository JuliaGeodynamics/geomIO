module geomIO

# Try forcing using conda python installation
using  Conda, PyCall
const pygeomio = PyNULL()

function __init__()

  # Load the required python packages (installs them in the local conda installation if required)
  module_name   = ["svgpathtools","numpy","stl","matplotlib","scipy","ipdb","vtk"]
  packages      = ["svgpathtools","numpy","numpy-stl","matplotlib","scipy","ipdb","vtk"]

  if PyCall.conda
    println("Employing the build-in julia conda python distribution")
    
    # Import packages
    chn           =   "conda-forge"
    for i=1:length(module_name)
      pyimport_conda(module_name[i], packages[i],  chn)
    end
   
  else  
    println("Employing the build-in python distribution, located in $(PyCall.libpython).")
    println("Ensure that the following packages are installed there: $(packages)")
    
    # Import packages
    for i=1:length(module_name)
      pyimport(module_name[i], packages[i],  chn)
    end

  end

  pushfirst!(PyVector(pyimport("sys")["path"]), (@__DIR__)*"/python")  # relative to /src of geomIO
  copy!(pygeomio, pyimport("geomio"))

end
export pygeomio

include("./geomIO_routines.jl")


end # module
