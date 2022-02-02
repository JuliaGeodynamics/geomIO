module geomIO

using PythonCall
const pygeomio = PythonCall.pynew()

function __init__()

  pth = (@__DIR__)*"/python"        # Path where the python routines are
  pyimport("sys").path.append(pth)  # append path

  # link geomio. Note that all python dependencies are listed in PythonCallDeps.toml
  PythonCall.pycopy!(pygeomio, pyimport("geomio")) 

end

export pygeomio

include("./geomIO_routines.jl")


end # module

