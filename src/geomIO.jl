module geomIO

ENV["PYTHON"]=""
using PyCall, Conda

# Load the required python packages (installs them in the local conda installation if required)
chn="conda-forge"
Conda.add("svgpathtools", :my_env, channel=chn)
Conda.add("numpy", :my_env, channel=chn)
Conda.add("numpy-stl", :my_env, channel=chn)
Conda.add("matplotlib", :my_env, channel=chn)
#Conda.add("math", :my_env)
#Conda.add("sys", :my_env)
#Conda.add("os", :my_env)
Conda.add("scipy", :my_env, channel=chn)
Conda.add("ipdb", :my_env, channel=chn)
Conda.add("vtk", :my_env, channel=chn)

#pyimport_conda("svgpathtools","svgpathtools")
#pyimport_conda("numpy","numpy")
#pyimport_conda("stl","numpy-stl")
#pyimport_conda("matplotlib","matplotlib")
#pyimport_conda("math","math")
#pyimport_conda("sys","sys")
#pyimport_conda("os","os")
##pyimport_conda("scipy","scipy")
#pyimport_conda("ipdb","ipdb")
#pyimport_conda("vtk","vtk")

pushfirst!(PyVector(pyimport("sys")."path"), "src/python")
pygeomio = pyimport("./src/geomio")

end # module
