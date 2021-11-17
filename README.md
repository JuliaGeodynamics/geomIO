# geomIO - creating 3D geometries from 2D input #

[![Build Status](https://github.com/JuliaGeodynamics/geomIO/workflows/CI/badge.svg)

This is the Julia and Python version of [geomIO](https://geomio.bitbucket.io), a free open source software to generate 3D volumes and surfaces from multiple 
2D cross sections. These are provided as .stl files which can be used e.g. in a 3D printer, or as input to geodynamic codes.

This is currently work-in-progress and not intended to be used yet. The python version is more mature and `run.py` in the /src directory can run and explain the code.


### Input ###

Input is an .svg file generated with Inkscape. Currently only one structure can be handled per file.
To generate a valid input file follow the following steps:

- Add a layer and rename it to a numerical value (representing the 3rd dimension coordinates)
- Draw your structure with the Beziertool (Shift + F6)
- Press F2 and select all nodes (ctrl +A), then press "Make selected nodes smooth" (top bar)
- Go back to the layer menu (shift + ctrl + L) and duplicate your layer via rightclick
- Rename the duplicated layer to another numerical value
- Use F2 to shift the control points according to your idea
- IMPORTANT: every layer must contain the same number of control points
- If you use a stencil layer you can begin its name with "$" which makes geomIO ignore the layer (this can lead to problems though, WIP)

### Getting started ###

To get started open the run.py file and edit it according to your needs.
The following paramters can be set:

- inFile: input .svg file
- name: name of the output file(must end with .stl)
- Volume : bool, wether your structure is a closed volume or an open surface
- numInterLayers: number of layers to interpolate between your given inkscape layers (larger number = finer mesh)
- nPrec = number of points that are computet per Bezier segment (larger number = finer mesh)

After setting your paramters you can call the geomioFront function, which does the rest for you and saves the output file.
If you have Open3D installed you can visualize your pointcloud directly in Python with plotCloud3D.
Since there is no scaling tool implemented yet you can get your structer boundary coordinates with getBounds.
For a test workflow you can use one of the .svg files from the input directory and test geomio4py with them.
Please note that geomio4py is still under development and various issues may occur.



