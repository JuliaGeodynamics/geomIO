
from curve_interpolations import *
from svgpathtools import svg2paths, real, imag, Line, svg2paths2, Document
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as Pt
import matplotlib.patches as patches
import math
import sys,os
import scipy as sc
from scipy import interpolate
#from numba import jit


#@jit(nopython=True) 
def indexingOpen(numLayers, lc, nInter):
    nQuads = nInter+numLayers -1 # number of quads to be split in two triangles
    tri1 = np.array([])#first set of quadhalfs
    tri2 = np.array([])#second set of quadhalfs
    for p in range(nQuads):
        numP = int(len(lc)/(nQuads+1))#number of ???
        nodeIdx = numP *p
        shape = (int(numP*2)-2,3)
        shapeS = (numP,3)
        ind = np.zeros(shapeS, dtype = int)
        ind2 = np.zeros(shapeS,dtype=int)
        arrshift = numP-1
        counter = 0
        if p == 0:
            for i in range(numP-1):
                if i == 0:
                    continue
                ind[i,0] = i
                ind[i,1] = i+1
                ind[i,2] = numP+i
                
                ind2[i,0] = numP +i
                ind2[i,1] = numP + i+1
                ind2[i,2] = i+1
                counter +=1
            tri1 = ind
            tri2 = ind2
        else:
    
            for i in range(numP-1):
                if i == 0:
                    continue
                ind[i,0] = nodeIdx + i
                ind[i,1] = nodeIdx + i+1
                ind[i,2] = nodeIdx + arrshift+i
                
                ind2[i,0] = nodeIdx + arrshift +i
                ind2[i,1] = nodeIdx + arrshift + i+1
                ind2[i,2] = nodeIdx+ i+1
                counter +=1
                
            tri2 = np.concatenate((ind2, tri2))    
            tri1 = np.concatenate((ind, tri1))
    
    triangles = np.concatenate((tri1,tri2))
    normals = triNormals(triangles, lc)

    mesh1 = np.array([])
    #print(counter)
    
    
            
    for p in range(len(triangles)):
        #first ind = first tri
        mesh1 = np.concatenate((mesh1, normals[p]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,0]]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,1]]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,2]]))
        
    mesh1 = np.reshape(mesh1,(len(triangles)*4,3))
    # CoorR = np.zeros_like(mesh1)
    # CoorR[0,:] = LayerCoors[1,:]
    # CoorR[1,:] = LayerCoors[2,:]
    # CoorR[2,:] = LayerCoors[0,:]
    
    
    return mesh1, triangles

      

#@jit(nopython=True) 
def triNormals(tri,lc):
    """
    Function to compute the normal of a triangle
    currently unused
    """
    normals = np.zeros((len(tri),3))
    for i in range(len(tri)):
        #print(i)
        p1 = lc[tri[i,0],:]
        p2 = lc[tri[i,1],:]
        p3 = lc[tri[i,2],:]
        norm = np.cross(p2-p1, p3-p1)
        normals[i,:] = norm
    return normals
def triSurfOpen(lc, numLayers, nInter, nPrec):
    """
    surface triangulation for open surfaces
    """
    #Layers, numLayers = getLayers(inFile)
    
    #lc = getPointCoords(inFile,nInter,nPrec) # get coordinates
    
    m1, t = indexingOpen(numLayers, lc, nInter)
    
    return m1, t
    
   

from scipy.spatial import ConvexHull, convex_hull_plot_2d,Delaunay 
    
#@jit(nopython=True) 
def layerTri(lc, numP, LLI, last = True):
    

    
    if (numP % 2) == 0:
        shape = (int((numP-2)/2),3)
        index1 = np.zeros(shape, dtype = int)
        index2 = np.zeros(shape, dtype = int)
        for i in range(len(index1)):
            index1[i,0] = i
            index1[i,1] = i+1
            index1[i,2] = numP-1-i
            
            index2[i,0] = i
            index2[i,1] = numP-1-i
            index2[i,2] = numP -2 -i
    
        triangles = np.concatenate((index1,index2))

    else:
        shape = (int((numP-3)/2),3)
        index1 = np.zeros(shape, dtype = int)
        index2 = np.zeros(shape, dtype = int)
        for i in range(len(index1)):
            index1[i,0] = i+1
            index1[i,1] = i+2
            index1[i,2] = numP-1-i
            
            index2[i,0] = i+1
            index2[i,1] = numP-i -1
            index2[i,2] = numP -2 -i
    
        first = np.array([[0,1,numP-1]])
        triangles = np.concatenate((index1,index2))
        triangles = np.concatenate((triangles, first))
        
    normals = triNormals(triangles, lc)

    mesh = np.array([])
    for p in range(len(triangles)):
        #first ind = first tri
        mesh = np.concatenate((mesh, normals[p]))
        mesh = np.concatenate((mesh, lc[triangles[p,0]]))
        mesh = np.concatenate((mesh, lc[triangles[p,1]]))
        mesh = np.concatenate((mesh, lc[triangles[p,2]]))
        
    mesh = np.reshape(mesh,(len(triangles)*4,3))
   
    if last and (numP % 2) == 0:
        trueIdx = np.concatenate((index1 + LLI, index2 + LLI))
        triangles = trueIdx
    elif last:
        trueIdx = np.concatenate((index1 + LLI, index2 + LLI))
        triangles = np.concatenate((trueIdx, first + LLI))

        
            
    return mesh, triangles
#@jit(nopython=True) 
def triangulateCover(lc):
    numP = len(lc)
    shape = len(lc)-2
    ind = np.zeros((shape,3))
    triangles = np.zeros((shape,3), dtype = int)
    for i in range(shape):
        triangles[i,0] = i 
        triangles[i,1] = i+1
        triangles[i,2] = numP-2-i
    normals = triNormals(triangles, lc)
    mesh = np.array([])
    for p in range(len(triangles)):
        #first ind = first tri
        mesh = np.concatenate((mesh, normals[p]))
        mesh = np.concatenate((mesh, lc[triangles[p,0]]))
        mesh = np.concatenate((mesh, lc[triangles[p,1]]))
        mesh = np.concatenate((mesh, lc[triangles[p,2]]))
        
    mesh = np.reshape(mesh,(len(triangles)*4,3))
    
    return  mesh, triangles

def triSurfClose(lc, numLayers, nInter, nPrec):
    """
    surface triangulation for closed volumes
    """
    #Layers, numLayers = getLayers(inFile)
    #lc = getPointCoords(inFile,nInter,nPrec)

    nQuads = nInter+numLayers -1
    tri1 = np.array([])
    tri2 = np.array([])
    
    #========Section where the "covers" are done : faulty
    #number of points per layer
    numP = int(len(lc)/(nQuads+1))
    #get x and y coordinates from first and last layer
    #cov1 = lc[0:numP,0:2]
    cov1Z = lc[0:numP,:]
    #last layer index, start of last layer
    LLI = len(lc)-numP
    #cov2 = lc[LLI:-1,0:2]
    cov2Z =lc[LLI:-1,:]
    dummy = np.array([lc[len(lc)-1]])
    cov2Z = np.concatenate((cov2Z,dummy))
    # startCov = Delaunay(cov1,furthest_site=True)
    # endCov = Delaunay(cov2, furthest_site=True)
    
    #cCF, cTF = layerTri(cov1Z,numP,LLI, False)
    #cCL, cTL = layerTri(cov2Z, numP, LLI, True)
    cCF, cTF = triangulateCover(cov1Z)
    cCL, cTL = triangulateCover(cov2Z)
    cTL += LLI
    #coordinatses of both covers
    cover = np.concatenate((cCF,cCL))
    #index of the triangles in lc
    covTri = np.concatenate((cTF,cTL))

    

    #=======section for inside: working    

    for p in range(nQuads):

        nodeIdx = numP *p
        shape = (int(numP*2)-2,3)
        shapeS = (numP,3)
        ind = np.zeros(shapeS, dtype = int)
        ind2 = np.zeros(shapeS,dtype=int)
        arrshift = numP-1
        counter = 0
        if p == 0:
            for i in range(numP):

                ind[i,0] = i
                ind[i,1] = i+1
                ind[i,2] = numP+i
                
                ind2[i,0] = numP +i
                ind2[i,1] = numP + i+1
                ind2[i,2] = i+1
                counter +=1
            tri1 = ind
            tri2 = ind2
        else:
    
            for i in range(numP):

                ind[i,0] = nodeIdx + i
                ind[i,1] = nodeIdx + i+1
                ind[i,2] = nodeIdx + arrshift+i
                
                ind2[i,0] = nodeIdx + arrshift +i
                ind2[i,1] = nodeIdx + arrshift + i+1
                ind2[i,2] = nodeIdx+ i+1
                counter +=1
                
            tri2 = np.concatenate((ind2, tri2))    
            tri1 = np.concatenate((ind, tri1))
    
    triangles = np.concatenate((tri1,tri2))
    normals = triNormals(triangles, lc)

    mesh1 = np.array([])
    #print(counter)
    
    
            
    for p in range(len(triangles)):
        #first ind = first tri
        mesh1 = np.concatenate((mesh1, normals[p]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,0]]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,1]]))
        mesh1 = np.concatenate((mesh1, lc[triangles[p,2]]))
        #vertex = Tri(lc[triangles[p,0]],lc[triangles[p,1]],lc[triangles[p,1]])
        
    mesh1 = np.reshape(mesh1,(len(triangles)*4,3))
    mesh1 = np.concatenate((mesh1,cover))
    triangles = np.concatenate((triangles,covTri))

    return mesh1, triangles







def wSTL(data, path, numInter, nPrec,name, volume = False, mode = "ASCII"):
    """
    placeholder function for writing stl files
    uses the lib np stl
    """
    
    import struct
    import stl
    from stl import mesh
    lc = getPointCoords(data, path, numInter,nPrec)

    if volume:
        triangles, face = triSurfClose(lc,len(path), numInter,nPrec)
    else:
            
        triangles, face = triSurfOpen(lc, len(path), numInter,nPrec)
    #nec for writing
    #os.chdir("output")
    cube = mesh.Mesh(np.zeros(face.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(face):
        for j in range(3):
            cube.vectors[i][j] = lc[f[j],:]
    
    print("Saved file " + str(name) + " to directory:")
    print(os.getcwd())
    #Write the mesh to file "cube.stl"
    if mode == "ASCII":
        
        cube.save(str(name),mode=stl.Mode.ASCII )
    elif mode == "BIN":
        cube.save(str(name),mode=stl.Mode.BINARY )
    #else:
    #    print("Please set writing mode to BIN or ASCII")
    #os.chdir("..")
    return


def wMulti(inFile, numInter, nPrec, names, volume):
    import struct
    import stl
    from stl import mesh
    

class Stl(object):
    dtype = np.dtype([
        ('normals', np.float32, (3, )),
        ('v0', np.float32, (3, )),
        ('v1', np.float32, (3, )),
        ('v2', np.float32, (3, )),
        ('attr', 'u2', (1, )),
    ])

    def __init__(self, header, data):
        self.header = header
        self.data = data

    @classmethod
    def from_file(cls, filename, mode='rb'):
        with open(filename, mode) as fh:
            header = fh.read(80)
            size, = struct.unpack('@i', fh.read(4))
            data = np.fromfile(fh, dtype=cls.dtype, count=size)
            return Stl(header, data)

    def to_file(self, filename, mode='wb'):
        with open(filename, mode) as fh:
            fh.write(self.header)
            fh.write(struct.pack('@i', self.data.size))
            self.data.tofile(fh)

# head = "head"
# head = bytes(head, 'utf-8')
# file = Stl(head,mesh)
# file.to_file('mesh.stl')

    
def write(ver):
    header = 0
    with open("slab.stl", mode="wb")as fh:
#        writer = wr.Binary_STL_Writer(fh)
 #       writer.add_faces(ver)
  #      writer.close()
         fh.write(header)
         fh.write(struct.pack)
  
    return
# sim = dlnSC.simplices
# vertex = dlnSC.vertices
# #vertex = lc[ind]
# from stl import mesh
# volume = mesh.Mesh(np.zeros(sim.shape[0], dtype=mesh.Mesh.dtype))
# for i, f in enumerate(sim):
#     for j in range(3):
#         volume.vectors[i][j] = vertex[f[j],:]
# volume.save('slab3.stl')

# dln = delny3D(lc)

# head = "head"
# head = bytes(head, 'utf-8')
# file = Stl(head,dlnSC.points)
# file.to_file('rick.stl')

# test = Stl.from_file('rick.stl')
# from vtk.util.np_support import vtk_to_np
# np_coordinates = vtk_to_np(dln)

# reader = vtk.vtkXMLUnstructuredGridReader()
# reader.SetFileName( "slabvtk.vtu" )
# reader.Update()


# Point_cordinates = reader.GetOutput().GetPoints().GetData()

# #write(vertex)
#     #file = stlWrite(dln.vertices)
# # fig = plt.figure()
# # ax = plt.axes(projection='3d')
# # ax.scatter3D(lc[:,0], lc[:,1], lc[:,2])
# # plt.show()

# # #dln.write()
# writer = vtk.vtkXMLUnstructuredGridWriter()
# #writer.SetFileType("stl")
# writer.SetFileName('slabvtk.vtu')
# writer.SetInputData(dln)
# writer.Write()

#pcd = o3d.geometry.PointCloud()
#pcd.points = o3d.utility.Vector3dVector(lc[:,:3])
#pcd.colors = o3d.utility.Vector3dVector(lc[:,3:6]/255)
#pcd.normals = o3d.utility.Vector3dVector(lc[:,6:9])
#o3d.visualization.draw_geometries([pcd])
#hull = pcd.compute_convex_hull()
#----------pivoting method
# distances = pcd.compute_nearest_neighbor_distance()
# avg_dist = np.mean(distances)
# radius = 3 * avg_dist

# #bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd,o3d.utility.DoubleVector([radius, radius * 2]))
