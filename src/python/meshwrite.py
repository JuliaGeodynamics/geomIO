# File: stlwrite.py
# Purpose: Export 3D objects built of faces with 3 or 4 vertices, as
# an ASCII or Binary STL file.
# License: MIT License
#
# Original source from http://code.activestate.com/recipes/578246-stl-writer/
# Modified by David Touretzky, November 2014, to compute surface normals.

import math
import struct

ASCII_FACET = """facet normal {normal[0]:.4f} {normal[1]:.4f} {normal[2]:.4f}
outer loop
vertex {face[0][0]:.4f} {face[0][1]:.4f} {face[0][2]:.4f}
vertex {face[1][0]:.4f} {face[1][1]:.4f} {face[1][2]:.4f}
vertex {face[2][0]:.4f} {face[2][1]:.4f} {face[2][2]:.4f}
endloop
endfacet
"""

BINARY_HEADER ="80sI"
BINARY_FACET = "12fH"

def crossproduct(u,v):
    s1 = u[1]*v[2] - u[2]*v[1]
    s2 = u[2]*v[0] - u[0]*v[2]
    s3 = u[0]*v[1] - u[1]*v[0]
    return [s1, s2, s3]

def norm(v):
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def normalize(v):
    m = norm(v)
    if m == 0:
        m = 1
    return [v[0]/m, v[1]/m, v[2]/m]

def normalto(u,v):
    return normalize(crossproduct(u,v))

def diff(u,v):
    return [u[0]-v[0], u[1]-v[1], u[2]-v[2]]

class ASCII_STL_Writer:
    """ Export 3D objects build of 3 or 4 vertices as ASCII STL file.
    """
    def __init__(self, stream):
        self.fp = stream
        self._write_header()

    def _write_header(self):
        self.fp.write("solid python\n")

    def close(self):
        self.fp.write("endsolid python\n")

    def _write(self, face):
        n = normalto(diff(face[0],face[1]),diff(face[1],face[2]))
        self.fp.write(ASCII_FACET.format(normal=n,face=face))

    def _split(self, face):
        p1, p2, p3, p4 = face
        return (p2, p1, p3), (p3, p1, p4)

    def add_face(self, face):
        """ Add one face with 3 or 4 vertices. """
        if len(face) == 4:
            face1, face2 = self._split(face)
            self._write(face1)
            self._write(face2)
        elif len(face) == 3:
            self._write(face)
        else:
            raise ValueError('only 3 or 4 vertices for each face')

    def add_faces(self, faces):
        """ Add many faces. """
        for face in faces:
            self.add_face(face)

class Binary_STL_Writer(ASCII_STL_Writer):
    """ Export 3D objects build of 3 or 4 vertices as binary STL file.
    """
    def __init__(self, stream):
        self.counter = 0
        ASCII_STL_Writer.__init__(self,stream)

    def close(self):
        self._write_header()

    def _write_header(self):
        self.fp.seek(0)
        self.fp.write(struct.pack(BINARY_HEADER, b'Python Binary STL Writer', self.counter))

    def _write(self, face):
        self.counter += 1
        n = normalto(diff(face[0],face[1]),diff(face[1],face[2]))
        data = [
            n[0], n[1], n[2],
            face[0][0], face[0][1], face[0][2],
            face[1][0], face[1][1], face[1][2],
            face[2][0], face[2][1], face[2][2],
            0
        ]
        self.fp.write(struct.pack(BINARY_FACET, *data))

# def example():
#     def get_triangles():
#         p1 = (0, 0, 0)
#         p2 = (1, 2, 3)
#         p3 = (-5, -7, 4)
#         p4 = (3, -2, 6)
#         return [ [p1, p3, p2], [p2, p3, p4], [p4, p3, p1], [p4, p1, p2] ]
#     filename = 'triangles.stl'
#     with open(filename, 'wb') as fp:
#         writer = ASCII_STL_Writer(fp)
#         writer.add_faces(get_triangles())
#         writer.close()
#     print 'Wrote ' + filename
