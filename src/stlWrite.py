#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 14:57:32 2021

@author: lucas
"""
import struct
import numpy

def stlWrite(ver):
    class Stl(object):
        dtype = numpy.dtype([
            ('normals', numpy.float32, (3, )),
            ('v0', numpy.float32, (3, )),
            ('v1', numpy.float32, (3, )),
            ('v2', numpy.float32, (3, )),
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
                data = numpy.fromfile(fh, dtype=cls.dtype, count=size)
                return Stl(header, data)
    
        def to_file(self, filename, mode='wb'):
            with open(filename, mode) as fh:
                fh.write(self.header)
                fh.write(struct.pack('@i', self.data.size))
                self.data.tofile(fh)
    file = Stl()
    return


# if __name__ == '__main__':
#     # Read from STL file
#     stl = Stl.from_file('test.stl')

#     # Increment the X axis by one
#     stl.data['v0'][:, 0] += 1
#     stl.data['v1'][:, 0] += 1
#     stl.data['v2'][:, 0] += 1

#     # Write back to file
#     stl.to_file('test.stl')