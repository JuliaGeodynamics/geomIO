import numpy as np

class rectLinGrid:
    xBnds       : tuple
    yBnds       : tuple
    zBnds       : tuple
    nx          : int
    ny          : int
    nz          : int
    X           : np.ndarray
    Y           : np.ndarray
    Z           : np.ndarray
    PHS         : np.ndarray
    grd         : list
    
    def __init__(self,nx : int ,ny: int, nz: int, xBnds:tuple, yBnds:tuple, zBnds:tuple):
        self.nx              = nx
        self.ny              = ny
        self.nz              = nx
        self.xBnds           = xBnds
        self.xBnds           = yBnds
        self.xBnds           = zBnds
        self.x               = np.linspace(xBnds[0],xBnds[1],nx)
        self.y               = np.linspace(yBnds[0],yBnds[1],ny)
        self.z               = np.linspace(zBnds[0],zBnds[1],nz)
        self.X,self.Y,self.Z = np.meshgrid(self.x,self.y,self.z, indexing ='ij') #matrix indexing{ij} array[i,j]
        self.PHS             = np.zeros_like(self.X)
        self.grd             = list((self.X,self.Y,self.Z))