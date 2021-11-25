# This contains various routines to read *.SVG files, with multiple layers and multiple curves
from tokenize import Number
from svgpathtools import svg2paths, real, imag, Line, svg2paths2, Document
from typing import NamedTuple
from curve_interpolations import *

class Scaling(NamedTuple):
    x0          : float
    y0          : float
    x0_SVG      : float
    y0_SVG      : float
    dx          : float
    dy          : float

class svgFileData(NamedTuple):
    CurveNames : list
    Curves     : list
    LayerNames : list
    numLayers  : int
    Commented  : list
    zCoord     : list
    Scaling    : Scaling
   
def readSVG(inFile, Verbose=True):
    """
    Reads an SVG file and returns all curves, layers and names (plus whether a layer is commented or not).
    The resulting layers/curves are sorted by depth and a scaling object is returned if a Reference layer was 
    present in the SVG file
    Has been tested for Inkscape & Affinity Design generated files.

    Parameters
    ----------
    inFile      : Input svg File
                    The .svg file used
    verbose     : Print output if True

    Returns
    -------
    Data        : Names Tuple with following fields: 
        CurveNames  : List with names of all curves found in the file (including on hidden layers)
        Curves      : List with curve paths of all curves found in file
        LayerNames  : List with layers on which the curves are
        numLayers   : Total number of layers
        Commented   : Boolean list that indicates whether the layer is commented out
        zCoord      : List with z-coordinates

    Note: all lists are sorted according to the z-coordinates (starting with the lowest value)   
    """

    # Read the names of the layers in case we have an inkscape file 
    LayerLabels, isInkscape = getLayerLabels_Inkscape(inFile)

    # Read all path's from file
    paths, attributes = svg2paths(inFile)
    
    # Initialize main output arrays
    LayerNames  = []
    CurveNames  = []
    Commented   = []
    Curves      = []
    
    # use svgpathtools to get the info from the file
    doc = Document(inFile)

    # Print names of layers:
    SVG_GROUP_TAG = 'svg:g'
    SVG_PATH_TAG = 'svg:path'
    SVG_NAMESPACE = {'svg': 'http://www.w3.org/2000/svg'}
    group = doc.tree.getroot()
    CurveNumber = 0
    numLayers= 0
    for elem in group.iterfind(SVG_GROUP_TAG, SVG_NAMESPACE):
        
        # 1) Extract name of layer.
        #    Note that Inkscape has a layer id, but the name you give it is stored under "inkscape:label" 
        #       Unfortunately, there appears to be a bug in svgpathtools (or a package it relies on), and 
        #       reading attributes that have an ":" in it does not work. 
        #       As a workaround, we use getLayerLabels_Inkscape to get the name of the Layer
        if elem.get('id')!=None:
            layer_str =  elem.get('id')

            if isInkscape:
                # In case we have inkscape the name of the label is stored under inkscape:label
                layer_str = LayerLabels.get(layer_str)
                
            numLayers += 1
        elif elem.get('label') != None:
            layer_str =  elem.get('label')
            numLayers += 1
        else:
            layer_str = None
        
        # 2) Determine if it is a commented layer or not
        comment_layer = False
        if layer_str!=None:
            if (layer_str[0]=="#") | (layer_str[0]=="$") | (layer_str=="Reference"):
                comment_layer = True

        # Print solution if requested
        if Verbose:    
            if comment_layer:
                print("Commented layer : " + layer_str )
            else:
                print("Found layer     : " + layer_str )
                
        # 3) Extract names of curves on the current layer  
        # if a curve is directly on this layer, you can find it like this
        for path in elem.iterfind(SVG_PATH_TAG, SVG_NAMESPACE):

            curve_str = None
            if path.get('id')!=None:
                curve_str =  path.get('id')

            if path.get('label') != None:
                curve_str =  path.get('label')
            
            if path.get('CoordRef') != None:    # in Inkscape we can specify CoordRef on a curve (not possible in AffinityDesign where we can only set id)
                curve_str =  "CoordRef_"+str(path.get('CoordRef'))
            
            # work-around for bug that doesn't find path.get('serif:id')
            lbl = attributes[CurveNumber].get('serif:id')
            if lbl != None:
                curve_str =  lbl
            
            # Store 
            LayerNames.append(layer_str)
            CurveNames.append(curve_str)
            Curves.append(paths[CurveNumber])
            Commented.append(comment_layer)
            CurveNumber += 1
            if Verbose:
                print("   Found PATH   : " + curve_str)
                
        # In some cases, Affinity Design puts another layer around the paths, so we go into that layer as well here
        for elem1 in elem.iterfind(SVG_GROUP_TAG, SVG_NAMESPACE):
            for path2 in elem1.iterfind(SVG_PATH_TAG, SVG_NAMESPACE):
                att = attributes[CurveNumber]   #
                lbl = att.get('serif:id')
                
                curve_str = None
                if elem1.get('id')!=None:
                    curve_str =  elem1.get('id')
                
                if elem1.get('label') != None:
                    curve_str =  elem1.get('label')
                
                lbl = attributes[CurveNumber].get('serif:id')
                if lbl != None:
                    curve_str =  lbl

                # Store info
                LayerNames.append(layer_str)
                CurveNames.append(curve_str)
                Curves.append(paths[CurveNumber])
                Commented.append(comment_layer)
                CurveNumber += 1
                if Verbose:
                    print("   Found PATH   : " + curve_str)
        
    if Verbose:
        print("Finished interpreting file : " + inFile + " --- ")

    # Filter names of curves, to take out "-"; 
    #   In inkscape, not specifying labels of curves will automatically append "-1","-2" etc. to the curve name
    #   This gets rid of these added parts. In generally, it is a better strategy to explicitly name the curves
    #   either by setting their id or by adding a label
    for i in range(len(CurveNames)):
        id = CurveNames[i].find('-')
        if id >= 0:
             CurveNames[i] =  CurveNames[i][0:id]

    # Interpret zLayers (get depth value from name)
    zCoord, sort_ind = get_zCoords(LayerNames, Commented)
    
    # Read Reference layer & determine scaling if it exist
    Scaling = get_Scaling(CurveNames, Curves, LayerNames)

    # Sort other lists accordingly, so it is all from lowest->highest coordinate 
    #  (with commented & reference layers listed afterwards)
    CurveNames_sort =   [CurveNames[i]  for i in sort_ind]
    Curves_sort     =   [Curves[i]      for i in sort_ind]
    LayerNames_sort =   [LayerNames[i]  for i in sort_ind]
    Commented_sort  =   [Commented[i]   for i in sort_ind]
    zCoord_sort     =   [zCoord[i]      for i in sort_ind]

    # Store data in Named Tuple (easier to handle later)
    Data = svgFileData(CurveNames_sort, Curves_sort, LayerNames_sort, numLayers, Commented_sort, zCoord_sort, Scaling)

    return Data


def getLayerLabels_Inkscape(inFile):
    """
    Parameters
    ----------
    inFile : Input svg File
        The .svg file used (can be inkscape or another file)

    Returns
    -------
    LayerLabels : dict
        Dictonary containing the name of the layer as well as the internal ID of the layer in case we have an inkscape *.svg file 
        If that is not the case the dict will be empty

    isInkscape : Bool
        Boolean that indicates if we have an inkscape SVG file or not    
    """

    LayerLabels = dict()  
    isInkscape  = False      
    
    f       = open(inFile)
    text    = f.readlines()
    
    # Step 1: determine the # of layers in the file (with the lines)
    index   = []
    numLayers = 0
    for i in range(len(text)):
        if "<g" in text[i]:     # start of a new layer
            numLayers += 1
            index.append(i)
    index.append(int(len(text)))

    #Step 2: Loop through each layer and extract curves that are present in the layer
    for iLayer in range(len(index)-1):
        layerLabel = None
        layerName  = None
        foundID    = False
        p = index[iLayer]
        for p in range(index[iLayer],index[iLayer+1]):    # loop over line-numbers of current layer
            
            # The label of the layer is indicated with "inkscape:layer" 
            if "inkscape:label=\""  in text[p]:
                layerSTR = text[p]
                layerSTR = layerSTR.split("\"")
                layerLabel = layerSTR[1]                 # Label of layer

            # Retrieve the ID of the layer. This is always specified before curves/paths, which is why we only read the 1th ID of a layer
            if ("id=\""  in text[p]) & (foundID==False):
                layerSTR = text[p]
                layerSTR = layerSTR.split("\"")
                layerName = layerSTR[1]                 # ID of layer
                foundID = True
                
        if (layerName!=None) & (layerLabel!=None):
            LayerLabels[layerName] = layerLabel
            isInkscape = True

    return LayerLabels, isInkscape


def get_zCoords(LayerNames, Commented):
    """
    Parameters
    ----------
    LayerNames  :   list
        Names of all the layers
    Commented   :   list
        Indicates whether the layer is commented or not   

    Returns
    -------
    zCoords : list
        List of z-values 
    """
    
    zCoord = []
    for i in range(len(LayerNames)):
     #   print(LayerNames[i] +"  "+ str(i))
        if (Commented[i]==True) :
            zCoord.append(None)    

        elif (LayerNames[i]=="Reference"):
            zCoord.append(None)    

        elif ((LayerNames[i][0]=="p") | (LayerNames[i][0]=="m")):
            number = Convert_To_Number(LayerNames[i])

            zCoord.append(number) 

        elif (LayerNames[i][0:2]=="HZ_"):                      # this is used in some old geomIO examples   
            number = Convert_To_Number(LayerNames[i][3:len(LayerNames[i])])
            zCoord.append(number) 

        else:
            try:    # try converting LayerName to float
                #number = float(LayerNames[i])
                number = Convert_To_Number(LayerNames[i])
            except: # set to 0 otherwise
                number = 0.0
            zCoord.append(number)     
            
    # Optional: sort the list with zCoord values from lowest->highest & compute indices     
    li=[]
    for i in range(len(zCoord)):
        if zCoord[i]==None:
            li.append([1e25,i])
        else:
            li.append([zCoord[i],i])
    #li.sort(reverse=True)               # sort from top->bottom to be consistent with rest of code
    li.sort()
    
    sort_ind    =   [x[1] for x in li]
    
    return zCoord, sort_ind


def Convert_To_Number(string):
    """
        This takes a string and converts it to a number

        The reason we need this is that names can be "10p25" where "p" implies a point
    """
    # plus @ beginning implies +
    if string[0] == "p":
        string = string[1:len(string)]

    string = string.replace('p','.')     # p in middle = '.'
    string = string.replace("m","-")     # minus sign

    number = float(string)
    return number


def get_CurveNames(svgFileData):
    """
    Parameters
    ----------
    Data  :   svgFileData structure 

    Returns
    -------
    UniqueCurveNames : List with unique Curve Names on "real" layers
    """

    Commented           =   svgFileData.Commented
    CurvesNames_init    =   svgFileData.CurveNames
    LayerNames          =   svgFileData.LayerNames
    CurveNames          =   []
    for i in range(len(CurvesNames_init)):
        if (Commented[i]==False) & (LayerNames[i]!="Reference"):
            CurveNames.append(CurvesNames_init[i])     
    
    CurveNamesUnique = list(set(CurveNames))
    
    return CurveNamesUnique



def get_Scaling(CurveNames, Curves, LayerNames):
    """
    Parameters
    ----------
    CurveNames  :   List with names of the curves
    Curves      :   Path info about the curves
    LayerNames  :   List with names of the layers on which the curves are
    

    Returns
    -------
    Scaling_data :  Scaling Named Tuple with scaling info: 
                    x0,         y0,       # lower left corner in real world coordinates
                    x0_SVG,     y0_SVG,   # lower left corner in SVG coordinates
                    dx,         dy        # spacing to go from SVG -> real world coordinates
                     
                    # x_real = (x_inkscape - x0_SVG)*dx + x0
                    # y_real = (y_inkscape - y0_SVG)*dy + y0

    """
    
    Scaling_data    =   Scaling(None,None,None,None,None,None)     
    for i in range(len(Curves)):
        if (LayerNames[i]=="Reference"):
            # A reference layer is present: 
            #   interpret the name of the curve to extract "real world" 
            #   x/y coordinates of begin & end-points
            CoordRef_str    =   CurveNames[i][10:len(CurveNames[i])-1]
            CoordRef_list   =   CoordRef_str.split(",")
            CoordRef        =   list(map(float, CoordRef_list))
            x0              =   CoordRef[0]
            y0              =   CoordRef[1]
            
            # Retrieve the bounding box of the reference curve
            curve           =   Curves[i]
            x0_SVG,x1_SVG,y0_SVG,y1_SVG =   curve.bbox()

            # Compute scaling 
            dx              =   (CoordRef[2]- CoordRef[0])/(x1_SVG - x0_SVG)
            dy              =   (CoordRef[3]- CoordRef[1])/(y1_SVG - y0_SVG)

            # Store scaling data in NamedTuple:
            Scaling_data    =   Scaling(x0,y0,x0_SVG,y0_SVG,dx,dy)              

    return Scaling_data
