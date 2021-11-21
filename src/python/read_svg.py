# This contains various routines to read *.SVG files, with multiple layers and multiple curves
from svgpathtools import svg2paths, real, imag, Line, svg2paths2, Document
from typing import NamedTuple
 
class svgFileData(NamedTuple):
    CurveNames: list
    Curves    : list
    LayerNames: list
    numLayers:  int
    Commented:  list
    zCoord:     list


def getLayers_General(inFile, Verbose=True):
    """
    Reads an SVG file and returns all curves, layers and names (plus whether a layer is commented or not) 
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
            if (layer_str[0]=="#") | (layer_str[0]=="$"):
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
            elif path.get('label') != None:
                curve_str =  path.get('label')
    
            # Store 
            LayerNames.append(layer_str)
            CurveNames.append(curve_str)
            Curves.append(paths[CurveNumber])
            Commented.append(comment_layer)
            CurveNumber += 1
            if Verbose:
                print("   Found PATH   : " + curve_str)
                  
                
        # in some cases, Affinity Design puts another layer around the paths, so we check that as well here
        for elem1 in elem.iterfind(SVG_GROUP_TAG, SVG_NAMESPACE):
            for path2 in elem1.iterfind(SVG_PATH_TAG, SVG_NAMESPACE):
                curve_str = None
                if elem1.get('id')!=None:
                    curve_str =  elem1.get('id')
                elif elem1.get('label') != None:
                    curve_str =  elem1.get('label')
                
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

    # Read Reference layer & determine scaling if it exists

    # Interpret zLayers
    zCoord, sort_ind = get_zCoords(LayerNames, Commented)

    # Sort other lists accordingly, so it is all from lowest->highest coordinate 
    #  (with commented & reference layers listed afterwards)
    CurveNames_sort =   [CurveNames[i]  for i in sort_ind]
    Curves_sort     =   [Curves[i]      for i in sort_ind]
    LayerNames_sort =   [LayerNames[i]  for i in sort_ind]
    Commented_sort  =   [Commented[i]   for i in sort_ind]
    zCoord_sort     =   [zCoord[i]      for i in sort_ind]

    # Store data in Named Tuple (easier to handle later)
    Data = svgFileData(CurveNames_sort, Curves_sort, LayerNames_sort, numLayers, Commented_sort, zCoord_sort)

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
    for i in range(len(LayerNames)-1):
     #   print(LayerNames[i] +"  "+ str(i))
        if (Commented[i]==True) :
            zCoord.append(None)    

        elif (LayerNames[i]=="Reference"):
            zCoord.append(None)    

        elif (LayerNames[i][0]=="p"):
            number = float(LayerNames[i][1:len(LayerNames[i])])
            zCoord.append(number) 

        elif (LayerNames[i][0]=="m"):
            number = -float(LayerNames[i][1:len(LayerNames[i])])
            zCoord.append(number) 
            
        else:
            try:    # try converting LayerName to float
                number = float(LayerNames[i])
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
    li.sort()

    sort_ind    =   [x[1] for x in li]
    
    return zCoord, sort_ind


def get_CurveNames(Data):
    """
    Parameters
    ----------
    Data  :   svgFileData structure 

    Returns
    -------
    UniqueCurveNames : List with unique Curve Names on "real" layers
    """

    Commented   =   Data[4]  
    Curves      =   Data[0]
    LayerNames  =   Data[2]
    CurveNames  =   []
    for i in range(len(Curves)):
        if (Commented[i]==False) & (LayerNames[i]!="Reference"):
            CurveNames.append(Curves[i])     
    
    CurveNamesUnique = set(CurveNames)
    
    return CurveNamesUnique