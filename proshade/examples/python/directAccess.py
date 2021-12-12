######################################################
######################################################
#   \file directAccess.py
#   \brief This code demonstrates the usage of the ProSHADE tool in the advanced mode.
#
#   This file should be the main source of wisdom when it comes to using ProSHADE in
#   Python in the advanced access mode. It demonstrates how many of the possible
#   tasks can be done, albeit one would not expect to be doing them all in one run...
#
#   Therefore, this should serve more as a "cook-book" rather than executable file. Generally,
#   the procedures shown here include creating the settings and structure objects, both from
#   file and from already existing Python array. It also shows how the map can be processed,
#   mapped onto spheres and how spherical harmonics are obtained, including how they can be
#   directly accessed from Python.
#
#   Moreover, structure re-boxing is shown, as well as the distances computation. The file also
#   demonstrates how the rotation function can be computed and how the E matrices, SO(3) coeff-
#   icients and the self-rotation map can be accessed, followed by a demonstration of how the
#   symmetry detection is called and results read.
#
#   Finally, the file shows and explains how the map overlay can be completed from Python, including
#   computation of the rotation function, access to it programatically and the same for the translation
#   function. Direct access to the optimal rotation Euler angles and rotation matrix is shown, as well
#   as direct access to the optimal translation vector.
#
#   Copyright by Michal Tykac and individual contributors. All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#   1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#   2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#   3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
#   This software is provided by the copyright holders and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In     no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data     or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility     of such damage.
#
#   \author    Michal Tykac
#   \author    Garib N. Murshudov
#   \version   0.7.6.2
#   \date      DEC 2021
######################################################
######################################################

######################################################
### Import modules
### ==============
###
### This is where Python modules are loaded.
###

### System modules
import sys
import numpy

### Import ProSHADE from non-system folder (local installation assumed)
sys.path.append                                       ( "/Users/mysak/BioCEV/proshade/experimental/install/pythonModule" )
import proshade

######################################################
### Create the ProSHADE_settings object
### ===================================
###
### This object contains all the settings that
### will be used throughout the ProSHADE exection.
### Therefore, the user should NOT change any
### of these settings once structures are read,
### as these values are assumend not the be
### changing throughout the run, except by the
### functions themselves.
###
### If you are running a different structure for
### a different purpose, please create a new
### settings object.
###

### Create the object
pSet                                                  = proshade.ProSHADE_settings ()
        
### Set settings values
pSet.task                                             = proshade.Distances
pSet.verbose                                          = 1
pSet.rotationUncertainty                              = 5.0
pSet.moveToCOM                                        = True
pSet.setMapResolutionChange                           ( False )

######################################################
### Create ProSHADE_structure object
### ================================
###
### This object is the basis of the advanced
### access interface. Most of the processing
### and computation functions are callable from
### this object. Each structure (map or co-ords)
### will need its own object. Also, these are
### pointers in C++ and so not deleting them
### before Python code termination will cause
### memory error message.
###
pStruct                                               = proshade.ProSHADE_data ( )

######################################################
### Read in structure
### =================
###
### This function reads in a molecular structure
### from a file. It internally uses the Gemmi library
### for reading in maps and co-ordinates, supporintg
### gunzipped files as well. It takes three arguments:
###
### 1) String: The path and filename of the structure
###            to be read.
### 2) Int: The order of the structure. This is used
###         to distinguish multiple file outputs.
###         Please used different value for each
###         structure.
### 3) ProSHADE_settings*: The settings object with
###                        all values required for
###                        the reading in.
###
### If the function fails, it will exit proshade and
### print error message explaining what happened.
###
### At this point, the following variables are
### meaningfully (in Python access meaning) filled in:
###
### string   pStruct.fileName                         < This is the original file from which the data were obtained.
### float    pStruct.xDimSize                         < This is the size of the map cell x dimension in Angstroms.
### float    pStruct.yDimSize                         < This is the size of the map cell y dimension in Angstroms.
### float    pStruct.zDimSize                         < This is the size of the map cell z dimension in Angstroms.
### float    pStruct.aAngle                           < This is the angle a of the map cell in degrees.
### float    pStruct.bAngle                           < This is the angle b of the map cell in degrees.
### float    pStruct.cAngle                           < This is the angle c of the map cell in degrees.
### int      pStruct.xDimIndices                      < This is the size of the map cell x dimension in indices.
### int      pStruct.yDimIndices                      < This is the size of the map cell y dimension in indices.
### int      pStruct.zDimIndices                      < This is the size of the map cell z dimension in indices.
### int      pStruct.xGridIndices                     < As far as I know, this is identical to the xDimIndices.
### int      pStruct.yGridIndices                     < As far as I know, this is identical to the yDimIndices.
### int      pStruct.zGridIndices                     < As far as I know, this is identical to the zDimIndices.
### int      pStruct.xAxisOrder                       < This is the order of the x axis.
### int      pStruct.yAxisOrder                       < This is the order of the y axis.
### int      pStruct.zAxisOrder                       < This is the order of the z axis.
### int      pStruct.xAxisOrigin                      < This is the origin position along the x axis.
### int      pStruct.yAxisOrigin                      < This is the origin position along the y axis.
### int      pStruct.zAxisOrigin                      < This is the origin position along the z axis.
### int      pStruct.xFrom                            < This is the starting index along the x axis.
### int      pStruct.yFrom                            < This is the starting index along the y axis.
### int      pStruct.zFrom                            < This is the starting index along the z axis.
### int      pStruct.xTo                              < This is the final index along the x axis.
### int      pStruct.yTo                              < This is the final index along the y axis.
### int      pStruct.zTo                              < This is the final index along the z axis.
###
### Please note that if you change any of these,
### then things may stop making any sense. So,
### only change things if you know enough about
### the underlying code or if you are prepared
### to experiment and possibly get crazy results.
###
pStruct.readInStructure                               ( "/Users/mysak/LMB/proshade/exp/demo/C3.pdb", 0, pSet )

######################################################
### Create ProSHADE_structure object from map
### =========================================
###
### An alternative approach to creating the
### ProSHADE_data object is not to use a structure
### saved on the drive, but instead to supply all
### the required information as well as the map.
###
### This allows the user to obtain the map from
### any source they like, but it requires them
### to know the map information and how it should
### be supplied. In this case, the ProSHADE_data
### constructure takes the following arguments:
###
### ProSHADE_settings   - object                      < This object contains all the settings for further processing.
### structureName       - string                      < A string to be used in naming any outout files from this structure.
### inputMap            - 1D or 3D float array        < Array containing the map values.
### xDimAngs            - float                       < The size of x dimension in Angstroms.
### yDimAngs            - float                       < The size of y dimension in Angstroms.
### zDimAngs            - float                       < The size of z dimension in Angstroms.
### xDimInds            - int                         < The size of x dimension in terms of number of indices.
### yDimInds            - int                         < The size of y dimension in terms of number of indices.
### zDimInds            - int                         < The size of z dimension in terms of number of indices.
### xFrom               - int                         < The initial index position along x axis.
### yFrom               - int                         < The initial index position along y axis.
### zFrom               - int                         < The initial index position along z axis.
### xTo                 - int                         < The last index position along x axis.
### yTo                 - int                         < The last index position along y axis.
### zTo                 - int                         < The last index position along z axis.
### ord                 - int                         < The order of the struct object in ProSHADE processing - important for multiple objects processing outputs.
###
### NOTE: There are two main conditions that need
### to be fullfilled for the constructor call to
### work. 1) The map dimensions needs to be the
### same as the x/y/zDimInds variables and 2)
### x/y/zTo - x/y/zFrom + 1 = x/y/zDimInds
###
### NOTE2: This function makes a lot of assumptions
### (all angles are 90 degrees, axis grids are
### equal to indices, axis order is XYZ and axis
### origin is the first index in all dimensions).
### If any of these are not true, the user is required
### to change the appropriate internal values after
### this function has returned the object.
###

######################################################
### Release the previous object
del pStruct

######################################################
### Set example values
xDimIndices                                           = 100
yDimIndices                                           = 120
zDimIndices                                           = 60
xDimAngstroms                                         = xDimIndices * 1.3
yDimAngstroms                                         = yDimIndices * 1.3
zDimAngstroms                                         = zDimIndices * 1.3
xFrom                                                 = int ( -xDimIndices/2 )
yFrom                                                 = int ( -yDimIndices/2 )
zFrom                                                 = int ( -zDimIndices/2 )
xTo                                                   = int ( (xDimIndices/2) )
yTo                                                   = int ( (yDimIndices/2) )
zTo                                                   = int ( (zDimIndices/2) )
ord                                                   = 0

if xDimIndices % 2 == 0:
    xTo                                               = xTo - 1
    
if yDimIndices % 2 == 0:
    yTo                                               = yTo - 1
    
if zDimIndices % 2 == 0:
    zTo                                               = zTo - 1

######################################################
### Create example map (this will be a ball in the middle of the map)
testMap = numpy.empty ( [ ( xDimIndices * yDimIndices * zDimIndices ) ] )
for xIt in range( 0, xDimIndices ):
    for yIt in range( 0, yDimIndices ):
        for zIt in range( 0, zDimIndices ):
            ind                                       = zIt + zDimIndices * ( yIt + yDimIndices * xIt )
            testMap[ind]                              = 1.0 / ( numpy.sqrt( numpy.power ( (xDimIndices/2) - xIt, 2.0 ) +
                                                                            numpy.power ( (yDimIndices/2) - yIt, 2.0 ) +
                                                                            numpy.power ( (zDimIndices/2) - zIt, 2.0 ) ) + 0.01 )

######################################################
### Create the ProSHADE_data object without structure file on drive
pStruct                                               = proshade.ProSHADE_data ( "python_map_test",
                                                                                 testMap,
                                                                                 xDimAngstroms,
                                                                                 yDimAngstroms,
                                                                                 zDimAngstroms,
                                                                                 xDimIndices,
                                                                                 yDimIndices,
                                                                                 zDimIndices,
                                                                                 xFrom,
                                                                                 yFrom,
                                                                                 zFrom,
                                                                                 xTo,
                                                                                 yTo,
                                                                                 zTo,
                                                                                 ord )

######################################################
### Should we ever need to use 3D map instead of 1D
### map, the same constructor is capable of dealing
### with this case as well.
testMap3D = numpy.empty ( ( xDimIndices, yDimIndices, zDimIndices ) )
for xIt in range( 0, xDimIndices ):
    for yIt in range( 0, yDimIndices ):
        for zIt in range( 0, zDimIndices ):
            testMap3D[xIt][yIt][zIt]                  = 1.0 / ( numpy.sqrt( numpy.power ( (xDimIndices/2) - xIt, 2.0 ) +
                                                                            numpy.power ( (yDimIndices/2) - yIt, 2.0 ) +
                                                                            numpy.power ( (zDimIndices/2) - zIt, 2.0 ) ) + 0.01 )

pStruct2                                              = proshade.ProSHADE_data ( "python_map_test",
                                                                                 testMap3D,
                                                                                 xDimAngstroms,
                                                                                 yDimAngstroms,
                                                                                 zDimAngstroms,
                                                                                 xDimIndices,
                                                                                 yDimIndices,
                                                                                 zDimIndices,
                                                                                 xFrom,
                                                                                 yFrom,
                                                                                 zFrom,
                                                                                 xTo,
                                                                                 yTo,
                                                                                 zTo,
                                                                                 ord )

######################################################
### Release the unnecessary object
del pStruct2

######################################################
### Write the internal map to disk
### ==============================
###
### This function writes the current internal
### map in CCP4 MAP format into a file given as
### its only argument. It uses all the internal
### values in the ProSHADE_data structure, so
### the user is responsible for these not being
### changed/or still making sense.
###
pStruct.writeMap                                      ( "initialMap.map" )

######################################################
### Get internal map representation
### ===============================
###
### ProSHADE allows access to the internal map that
### it uses for all its computations. This allows the
### user to check any particular ProSHADE's functionality
### effect or simply running only the part of ProSHADE
### the user likes and then retrieving the map for
### other processing.
###
### The maps are outputted in the numpy.ndarray format.
###
### This function can be called at any time on
### the ProSHADE_data object to get the current
### internal map representation.
###
initialMap                                            = pStruct.getMap ( )

print                                                 ( "Mean map value: " + str( numpy.mean ( initialMap ) ) )
# Expected output: Mean map value: 0.025628592352769597

######################################################
### Process internal map
### ====================
###
### This function is where all the initial map
### manipulation happens. If requested (i.e. set
### in the settings object), it can do the following
### modifications of the internal map:
###
### 1) Map invertion: This switches all XYZ positions
###    to -X -Y -Z positions.
### 2) Map normalisation: This changes all density to
###    have mean zero and standard deviation one.
### 3) Map masking: Here, the map is blurred by a factor
###    and then a threshold is computed from the blurred
###    map. All passing points are left, non-passing points
###    become zeroes.
### 4) Map centering: The map will be moved to have its
###    centre of mass at the co-ordinate a/2, b/2, c/2.
### 5) Add extra space: This will add specified number of
###    Angstroms before and after the data along all
###    dimensions. This is useful to avoid unwanted inter-
###    actions from periodic cells.
### 6) Removing phase information: If you want molecular
###    replacement type of search instead, this option is
###    available.
###
### NOTE: All of these modifications can be done on
### co-ordinates originating internal maps as well.
pStruct.processInternalMap                            ( pSet )

######################################################
### Map re-boxing
### =============
###
### This is the first task to be discussed as
### such. In the re-boxing functionality, the
### user first needs to determine the boundaries
### from which the new, re-boxed map should be
### created and then he needs to create a new
### empty structure, finally having it filled
### from the original structure using the already
### defined boundaries. The following three
### steps will accomplish just that using the
### ProSHADE map masking approach, but the users
### are free to change the bounds or determine their
### own boundaries, if they so please.
###

######################################################
### Determine new boundaries from mask (or custom)
### ==============================================
###
### If the user requires ProSHADE to determine
### new boundaries based on the ProSHADE masking
### procedure (and using the settings object
### values), the following function does just that.
###
### NOTE: If the users want to supply their own
### boundaries, this step can be skipped. Just
### make sure your custom boundaries are in the
### numpy.ndarray format, have length of 6 and
### dtype = int64. Also, the 6 numbers have
### meaning as follows:
###
###    [0] = min X-axis index
###    [1] = max X-axis index
###    [2] = min Y-axis index
###    [3] = max Y-axis index
###    [4] = min Z-axis index
###    [5] = max Z-axis index
###
minimalBounds                                         = pStruct.getReBoxBoundaries ( pSet )

print                                                 ( minimalBounds )
# Expected output: [  6 109   6 129   6  69]

######################################################
### Create new structure to hold the new map
### ========================================
###
### Create a new structure, which will have the
### re-boxed map and values. This is so that the
### re-boxing would not be done in place.
###
reBoxStr                                              = proshade.ProSHADE_data ( )

######################################################
### Set the re-boxed structure values and map
### =========================================
###
### Fill the new structure with the calling
### structure's map values in the boudaries
### supplied and set all of its required fields
### as well. The syntax of the call is that the
### calling structure is the one which is the
### source of map, while the second argument
### structure is the empty one.
###
pStruct.createNewMapFromBounds                        ( minimalBounds, reBoxStr, pSet )

######################################################
### Map internal map to spheres
### ===========================
###
### This function does the automatic spherical
### harmonics settings determination (unless
### these are already set in the settings object)
### and then it creates the required number of
### concentric spheres, finally mapping the inter-
### nal map onto the spheres using tri-linear
### interpolation.
###
### This will fill the following variables properly
###
### int           pStruct.noSpheres                   < The number of spheres with map projected onto them.
### list          pStruct.spherePos                   < List of sphere radii from the centre of the map.
###
###
pStruct.mapToSpheres                                  ( pSet )

######################################################
### Compute spherical harmonics
### ===========================
###
### This function takes the shells with the mapped
### data and proceeds to compute the spherical
### harmonics for each of them. This may take some
### time depending on the bandwidth and number of
### shells.
###
pStruct.computeSphericalHarmonics                     ( pSet )

######################################################
### Accessing spherical harmonics
### =============================
###
### In order to access the spherical harmonics
### values for each sphere, please use the
### following functions:
###
### 1) getSphericalHarmonics()
###     This function returns a 2D numpy.ndarray of
###     complex numbers. This 2D array is organised so
###     that the first dimension is the sphere number
###     and the second dimension is a 1D array of all
###     the spherical harmonics values. This 1D array
###     is complicated by the fact that each sphere
###     can have different number of bands (depending
###     on the settings.progressiveSphereMapping value).
###     This means that some values for some spheres
###     are empty. Best way of accessing this array is
###     to use the second function.
###
### 2) findSHIndex()
###     This function takes the shell, band and order
###     values and returns the correct index for these
###     in the 1D array. Please note that the order
###     value needs to be in range -band <= order <= band
###     and NOT in the 0 <= order <= (2 * band) + 1.
###
sphericalHarmonics                                    = pStruct.getSphericalHarmonics ( )

shell =  3
band  =  4
order = -2
Shell3Band4OrderMin2Value                             = sphericalHarmonics[shell][pStruct.findSHIndex(shell, band, order)]
print                                                 ( Shell3Band4OrderMin2Value )
# Expected output: (6.14462905753556e-05+0.0004590295017376565j)

######################################################
### Computing distances between two structures
### ==========================================
###
### In order to compute shape distances between
### two structures, two structures need to exist
### :-). All of the structures for which the distances
### are to be computed need to have their spherical
### harmonics computed.
###
### Next, the pairwise distances can simply be obtained
### by calling the appropriate distances computation
### functions as shown below.
###

### Create a second structure to have someting to compute distances to
pStruct2                                              = proshade.ProSHADE_data ( )
pStruct2.readInStructure                              ( "/Users/mysak/BioCEV/proshade/playground/emd_6324.map", 1, pSet )
pStruct2.processInternalMap                           ( pSet )
pStruct2.mapToSpheres                                 ( pSet )
pStruct2.computeSphericalHarmonics                    ( pSet )

### Get the three descriptors
energyLevelsDescriptor                                = proshade.computeEnergyLevelsDescriptor     ( pStruct, pStruct2, pSet )
traceSigmaDescriptor                                  = proshade.computeTraceSigmaDescriptor       ( pStruct, pStruct2, pSet )
fullRotationFunctionDescriptor                        = proshade.computeRotationFunctionDescriptor ( pStruct, pStruct2, pSet )

print                                                 ( energyLevelsDescriptor )
# Expected output: 0.08404706339177564
print                                                 ( traceSigmaDescriptor )
# Expected output: 0.2703318898591193
print                                                 ( fullRotationFunctionDescriptor )
# Expected output: 0.24178888432250034

######################################################
### Delete the C++ pointers
### =======================
###
del pStruct
del pStruct2
del pSet

######################################################
### Create new structure for symmetry detection
### ===========================================
###
### As the tasks required different settings, a new
### structure will be created, this time with
### settings optimised for symmetry detection.
###
### To get settings optimal for symmetry detection,
### a new ProSHADE_settings object will be created,
### this time with the symmetry task as its constructor
### argument. This will automatically set the default
### values towards the task given.
###
pSet                                                  = proshade.ProSHADE_settings ( proshade.Symmetry )
pSet.requestedResolution                              = 10.0  ## This is to make the example run faster as here we are not really interested in accuracy of results
pSet.changeMapResolution                              = True  ## but rather in showing how to use ProSHADE.

pStruct                                               = proshade.ProSHADE_data ( )
pStruct.readInStructure                               ( "/Users/mysak/BioCEV/proshade/playground/emd_0116.map.gz", 1, pSet )
pStruct.processInternalMap                            ( pSet )
pStruct.mapToSpheres                                  ( pSet )
pStruct.computeSphericalHarmonics                     ( pSet )

######################################################
### Computing self-rotation function
### ================================
###
### This function computes the self-rotation
### function by firstly computing and normalising
### the E matrices, then combining these into
### the SO(3) coefficients and finally calcula-
### ting the inverse SO(3) Fourier Transform
### (SOFT) from them. It therefore needs to be
### called before any symmetry detection can be
### attempted.
###
pStruct.computeRotationFunction                       ( pSet )

######################################################
### Accessing E Matrices
### ====================
###
### ProSHADE allows access to the E matrices
### ( Integral _0 ^rMAX ( c1^lm * c2*^lm' ) of
### the structure combination. These are
### returned as a 3D Numpy array with indices being
### band of the E matrix, order1 of the E matrix
### and order2 of the E matrix. Note that because
### indices need to go from zero, the order of indexing
### goes like, for example, this:
###
### BAND = 2 || ORDER = -2  || INDEX [2][0]
### BAND = 2 || ORDER = -1  || INDEX [2][1]
### BAND = 2 || ORDER =  0  || INDEX [2][2]
### BAND = 2 || ORDER =  1  || INDEX [2][3]
### BAND = 2 || ORDER =  2  || INDEX [2][4]
### i.e. order index = ORDER + BAND
###
### NOTE: As Numpy arrays have single shape, the
### lower bands (which will have less orders)
### are padded with zeroes to have the same length
### as the largest band.
###
eMat                                                  = pStruct.getEMatrix ( )


Band4OrderOneMin2OrderTwo3EMatrixValue                = eMat[4][2][7] # Band = 4, Order1 = -2 and Order2 = 3

print ( Band4OrderOneMin2OrderTwo3EMatrixValue )
# Expected output: (2.0623081672724683e-11-2.388437874366713e-12j)

######################################################
### Accessing SO(3) coefficients
### ============================
###
### ProSHADE also allows access to the SO(3)
### coefficients computed by normalising the E
### matrix values and dealing with the signs.
### The inverse SO(3) Fourier Transform (SOFT)
### of these values then results in the rotation
### function. The complete complex array of
### these values can be accessed as shown,
### however, the organisation of the array is
### done by the SOFT library and reflects
### internal value symmetries. Therefore, to
### access a specific value, please use the
### so3CoeffsArrayIndex() function as shown.
###
### Given that the SO(3) coefficients have the same
### indexing as the E matrices, ProSHADE provides them
### to the user in the same format as the E matrices.
### This does, however, mean that the same caveats do
### apply to the SO(3) coefficients array as did to the
### E matrix array.
###
so3Coeffs                                             = pStruct.getSO3Coefficients ( )


Band4OrderOneMin2OrderTwo3SO3CoeffsValue              = so3Coeffs[4][2][7] # Band = 4, Order1 = -2 and Order2 = 3

print ( Band4OrderOneMin2OrderTwo3SO3CoeffsValue )
# Expected output: (-6.108395846399666e-11+7.07436658725007e-12j)

######################################################
### Accessing self-rotation function
### ================================
###
### ProSHADE also gives access to the self -
### rotation function as shown next. The returned 3D complex
### array has dimensions 2 * bandwidth. Please note that
### the indices have nothing to do with the angle values,
### if you want to know the Euler angle values for a
### particular index, you need to convert it
### yourself. Alternatively (and a recommended
### approach is), you can use the
### getRotationMatrixFromEulerIndices() function
### also demonstrated here.
###
selfRotationFunction                                  = pStruct.getRotationFunctionMap ( )

### Find the map maximum
rotMapMax                                             = numpy.where ( selfRotationFunction == numpy.amax ( selfRotationFunction ) )

### Find maximum value
print                                                 ( "Rotation map maximum is: " + str( selfRotationFunction[rotMapMax[0][0]][rotMapMax[1][0]][rotMapMax[2][0]] ) )

### Expected output: Rotation map maximum is: (0.9960890372115769-4.244188357793958e-18j)

### Find rotation matrix for the maximum
rotMatMaxVal                                          = pStruct.getRotationMatrixFromSOFTCoordinates ( rotMapMax[0][0], rotMapMax[1][0], rotMapMax[2][0] )
print                                                 ( rotMatMaxVal )

# Expected output: [[ 1.00000000e+00 -2.46330734e-16 -0.00000000e+00]
# Expected output:  [ 2.46330734e-16  1.00000000e+00  0.00000000e+00]
# Expected output:  [ 0.00000000e+00 -0.00000000e+00  1.00000000e+00]]

######################################################
### Run symmetry detection
### ======================
###
### Once the self-rotation computation is done
### (does not need getting it into python), the
### symmetry detection algorithm can be run as
### shown. Once the detectSymmetryInStructurePython()
### function is complete, the detected symmetry
### values can be obtained as demonstrated.
###

### Detect symmetry
pSet.requestedSymmetryFold                            = 3
pSet.requestedSymmetryType                            = "D"
pStruct.detectSymmetryInStructure                     ( pSet )

### Retrieve results
recSymmetryType                                       = pStruct.getRecommendedSymmetryType ( pSet )
recSymmetryFold                                       = pStruct.getRecommendedSymmetryFold ( pSet )
recSymmetryAxes                                       = pStruct.getRecommendedSymmetryAxes ( pSet )

### Print results
print                                                 ( "Detected " + str( recSymmetryType ) + "-" + str( recSymmetryFold ) + " symetry." )

# Expected output: Detected D-3 symetry.

### Print more results
print                                                 ( "Fold      x         y         z       Angle     Height    Average FSC" )
for iter in range ( 0, len( recSymmetryAxes ) ):
     print                                            ( "  %s    %+1.3f    %+1.3f    %+1.3f    %+1.3f    %+1.4f      %+1.4f" % ( recSymmetryAxes[iter][0], recSymmetryAxes[iter][1], recSymmetryAxes[iter][2], recSymmetryAxes[iter][3], recSymmetryAxes[iter][4], recSymmetryAxes[iter][5], recSymmetryAxes[iter][6] ) )
     
# Expected output: Fold      x         y         z       Angle     Height    Average FSC
# Expected output:   3.0    +0.000    +0.000    +1.000    +2.094    +0.9726      +0.9206
# Expected output:   2.0    +1.000    +0.000    +0.000    +3.142    +0.9920      +1.0000

######################################################
### Get more symmetry results
### =========================
###
### Once the symmetry detection has been run,
### it is now possible to access more detailed
### results. These include the list of all
### detected cyclic point groups, a list of indices
### of these cyclic point groups which form any
### particular non-cyclic point group as well as
### list of all group elements for any point group
### comprised from detected cyclic point groups. For
### more details, please see the advancedAccess_symmetry.py
### file.
###

### Get list of all cyclic point groups
allCAxes                                              = pStruct.getAllCSyms ( pSet )

print                                                 ( "Found a total of " + str( len ( allCAxes ) ) + " cyclic point groups." )
# Expected output: Found a total of 9 cyclic point groups.

### Get a list of all non-cyclic point groups detected and indices of the cyclic groups forming them
allNonCAxesIndices                                    = pStruct.getNonCSymmetryAxesIndices ( pSet )

print                                                 ( "Found a total of " + str( len ( allNonCAxesIndices["D"] ) ) + " dihedral point groups." )
# Expected output: Found a total of 21 dihedral point groups.

### Get group elements for the best dihedral symmetry (indices 2 and 26)
bestDCombination                                      = []
bestDCombination.append                               ( allNonCAxesIndices["D"][2][0] )
bestDCombination.append                               ( allNonCAxesIndices["D"][2][1] )
allGroupElements                                      = pStruct.getAllGroupElements ( pSet, bestDCombination, "D" )

print                                                 ( "Found a total of " + str( len ( allGroupElements ) ) + " elements." )
# Expected output: Found a total of 12 elements.

print                                                 ( allGroupElements[1] )
# Expected output: [[ 0.5        0.8660254  0.       ]
# Expected output:  [-0.8660254  0.5        0.       ]
# Expected output:  [ 0.         0.         1.       ]]

######################################################
### Computing and combining group elements
### ======================================
###
### Assuming you have modified / created your own
### cyclic groups and want to obtain their group
### elements, ProSHADE can compute these as follows:
gr1Elements                                           = proshade.computeGroupElementsForGroup ( 1,      # X-axis element of the group rotation axis
                                                                                                0,      # Y-axis element of the group rotation axis
                                                                                                0,      # Z-axis element of the group rotation axis
                                                                                                4 )     # Fold
                                                      
gr2Elements                                           = proshade.computeGroupElementsForGroup ( 0,      # X-axis element of the group rotation axis
                                                                                                1,      # Y-axis element of the group rotation axis
                                                                                                0,      # Z-axis element of the group rotation axis
                                                                                                2 )     # Fold
                                                  
### Moreover, if you have computed group elements for
### multiple groups and would like to combine these
### into a single group, possibly discovering new
### group elements, this can also be done by ProSHADE
### by using the following function:
combinedGroupElements                                 = proshade.joinElementsFromDifferentGroups ( gr1Elements, # Elements of the first group
                                                                                                   gr2Elements, # Elements of the second group
                                                                                                   0.0001,      # Matrix tolerance (i.e. if the abs ( trace ( Mat1 * Mat2^(T) ) - 3.0 ) must be below this number for two matrices to be considered identical )
                                                                                                   True )       # Should new group elements resulting from multiplication of elements from group 1 with elements from group 2 be computed?

print ( len ( combinedGroupElements ) )
# Expected output: 8

##############################################
### Delete the C++ pointers
### =======================
###
del reBoxStr
del pStruct
del pSet

##############################################
### Computing the map overlay
### =========================
###
### As the overlay code does require the phase
### to be removed for optimal rotation compu-
### tation, most of the already demonstrated
### functions will need to be called again.
### However, they will not be described in much
### detail, as the descriptions are above.
###
### More specifically to the overlay computation,
### the user here should know how ProSHADE does
### this in order to call the functions properly.
### If you do not want to know the details, I
### suggest using the simpleAccess files instead,
### as they are almost as fast and do not require
### any internal knowledge (albeit they do not have
### as much flexibility as the advancedAccess)
###
### Now, the overlay mode is divided into two
### separate steps. Firstly, the phase is removed
### from the two internal maps and the resulting
### Patterson maps are subjected to the spherical
### harmonics decomposition. This makes sure the
### centering is precise, while it still allows
### for computing the rotation function (from the
### spherical harmonics values and inverse SOFT
### transform). The highest peak of the rotation
### function then gives the global optimal rotation
### angles.
###
### Then, all the internal data are released and the
### two structures are read again, this time with
### phase. The moving structure then has the optimal
### rotation applied (it must be retained from the
### phase-less step); now, the translation function
### can be computed for two optimally rotated structures.
### The results from the translation function then
### form the optimal translation vector. To do this,
### please follow the next section.
###

##############################################
######################################################
### Create the settings object
pSet                                                  = proshade.ProSHADE_settings ( )

######################################################
### Set basic settings values
pSet.task                                             = proshade.OverlayMap
pSet.verbose                                          = 4
pSet.requestedResolution                              = 4.0
pSet.usePhase                                         = False
pSet.changeMapResolution                              = True
pSet.maskMap                                          = False
pSet.moveToCOM                                        = False
pSet.normaliseMap                                     = False
pSet.reBoxMap                                         = False

######################################################
### Create data objects
pStruct_static                                        = proshade.ProSHADE_data ( )
pStruct_moving                                        = proshade.ProSHADE_data ( )

######################################################
### Read in the structures
pStruct_static.readInStructure                        ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/bf/1BFO_A_dom_1.pdb", 0, pSet ) # This is a BALBES domain 1BFO_A_dom_1.
pStruct_moving.readInStructure                        ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/h8/1H8N_A_dom_1.pdb", 1, pSet ) # This is a BALBES domain 1H8N_A_dom_1.

######################################################
### Get spherical harmonics for both structures
pStruct_static.processInternalMap                     ( pSet )
pStruct_moving.processInternalMap                     ( pSet )
        
pStruct_static.mapToSpheres                           ( pSet )
pStruct_moving.mapToSpheres                           ( pSet )
        
pStruct_static.computeSphericalHarmonics              ( pSet )
pStruct_moving.computeSphericalHarmonics              ( pSet )

######################################################
### Create structure objects (phase-less)
### =====================================
###
### This is the first step not already described
### above. Albeit this step is similar to the
### symmetry detection self-rotation function,
### here the rotation function is computed by
### combining the spherical harmonics coefficients
### from two different structures rather than from
### the same structure. The combination then results
### in the SO(3) group coefficients, which can be
### converted to the rotation function by the
### Fourier transform on the SO(3) group.
###
pStruct_moving.getOverlayRotationFunction             ( pSet, pStruct_static )

######################################################
### The E matrix and the SO(3) coefficients are
### still accessible by the same functions.
eMat                                                  = pStruct_moving.getEMatrix ( )
Band4OrderOneMin2OrderTwo3EMatrixValue                = eMat[4][2][7] # Band = 4, Order1 = -2 and Order2 = 3

print ( Band4OrderOneMin2OrderTwo3EMatrixValue )
# Expected output: (-0.0007796590626317366-0.0010195435475450075j)

so3Coeffs                                             = pStruct_moving.getSO3Coefficients ( )
Band4OrderOneMin2OrderTwo3SO3CoeffsValue              = so3Coeffs[4][2][7] # Band = 4, Order1 = -2 and Order2 = 3

print ( Band4OrderOneMin2OrderTwo3SO3CoeffsValue )
# Expected output: (0.002309289297964725+0.0030198084213981114j)

######################################################
### Acceasing rotation function and etc.
### ====================================
###
### Similarly to the self-rotation function used
### by the symmetry detection part above, ProSHADE
### allows access to the rotation function (as well
### as the E matrices and the SO(3) coefficients)
### using the same function as above. The only
### caveat that the user needs to be aware of is
### that these values are always saved in the
### moving structure class and never in the static
### structure class (static can be compared to
### multiple moving without the need for new
### static class).
###
### Get rotation function
rotationFunction                                      = pStruct_moving.getRotationFunctionMap ( )

### Find the map maximum
rotMapMax                                             = numpy.where ( rotationFunction == numpy.amax ( rotationFunction ) )

### Find maximum value
print                                                 ( "Rotation map maximum is: " + str( rotationFunction[rotMapMax[0][0]][rotMapMax[1][0]][rotMapMax[2][0]] ) )

### Expected output: Rotation map maximum is: (0.9774796595100823-3.358472249513311e-18j)

### Find rotation matrix for the maximum
rotMatMaxVal                                          = pStruct_moving.getRotationMatrixFromSOFTCoordinates ( rotMapMax[0][0], rotMapMax[1][0], rotMapMax[2][0] )
print                                                 ( rotMatMaxVal )

# Expected output: [[-0.9033352   0.11520062 -0.41317591]
# Expected output:  [ 0.11520062 -0.86270924 -0.49240388]
# Expected output:  [-0.41317591 -0.49240388  0.76604444]]

######################################################
### Finding optimal rotation
### ========================
###
### Finally, once the rotation function is available,
### the optimal rotation Euler angles (ZXZ form) can
### be obtained as shown. Technically, this is done
### by finding the highest peak in the rotation map
### and then finding the Euler angles from its co-
### ordinates, but I assume this will be faster in
### C++ rather than in Python.
###
optimalRotationAngles                                 = pStruct_moving.getBestRotationMapPeaksEulerAngles ( pSet )
optimalRotationMatrix                                 = pStruct_moving.getBestRotationMapPeaksRotationMatrix ( pSet )

print ( optimalRotationMatrix )
# Expected output: [[-0.90327969  0.11682726 -0.41284041]
# Expected output:  [ 0.11355349 -0.86280879 -0.49261201]
# Expected output:  [-0.41375284 -0.49184589  0.76609151]]

######################################################
### Delete the phase-less data
### ==========================
###
### As described above for the Overlay procedure,
### this is what needs to be done.
###
del pStruct_static
del pStruct_moving

######################################################
### Changing the settings
### =====================
###
### As described above for the Overlay procedure,
### this is what needs to be done - otherwise the
### translation cannot be obtained in the real
### space. Also, the map changing is required.
###
pSet.usePhase                                         = True
pSet.changeMapResolution                              = True

##############################################
### Read in the structures with phase
### =================================
###
### As described above for the Overlay procedure,
### this is what needs to be done. However, this
### time, the spherical harmonics decomposition
### is required only for the moving structure.
###

pStruct_static                                        = proshade.ProSHADE_data ( )
pStruct_moving                                        = proshade.ProSHADE_data ( )

######################################################
### Read in the structures
pStruct_static.readInStructure                        ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/bf/1BFO_A_dom_1.pdb", 0, pSet ) # This is a BALBES domain 1BFO_A_dom_1.
pStruct_moving.readInStructure                        ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/h8/1H8N_A_dom_1.pdb", 1, pSet ) # This is a BALBES domain 1H8N_A_dom_1.

######################################################
### Get spherical harmonics for moving structure only
pStruct_static.processInternalMap                     ( pSet )
        
pStruct_moving.processInternalMap                     ( pSet )
pStruct_moving.mapToSpheres                           ( pSet )
pStruct_moving.computeSphericalHarmonics              ( pSet )

######################################################
### Rotate map
### ==========
###
### This function does the map rotation. In ProSHADE,
### this is done by computing the spherical harmonics,
### computing the Wigner D matrices for the required
### rotation (in Euler ZXZ angles) and then multi-
### plying the coefficients. This results in spherical
### harmonics coefficients of a rotated structure,
### which can subsequently be converted to the structure
### itself by inverse spherical harmonics decomposition
### calculation. This approach, however, uses interpolation
### to get the Cartesian map positions, so the resulting
### map tends to have artefacts; therefore, the user is
### discouraged from using such maps directly. They are
### good enough to get the translation map in the next
### steps, but they are not good enough for further
### processing by ProSHADE or any other software. The
### recommended approach to computing the rotated map
### from the optimal angles (or rotation matrix) is to
### use EMDA.
###
### Nonetheless, if the user so desires, the rotated map
### can be obtained and processed as any other map and as
### described above.
###
pStruct_moving.rotateMapRealSpaceInPlace              ( optimalRotationAngles[0], optimalRotationAngles[1], optimalRotationAngles[2] )

######################################################
### Zero padd the maps
### ==================
###
### The following code adds zeroes to both maps so that
### they will have the same dimensions (which will,
### in turn, be the maximum dimensions of the two
### structures). This is required for the Fourier
### coefficients used to compute the translation
### map to be of the same orders.
###
pStruct_static.zeroPaddToDims                         ( int ( numpy.max ( [ pStruct_static.getXDim(), pStruct_moving.getXDim() ] ) ),
                                                        int ( numpy.max ( [ pStruct_static.getYDim(), pStruct_moving.getYDim() ] ) ),
                                                        int ( numpy.max ( [ pStruct_static.getZDim(), pStruct_moving.getZDim() ] ) ) )
pStruct_moving.zeroPaddToDims                         ( int ( numpy.max ( [ pStruct_static.getXDim(), pStruct_moving.getXDim() ] ) ),
                                                        int ( numpy.max ( [ pStruct_static.getYDim(), pStruct_moving.getYDim() ] ) ),
                                                        int ( numpy.max ( [ pStruct_static.getZDim(), pStruct_moving.getZDim() ] ) ) )
                                                
######################################################
### Computing translation map
### =========================
###
### Now that all the preparations are done, the following
### function can be used to actually compute the
### translation function.
###
### NOTE: This function will fail if the two structures
### do not have the same dimensions and sampling, so
### for this reason it is required that the settings
### option for re-sampling the internal map to the
### same resolution and using the zero-padding above
### are used.
###
pStruct_moving.computeTranslationMap                  ( pStruct_static )

######################################################
### Accessing the translation function
### ==================================
###
### The translation function can be accessed
### from python in a very similar fashion as the
### rotationFunction or the internal map - in the
### three-dimensional numpy.ndarray format.
###
### The highest value of this translation function
### is the optimal global overlay translation
### vector for the two structures used to produce
### it.
###
translationFunction                                   = pStruct_moving.getTranslationFunctionMap ( )

### Find the map maximum
rotMapMax                                             = numpy.where ( translationFunction == numpy.amax ( translationFunction ) )

### Find maximum value
print                                                 ( "Translation map maximum is: " + str( translationFunction[rotMapMax[0][0]][rotMapMax[1][0]][rotMapMax[2][0]] ) )

### Expected output: Translation map maximum is: (2.160945931840744+0j)

######################################################
### Obtaining the optimal translation
### =================================
###
### The optimal translation vector is deeply connected
### to the centre of rotation, especially when we consider
### that the same vector needs to work for input maps
### as well as input co-ordinates. Therefore, ProSHADE
### reports the optimal translation as a dictionary with
### two entries:
###
### translationVecs["centreOfRotation"]
###     This vector is the translation required to move
###     the rotation centre to the origin (alternatively,
###     it can be seen as negative values of the rotation
###     centre itself).
###
### translationVecs["rotCenToOverlay"]
###     This is the translation vector from the centre
###     of rotation to the optimal overlay position.
###     This means that if the structure was moved to
###     origin for rotation, it needs to be moved back
###     to the original position and THEN this vector
###     needs to be applied.
###
translationVecs                                       = pStruct_moving.getOverlayTranslations ( pStruct_static  )

### Print the results
print                                                 ( "The centre of rotation is:                                " + str( -translationVecs["centreOfRotation"][0] ) + " ; " + str( -translationVecs["centreOfRotation"][1] ) + " ; " + str( -translationVecs["centreOfRotation"][2] ) )
print                                                 ( "The centre of rotation to optimal overlay translation is: " + str( translationVecs["rotCenToOverlay"][0] ) + " ; " + str( translationVecs["rotCenToOverlay"][1] ) + " ; " + str( translationVecs["rotCenToOverlay"][2] ) )

### Expected output
#   The centre of rotation is:                                -9.14285659790039 ; -17.77777862548828 ; -17.77777862548828
#   The centre of rotation to optimal overlay translation is: 16.0 ; 16.0 ; -8.0

######################################################
### Writing out the final structures
### ================================
###
### The rotated and translated map can be written out
### using the same function for writing out the internal
### map as the internal map is being modified by the
### functions.
###
### The case of co-ordinates, however, is different. These
### are outputted by the second showcased function, which
### firstly reads the co-ordinates from the original disk
### file, then applies the rotation around the rotation
### centre and then applies the translation as computed
### by the translation function. Finally, these new
### co-oridnates are written out.
###
pStruct_moving.writeMap                               ( "/Users/mysak/Desktop/movedPy.map" )
pStruct_moving.writePdb                               ( "/Users/mysak/Desktop/movedPy.pdb",
                                                        optimalRotationAngles[0],
                                                        optimalRotationAngles[1],
                                                        optimalRotationAngles[2],
                                                        translationVecs["rotCenToOverlay"][0],
                                                        translationVecs["rotCenToOverlay"][1],
                                                        translationVecs["rotCenToOverlay"][2],
                                                        translationVecs["centreOfRotation"][0],
                                                        translationVecs["centreOfRotation"][1],
                                                        translationVecs["centreOfRotation"][2],
                                                        pSet.firstModelOnly )

######################################################
### Clean up!
### =========
###
del pStruct_static
del pStruct_moving

######################################################
### Done
### ====
###
print ( "The end." )
