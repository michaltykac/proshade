##############################################
##############################################
# \file advancedAccess_reboxing.py
# \brief This code demonstrates the usage of the ProSHADE tool in the advanced mode.
#
# This file shows a fast demonstration of how the advanced access interfacte can be used to find the
# minimal boundaries containing all substantial density and how the results can be obtained. This file
# does not contain all the explanations and possible settings, for complete documentation, please see
# the advancedAccess.py file instead.
#
# Please note that the input structure paths are hard coded and you will need
# to change these before this file runs.
#
# This file is part of the ProSHADE library for calculating
# shape descriptors and symmetry operators of protein structures.
# This is a prototype code, which is by no means complete or fully
# tested. Its use is at your own risk only. There is no quarantee
# that the results are correct.
#
# \author    Michal Tykac
# \author    Garib N. Murshudov
# \version   0.7.3
# \date      MAY 2020
##############################################
##############################################


### System modules
import sys
import numpy

### Import ProSHADE from non-system folder (local installation assumed)
#sys.path.append                               ( "/Users/mysak/BioCEV/proshade/experimental/install/python2" )
import proshade

### Create the settings object
pSet                                          = proshade.ProSHADE_settings ()

### Set settings values
pSet.task                                     = proshade.MapManip
pSet.verbose                                  = 1

pSet.setMasking                               ( True )                                 # Should maps be masked by blurring?
pSet.setMapReboxing                           ( True )                                 # Should the structure be re-boxed? Required masking to be done in order to be meaningful.
pSet.setNormalisation                         ( False )                                # Should internal map representation be normalised to mean 0 and standard deviation 1?
pSet.setMapCentering                          ( False )                                # Move structure COM to the centre of map box?
pSet.setOutputFilename                        ( str ( "reBoxed" ) )                    # Filename to where re-boxed structure will be written to.
pSet.setMapResolutionChange                   ( False )                                # Should maps be re-sample to the computation resolution?
pSet.setBoundsSpace                           ( 3.0 )                                  # The extra space in Angs to add to the minimal boundaries when re-boxing.
pSet.setBoundsThreshold                       ( 5 )                                    # If two boundaries are within this threshold, the smaller one will be increased to have the same value as the larger one.
pSet.setSameBoundaries                        ( False )                                # Make multiple structures have the same boundaries. This is useful for half-maps.
pSet.setPDBBFactor                            ( -1.0 )                                 # Should all B-factors in a PDB file changed to this value? If no, set to negative value.
pSet.setMaskBlurFactor                        ( 350.0 )                                # If masking, what blur factor should be used? 350 seems to work for most maps.
pSet.setMaskIQR                               ( 3.0 )                                  # Number of inter-quartile ranges from median to use as the masking threshold.
pSet.setMaskSaving                            ( False )                                # Should map mask be saved?
pSet.setMaskFilename                          ( str ( "maskFile" ) )                   # The filename (no extension) to which the map masks will be saved into.


### Set example values
xDimIndices                                   = 100
yDimIndices                                   = 120
zDimIndices                                   = 60
xDimAngstroms                                 = xDimIndices * 1.3
yDimAngstroms                                 = yDimIndices * 1.3
zDimAngstroms                                 = zDimIndices * 1.3
xFrom                                         = int ( -xDimIndices/2 )
yFrom                                         = int ( -yDimIndices/2 )
zFrom                                         = int ( -zDimIndices/2 )
xTo                                           = int ( (xDimIndices/2)-1 )
yTo                                           = int ( (yDimIndices/2)-1 )
zTo                                           = int ( (zDimIndices/2)-1 )
ord                                           = 0

### Create example map (this will be a ball in the middle of the map)
testMap = numpy.empty ( [ ( xDimIndices * yDimIndices * zDimIndices ) ] )
for xIt in range( 0, xDimIndices ):
    for yIt in range( 0, yDimIndices ):
        for zIt in range( 0, zDimIndices ):
            ind = zIt + zDimIndices * ( yIt + yDimIndices * xIt )
            testMap[ind] = 1.0 / ( numpy.sqrt( numpy.power ( (xDimIndices/2) - xIt, 2.0 ) + numpy.power ( (yDimIndices/2) - yIt, 2.0 ) + numpy.power ( (zDimIndices/2) - zIt, 2.0 ) ) + 0.01 )

### Create the ProSHADE_data object without structure file on drive
pStruct                                       = proshade.ProSHADE_data ( pSet, "python_map_test", testMap, xDimAngstroms, yDimAngstroms, zDimAngstroms, xDimIndices, yDimIndices, zDimIndices, xFrom, yFrom, zFrom, xTo, yTo, zTo, ord )

### Process internal map
pStruct.processInternalMap                    ( pSet )

### Find the new boundaries (the 6 is required, but useless from the users point of view)
minimalBounds                                 = pStruct.getReBoxBoundariesPy ( pSet, 6 )

### Create a new structure, which will have the re-boxed boundaries
reBoxStr                                      = proshade.ProSHADE_data ( pSet )

### Fill the new structure with the calling structures map in the boudaries supplied
pStruct.createNewMapFromBoundsPy              ( pSet, reBoxStr, minimalBounds )

### Get the map into Python (as 1D array)
initialMapArray                               = reBoxStr.getMapPython ( reBoxStr.getMapArraySizePython() )

### How did the re-boxing go?
print ( "Original dimensions were: %+1.0f  x  %+1.0f  x  %+1.0f" % ( xDimIndices, yDimIndices, zDimIndices ) )
print ( "Current dimensions are:   %+1.0f  x  %+1.0f  x  %+1.0f" % ( reBoxStr.getXDim(), reBoxStr.getYDim(), reBoxStr.getZDim() ) )
perc                                          = ( ( reBoxStr.getXDim() * reBoxStr.getYDim() * reBoxStr.getZDim() ) / ( xDimIndices * yDimIndices * zDimIndices ) ) * 100.0
faster                                        = 100.0 / perc
print ( "Saved %+1.3f percents of indices and thus made the processing of the map %+1.3f times faster." % ( 100 - perc, faster ) )

### Print map average (just to show that the map in now in python)
print ( "Map average is %+1.3f" % ( numpy.mean ( initialMapArray ) ) )

### Release C++ pointers
del pStruct
del reBoxStr
del pSet

### Done
