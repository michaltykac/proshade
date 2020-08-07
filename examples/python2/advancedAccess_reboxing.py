##############################################
##############################################
#   \file advancedAccess_reboxing.py
#   \brief This code demonstrates the usage of the ProSHADE tool in the advanced mode.
#
#   This file shows a fast demonstration of how the advanced access interfacte can be used to find the
#   minimal boundaries containing all substantial density and how the results can be obtained. This file
#   does not contain all the explanations and possible settings, for complete documentation, please see
#   the advancedAccess.py file instead.
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
#   \version   0.7.3
#   \date      AUG 2020
##############################################
##############################################


### System modules
import sys
import numpy

### Import ProSHADE from non-system folder (local installation assumed)
#sys.path.append                               ( "/Users/mysak/BioCEV/proshade/master/install/python2" )
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
perc                                          = ( float ( reBoxStr.getXDim() * reBoxStr.getYDim() * reBoxStr.getZDim() ) / ( xDimIndices * yDimIndices * zDimIndices ) ) * 100.0
faster                                        = 100.0 / perc
print ( "Saved %+1.3f percents of indices and thus made the processing of the map %+1.3f times faster." % ( 100 - perc, faster ) )

### Print map average (just to show that the map in now in python)
print ( "Map average is %+1.3f" % ( numpy.mean ( initialMapArray ) ) )

### Expected output
#   Original dimensions were: +100  x  +120  x  +60
#   Current dimensions are:   +42  x  +42  x  +42
#   Saved +89.710 percents of indices and thus made the processing of the map +9.718 times faster.
#   Map average is +0.031

### Release C++ pointers
del pStruct
del reBoxStr
del pSet

### Done
