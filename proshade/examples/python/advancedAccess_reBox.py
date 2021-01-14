######################################################
######################################################
#   \file advancedAccess_reBox.py
#   \brief This code demonstrates the usage of the ProSHADE tool in the advanced mode for re-boxing functionality.
#
#   This file shows a fast demonstration of how the advanced access interface can be used to re-box a density map.
#   It starts by loading and setting up the ProSHADE module and then it proceeds to create a density map using
#   the formula for three-dimensional ball. This map is then loaded into ProSHADE_data class and this object is then
#   processed to find the minimum boundaries which contain all the "significant" density (as defined by the ProSHADE
#   procedure). Finally, a new ProSHADE_data object with the smaller map is created and the map is transferred back
#   to python, concluding the example.
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
#   \version   0.7.5.1
#   \date      JAN 2021
######################################################
######################################################

######################################################
### Import system modules
import sys
import numpy

######################################################
### Import ProSHADE from non-system folder (local installation assumed)
sys.path.append                                       ( "/Users/mysak/BioCEV/proshade/experimental/install/pythonModule" )
import proshade

######################################################
### Create the settings object
pSet                                                  = proshade.ProSHADE_settings ()

######################################################
### Set basic settings values
pSet.task                                             = proshade.MapManip
pSet.verbose                                          = -1
        
pSet.setMasking                                       ( True )
pSet.setMapReboxing                                   ( True )
pSet.setNormalisation                                 ( False )
pSet.setMapCentering                                  ( False )
pSet.setOutputFilename                                ( str ( "reBoxed" ) )
pSet.setMapResolutionChange                           ( False )

######################################################
### Create example map (this will be a ball in the middle of the map)
xDimIndices                                           = 100
yDimIndices                                           = 120
zDimIndices                                           = 60
xDimAngstroms                                         = xDimIndices * 1.3
yDimAngstroms                                         = yDimIndices * 1.3
zDimAngstroms                                         = zDimIndices * 1.3
xFrom                                                 = int ( -xDimIndices/2 )
yFrom                                                 = int ( -yDimIndices/2 )
zFrom                                                 = int ( -zDimIndices/2 )
xTo                                                   = int ( (xDimIndices/2)-1 )
yTo                                                   = int ( (yDimIndices/2)-1 )
zTo                                                   = int ( (zDimIndices/2)-1 )
ord                                                   = 0
testMap = numpy.empty ( [ ( xDimIndices * yDimIndices * zDimIndices ) ] )
for xIt in range( 0, xDimIndices ):
    for yIt in range( 0, yDimIndices ):
        for zIt in range( 0, zDimIndices ):
            ind = zIt + zDimIndices * ( yIt + yDimIndices * xIt )
            testMap[ind] = 1.0 / ( numpy.sqrt( numpy.power ( (xDimIndices/2) - xIt, 2.0 ) + numpy.power ( (yDimIndices/2) - yIt, 2.0 ) + numpy.power ( (zDimIndices/2) - zIt, 2.0 ) ) + 0.01 )

######################################################
### Create the ProSHADE_data object without structure file on drive
pStruct                                               = proshade.ProSHADE_data ( pSet,
                                                                                 "python_map_test",
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
### Process internal map
pStruct.processInternalMap                            ( pSet )

######################################################
### Find the new boundaries
minimalBounds                                         = pStruct.getReBoxBoundaries ( pSet )

######################################################
### Create a new structure, which will have the re-boxed boundaries
reBoxStr                                              = proshade.ProSHADE_data ( pSet )

######################################################
### Fill the new structure with the calling structures map in the boudaries supplied
pStruct.createNewMapFromBounds                        ( minimalBounds, reBoxStr, pSet )

######################################################
### Get the map into Python
initialMap                                            = reBoxStr.getMap ( )

######################################################
### How did the re-boxing go?
print                                                 ( "Original dimensions were: %+1.0f  x  %+1.0f  x  %+1.0f" % ( xDimIndices, yDimIndices, zDimIndices ) )
print                                                 ( "Current dimensions are:   %+1.0f  x  %+1.0f  x  %+1.0f" % ( reBoxStr.getXDim(), reBoxStr.getYDim(), reBoxStr.getZDim() ) )
perc                                                  = ( ( reBoxStr.getXDim() * reBoxStr.getYDim() * reBoxStr.getZDim() ) / ( xDimIndices * yDimIndices * zDimIndices ) ) * 100.0
faster                                                = 100.0 / perc
print                                                 ( "Saved %+1.3f percents of indices and thus made the processing of the map %+1.3f times faster." % ( 100 - perc, faster ) )

######################################################
### Print map type and average (just to show that the map in now in python)
print ( "Map type is " + str( initialMap.ndim ) + "D " + str( type ( initialMap ) ) + " with average being %+1.3f" % ( numpy.mean ( initialMap ) ) )

### Expected output
#   Original dimensions were: +100  x  +120  x  +60
#   Current dimensions are:   +42  x  +42  x  +42
#   Saved +89.710 percents of indices and thus made the processing of the map +9.718 times faster.
#   Map type is 3D <class 'numpy.ndarray'> with average being +0.031

### Release C++ pointers
del initialMap
del pStruct
del reBoxStr
del pSet
