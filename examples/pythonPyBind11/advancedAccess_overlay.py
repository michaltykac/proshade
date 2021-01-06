##############################################
##############################################
#   \file advancedAccess_overlay.py
#   brief This code demonstrates the usage of the ProSHADE tool in the advanced mode for overlay functionality.
#
#   This file demonstrates how the ProSHADE Overlay mode can be used from the Python module in the advanced mode. In general, the user needs to
#   create the two structures (static and moving) and process them as shown, making sure that the phase is removed. Then, it is possible to call
#   the getOverlayRotationFunction(), leading to optimal rotation matrix being returned. Next, the user needs to delete all the structural info,
#   re-load the structures and process them again, this time with phase information. From here, it is possible to rotate the moving map and obtain
#   the optimal translation vector using getBestTranslationMapPeaksAngstrom() function. This file does not contain all the explanations and possible
#   settings, for complete documentation, please see the advancedAccess.py file instead.
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
#   \version   0.7.5.0
#   \date      DEC 2020
##############################################
##############################################

### System modules
import sys
import numpy

### Import ProSHADE from non-system folder (local installation assumed)
sys.path.append                               ( "/Users/mysak/BioCEV/proshade/build" )
import pyproshade as proshade

### Create the object
pSet                                          = proshade.ProSHADE_settings ()

### Set settings values for optimal overlay
pSet.task                                     = proshade.OverlayMap
pSet.verbose                                  = 4
pSet.requestedResolution                      = 4.0;
pSet.usePhase                                 = False;
pSet.changeMapResolution                      = True;
pSet.maskMap                                  = False;
pSet.moveToCOM                                = False;
pSet.normaliseMap                             = False;
pSet.reBoxMap                                 = False;

### Create objects
pStruct_static                                = proshade.ProSHADE_data ( pSet )
pStruct_moving                                = proshade.ProSHADE_data ( pSet )

### Read in the structures
pStruct_static.readInStructure                ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/bf/1BFO_A_dom_1.pdb", 0, pSet ) # This is a BALBES domain 1BFO_A_dom_1.
pStruct_moving.readInStructure                ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/h8/1H8N_A_dom_1.pdb", 1, pSet ) # This is a BALBES domain 1H8N_A_dom_1.

### Get spherical harmonics for both structures
pStruct_static.processInternalMap             ( pSet )
pStruct_moving.processInternalMap             ( pSet )

pStruct_static.mapToSpheres                   ( pSet )
pStruct_moving.mapToSpheres                   ( pSet )

pStruct_static.computeSphericalHarmonics      ( pSet )
pStruct_moving.computeSphericalHarmonics      ( pSet )

### Compute rotation function
pStruct_moving.getOverlayRotationFunction     ( pSet, pStruct_static )

### Get the Euler angles
optimalRotationAngles                         = pStruct_moving.getBestRotationMapPeaksEulerAngles ( pSet )
optimalRotationMatrix                         = pStruct_moving.getBestRotationMapPeaksRotationMatrix ( pSet )

### Delete phase-less data
del pStruct_static
del pStruct_moving

### Change the settings - this is mandatory for Overlay task. Use the simple access if you do not want to deal with such things
pSet.usePhase                                 = True
pSet.changeMapResolution                      = True

### Create objects
pStruct_static                                = proshade.ProSHADE_data ( pSet )
pStruct_moving                                = proshade.ProSHADE_data ( pSet )

### Read in the structures
pStruct_static.readInStructure                ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/bf/1BFO_A_dom_1.pdb", 0, pSet ) # This is a BALBES domain 1BFO_A_dom_1.
pStruct_moving.readInStructure                ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/h8/1H8N_A_dom_1.pdb", 1, pSet ) # This is a BALBES domain 1H8N_A_dom_1.

### Get spherical harmonics for moving structure only
pStruct_static.processInternalMap             ( pSet )

pStruct_moving.processInternalMap             ( pSet )
pStruct_moving.mapToSpheres                   ( pSet )
pStruct_moving.computeSphericalHarmonics      ( pSet )

### Rotate the moving structure (the signs with the Euler angles need to be like this!)
pStruct_moving.rotateMap                      ( pSet, optimalRotationAngles[0], optimalRotationAngles[1], optimalRotationAngles[2] )

### Zero padding for both structures (only really applied to the smaller one, as nothing is added if the dimensions already have the requested size)
pStruct_static.zeroPaddToDims                 ( int ( numpy.max ( [ pStruct_static.getXDim(), pStruct_moving.getXDim() ] ) ),
                                                int ( numpy.max ( [ pStruct_static.getYDim(), pStruct_moving.getYDim() ] ) ),
                                                int ( numpy.max ( [ pStruct_static.getZDim(), pStruct_moving.getZDim() ] ) ) )
pStruct_moving.zeroPaddToDims                 ( int ( numpy.max ( [ pStruct_static.getXDim(), pStruct_moving.getXDim() ] ) ),
                                                int ( numpy.max ( [ pStruct_static.getYDim(), pStruct_moving.getYDim() ] ) ),
                                                int ( numpy.max ( [ pStruct_static.getZDim(), pStruct_moving.getZDim() ] ) ) )

### Find the translation function
pStruct_moving.computeTranslationMap          ( pStruct_static )

### Get optimal translation vector
optimalTranslationVector                      = pStruct_moving.getBestTranslationMapPeaksAngstrom ( pStruct_static )
allTranslations                               = pStruct_moving.getOverlayTranslations ( pStruct_static )
toOrigin                                      = allTranslations["centreOfRotation"]
toMapCen                                      = allTranslations["withinBoxTranslations"]
toOverlay                                     = allTranslations["originToOverlay"]

### Print results
print ( "Optimal Euler angles         :  %+1.3f    %+1.3f    %+1.3f\n" % ( optimalRotationAngles[0], optimalRotationAngles[1], optimalRotationAngles[2] ) )
print ( "Optimal Euler rotation matrix:  %+1.3f    %+1.3f    %+1.3f" % ( optimalRotationMatrix[0][0], optimalRotationMatrix[0][1], optimalRotationMatrix[0][2] ) )
print ( "                             :  %+1.3f    %+1.3f    %+1.3f" % ( optimalRotationMatrix[1][0], optimalRotationMatrix[1][1], optimalRotationMatrix[1][2] ) )
print ( "                             :  %+1.3f    %+1.3f    %+1.3f\n" % ( optimalRotationMatrix[2][0], optimalRotationMatrix[2][1], optimalRotationMatrix[2][2] ) )
print ( "To origin translation        :  %+1.3f    %+1.3f    %+1.3f" % ( toOrigin[0], toOrigin[1], toOrigin[2] ) )
print ( "To map centre translation    :  %+1.3f    %+1.3f    %+1.3f" % ( toMapCen[0], toMapCen[1], toMapCen[2] ) )
print ( "To overlat translation       :  %+1.3f    %+1.3f    %+1.3f" % ( toOverlay[0], toOverlay[1], toOverlay[2] ) )

### Expected output
#   Optimal Euler angles         :  +5.433    +0.753    +3.927
#
#   Optimal Euler rotation matrix:  -0.872    -0.191    +0.451
#                                :  -0.078    -0.854    -0.514
#                                :  +0.483    -0.483    +0.730
#
#   To origin translation        :  +0.000    +0.000    +0.000
#   To map centre translation    :  +0.000    +0.000    +0.000
#   To overlat translation       :  +8.000    +8.000    +8.000

### Write out the map - it has some artefacts, this is caused by the double interpolation - I recommend applying the rotation and translation in EMDA instead of using this map.
pStruct_moving.translateMap                   ( pSet, optimalTranslationVector[0], optimalTranslationVector[1], optimalTranslationVector[2] );
pStruct_moving.writeMap                       ( "/Users/mysak/Desktop/movedPy.map" )
pStruct_moving.writePdb                       ( "/Users/mysak/Desktop/movedPy.pdb", optimalRotationAngles[0], optimalRotationAngles[1], optimalRotationAngles[2], optimalTranslationVector[0], optimalTranslationVector[1], optimalTranslationVector[2] )

### Delete the data
del pStruct_static
del pStruct_moving

### Done