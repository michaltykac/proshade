######################################################
######################################################
#   \file advancedAccess_overlay.py
#   brief This code demonstrates the usage of the ProSHADE tool in the advanced mode for overlay functionality.
#
#   This file demonstrates how the ProSHADE Overlay mode can be used from the Python module in the advanced mode. In general, the user needs to
#   create the two structures (static and moving) and process them as shown, making sure that the phase is removed. Then, it is possible to call
#   the getOverlayRotationFunction(), leading to optimal rotation matrix being returned. Next, the user needs to delete all the structural info,
#   re-load the structures and process them again, this time with phase information. From here, it is possible to rotate the moving map and obtain
#   the optimal translation vector using getOverlayTranslations() function. This file does not contain all the explanations and possible
#   settings, for complete documentation, please see the directAccess.py file instead.
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
pSet.task                                             = proshade.OverlayMap
pSet.verbose                                          = 4
pSet.requestedResolution                              = 2.0
pSet.usePhase                                         = False
pSet.changeMapResolution                              = True
pSet.maskMap                                          = False
pSet.moveToCOM                                        = False
pSet.normaliseMap                                     = False
pSet.reBoxMap                                         = False

######################################################
### Create data objects
pStruct_static                                        = proshade.ProSHADE_data ( pSet )
pStruct_moving                                        = proshade.ProSHADE_data ( pSet )

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
### Compute rotation function
pStruct_moving.getOverlayRotationFunction             ( pSet, pStruct_static )

######################################################
### Get the Euler angles
optimalRotationAngles                                 = pStruct_moving.getBestRotationMapPeaksEulerAngles ( pSet )
optimalRotationMatrix                                 = pStruct_moving.getBestRotationMapPeaksRotationMatrix ( pSet )

######################################################
### Delete phase-less data
del pStruct_static
del pStruct_moving

######################################################
### Change the settings - this is mandatory for Overlay task. Use the simple access if you do not want to deal with such things
pSet.usePhase                                         = True
pSet.changeMapResolution                              = True

######################################################
### Create objects
pStruct_static                                        = proshade.ProSHADE_data ( pSet )
pStruct_moving                                        = proshade.ProSHADE_data ( pSet )

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
### Rotate the moving structure (the signs with the Euler angles need to be like this!)
pStruct_moving.rotateMap                              ( pSet, optimalRotationAngles[0], optimalRotationAngles[1], optimalRotationAngles[2] )

######################################################
### Zero padding for both structures (only really applied to the smaller one, as nothing is added if the dimensions already have the requested size)
pStruct_static.zeroPaddToDims                         ( int ( numpy.max ( [ pStruct_static.getXDim(), pStruct_moving.getXDim() ] ) ),
                                                        int ( numpy.max ( [ pStruct_static.getYDim(), pStruct_moving.getYDim() ] ) ),
                                                        int ( numpy.max ( [ pStruct_static.getZDim(), pStruct_moving.getZDim() ] ) ) )
pStruct_moving.zeroPaddToDims                         ( int ( numpy.max ( [ pStruct_static.getXDim(), pStruct_moving.getXDim() ] ) ),
                                                        int ( numpy.max ( [ pStruct_static.getYDim(), pStruct_moving.getYDim() ] ) ),
                                                        int ( numpy.max ( [ pStruct_static.getZDim(), pStruct_moving.getZDim() ] ) ) )

######################################################
### Compute the translation map on the now identically sized and sampled maps
pStruct_moving.computeTranslationMap                  ( pStruct_static )

######################################################
### Find the translation vectors
translationVecs                                       = pStruct_moving.getOverlayTranslations ( pStruct_static,
                                                                                                optimalRotationAngles[0],
                                                                                                optimalRotationAngles[1],
                                                                                                optimalRotationAngles[2] )

######################################################
### Print results
print                                                 ( "Optimal Euler angles                :  %+1.3f    %+1.3f    %+1.3f" % ( optimalRotationAngles[0],
                                                                                                                                optimalRotationAngles[1],
                                                                                                                                optimalRotationAngles[2] ) )
print                                                 ( "Optimal Euler rotation matrix       :  %+1.3f    %+1.3f    %+1.3f" % ( optimalRotationMatrix[0][0],
                                                                                                                                optimalRotationMatrix[0][1],
                                                                                                                                optimalRotationMatrix[0][2] ) )
print                                                 ( "                                    :  %+1.3f    %+1.3f    %+1.3f" % ( optimalRotationMatrix[1][0],
                                                                                                                                optimalRotationMatrix[1][1],
                                                                                                                                optimalRotationMatrix[1][2] ) )
print                                                 ( "                                    :  %+1.3f    %+1.3f    %+1.3f" % ( optimalRotationMatrix[2][0],
                                                                                                                                optimalRotationMatrix[2][1],
                                                                                                                                optimalRotationMatrix[2][2] ) )
print                                                 ( "Rot. centre to origin translation   :  %+1.3f    %+1.3f    %+1.3f" % ( translationVecs["centreOfRotation"][0],
                                                                                                                                translationVecs["centreOfRotation"][1],
                                                                                                                                translationVecs["centreOfRotation"][2] ) )
print                                                 ( "Rot. centre to overlay translation  :  %+1.3f    %+1.3f    %+1.3f" % ( translationVecs["rotCenToOverlay"][0],
                                                                                                                                translationVecs["rotCenToOverlay"][1],
                                                                                                                                translationVecs["rotCenToOverlay"][2] ) )

######################################################
### Expected output
#   Optimal Euler angles                :  +5.465    +0.753    +3.894
#   Optimal Euler rotation matrix       :  -0.863    -0.079    +0.499
#                                       :  -0.192    -0.863    -0.467
#                                       :  +0.467    -0.499    +0.730
#   Rot. centre to origin translation   :  +18.000    +22.000    +24.000
#   Rot. centre to overlay translation  :  +4.000    +4.000    -5.000


######################################################
### Write out the map - it has some artefacts, this is caused by the double interpolation - I recommend applying the rotation and translation in EMDA instead of using this map.
pStruct_moving.writeMap                               ( "/Users/mysak/Desktop/movedPy.map" )
pStruct_moving.writePdb                               ( "/Users/mysak/Desktop/movedPy.pdb",
                                                        optimalRotationAngles[0],
                                                        optimalRotationAngles[1],
                                                        optimalRotationAngles[2],
                                                        translationVecs["rotCenToOverlay"][0],
                                                        translationVecs["rotCenToOverlay"][1],
                                                        translationVecs["rotCenToOverlay"][2],
                                                        pSet.firstModelOnly )

######################################################
### Delete the data
del pStruct_static
del pStruct_moving
