##############################################
##############################################
#   \file advancedAccess_symmetry.py
#   \brief This code demonstrates the usage of the ProSHADE tool in the advanced mode.
#
#   This file shows a fast demonstration of how the advanced access interfacte can be used to compute the
#   symmetry of a particular structure and how the results can be obtained. This file does not contain all
#   the explanations and possible settings, for complete documentation, please see the advancedAccess.py
#   file instead.
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
sys.path.append                               ( "/Users/mysak/BioCEV/proshade/development/install/python3" )
import proshade

### Create the settings object
pSet                                          = proshade.ProSHADE_settings ()

### Set settings values
pSet.task                                     = proshade.Symmetry
pSet.verbose                                  = 1
pSet.setResolution                            ( 12.0 )
pSet.moveToCOM                                = False
pSet.changeMapResolution                      = True
pSet.changeMapResolutionTriLinear             = False
pSet.requestedSymmetryType                    = "C"
pSet.requestedSymmetryFold                    = 12

### Create the structure object
pStruct                                       = proshade.ProSHADE_data ( pSet )

### Read in the structure
pStruct.readInStructure                       ( "./emd_6324.map", 0, pSet ) # The path to the structure to be processed. This example uses EMD 6324 (PDB 3JA7)

### Process map
pStruct.processInternalMap                    ( pSet )

### Map to spheres
pStruct.mapToSpheres                          ( pSet )

### Compute spherical harmonics
pStruct.computeSphericalHarmonics             ( pSet )

### Compute self-rotation function
pStruct.getRotationFunction                   ( pSet )

### Detect symmetry
pStruct.detectSymmetryInStructurePython       ( pSet )
symmetryType                                  = pStruct.getRecommendedSymmetryType ( pSet )
symmetryFold                                  = pStruct.getRecommendedSymmetryFold ( pSet )
symmetryAxes                                  = proshade.getSymmetryAxesPython ( pStruct, pSet )

print ( "Detected " + str( symmetryType ) + "-" + str( symmetryFold ) + " symetry." )
print ( "Fold      x         y         z       Angle     Height" )
for iter in range ( 0, len( symmetryAxes ) ):
     print ( "  %s    %+1.3f    %+1.3f    %+1.3f    %+1.3f    %+1.4f" % ( symmetryAxes[iter][0], symmetryAxes[iter][1], symmetryAxes[iter][2], symmetryAxes[iter][3], symmetryAxes[iter][4], symmetryAxes[iter][5] ) )

### Expected output
#   Detected C-12 symetry.
#   Fold      x         y         z       Angle     Height
#     12    -0.012    +0.004    +1.000    +0.524    +0.9621

### Release C++ pointers
del pStruct
del pSet

### Done
