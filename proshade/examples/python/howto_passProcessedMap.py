######################################################
######################################################
#   \file howto_passProcessedMap.py
#   \brief This code demonstrates how used obtained numpy array map can be supplied to ProSHADE.
#
#   This file starts by importing all the required modules and reading in a map (as a 3D numpy array), its dimensions in Angstroms and the starting indices. It then
#   proceeds to read in the mask for this map, making the assumption that these have the same dimensions (this is an example, so no sanity check is done). The mask and
#   map are then simply multiplied together to simulate masking.
#
#   The next section then shows how ProSHADE can be set (using the settings object) and supplied with the resulting masked map along with some of the map information read
#   in previously to create a full-fledged data object. This object can then be used as if the readInStructure () function was called - in this case, the symmetry detection
#   is demonstrated with the results being simply outputted out.
#
#   Copyright by Michal Tykac and individual contributors. All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#   1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#   2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#   3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
#   This software is provided by the copyright holders and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
#
#   \author    Michal Tykac
#   \author    Garib N. Murshudov
#   \version   0.7.6.6
#   \date      JUL 2022
######################################################
######################################################

######################################################
### Import system modules
import sys
import numpy

######################################################
### Import support modules
import mrcfile

######################################################
### Import ProSHADE
import proshade

######################################################
### Read in map from file (to have some example
### values, replace this with any other way you like)
with mrcfile.open('/Users/mysak/BioCEV/proshade/playground/emd_0011.map') as mrc:
    mapArr   = mrc.data
    mapSizes = numpy.array ( [ mrc.header.cella['x'], mrc.header.cella['y'], mrc.header.cella['z'] ] )
    mapFroms = numpy.array ( [ mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart ] )
    mrc.close()

######################################################
### Create mask array. Here we assume same size as map
with mrcfile.open('/Users/mysak/BioCEV/proshade/playground/emd_0011_mask.map') as mrc:
    maskArr = mrc.data
    mrc.close()

######################################################
### Process map to your liking, here, I will just
### multiply map with mask
maskedMap = mapArr * maskArr

######################################################
### Create the settings object
pSet                                                  = proshade.ProSHADE_settings ( proshade.Symmetry )

######################################################
### Set basic settings values and turn all map
### processing off.
pSet.requestedResolution                              = 8.0   ## This is the computation limit and is required.
pSet.changeMapResolution                              = True  ## Should the map be re-sampled to the resolution? If True, computation will be faster.
pSet.invertMap                                        = False ## Map has the correct hand orientation
pSet.normaliseMap                                     = False ## Map will be used as it is, no normalisation required.
pSet.maskMap                                          = False ## Map will be used as it is, no extra internal masking required.
pSet.moveToCOM                                        = True  ## Map will be used as it is, no centering required.
pSet.usePhase                                         = True  ## Map will be used as it is, not interested in Patterson.
pSet.verbose                                          = -1    ## How loud should the run be?

######################################################
### Set all the values required to pass map to
### proshade directly
xDimIndices                                           = maskedMap.shape[0]
yDimIndices                                           = maskedMap.shape[1]
zDimIndices                                           = maskedMap.shape[2]
xDimAngstroms                                         = mapSizes[0]
yDimAngstroms                                         = mapSizes[1]
zDimAngstroms                                         = mapSizes[2]
xFrom                                                 = mapFroms[0]
yFrom                                                 = mapFroms[1]
zFrom                                                 = mapFroms[2]
xTo                                                   = mapFroms[0] + xDimIndices - 1
yTo                                                   = mapFroms[1] + yDimIndices - 1
zTo                                                   = mapFroms[2] + zDimIndices - 1
ord                                                   = 0

######################################################
### Create the ProSHADE_data object from the array and
### all the other information
pStruct                                               = proshade.ProSHADE_data ( "python_map_test",
                                                                                 maskedMap,
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


#pStruct.reSampleMap ( pSet )

######################################################
### Proceed with the structure as usual - this example
### is for symmetry detection
pStruct.processInternalMap                            ( pSet )
pStruct.mapToSpheres                                  ( pSet )
pStruct.computeSphericalHarmonics                     ( pSet )
pStruct.computeRotationFunction                       ( pSet )
pStruct.detectSymmetryInStructure                     ( pSet )
    
######################################################
### Retrieve results
recSymmetryType                                       = pStruct.getRecommendedSymmetryType ( )
recSymmetryFold                                       = pStruct.getRecommendedSymmetryFold ( )
recSymmetryAxes                                       = pStruct.getRecommendedSymmetryAxes ( )
allCAxes                                              = pStruct.getAllCSyms ( )

######################################################
### Print results
print ( "Recommended symmetry:" )
print ( str( recSymmetryType ) + str( recSymmetryFold ) )
print ( "" )
print ( "Recommended axes:" )
for i in range( 0, len( recSymmetryAxes ) ):
    print ( str( recSymmetryAxes[i][0] ) + " | " + str( recSymmetryAxes[i][1] ) + " x " + str( recSymmetryAxes[i][2] ) + " x " + str( recSymmetryAxes[i][3] ) + " || " + str( recSymmetryAxes[i][6] ) )
print ( "" )
print ( "All axes:" )
for i in range( 0, len( allCAxes ) ):
    print ( str( allCAxes[i][0] ) + " | " + str( allCAxes[i][1] ) + " x " + str( allCAxes[i][2] ) + " x " + str( allCAxes[i][3] ) + " || " + str( allCAxes[i][6] ) )

### EXPECTED OUTPUT: Recommended symmetry:
### EXPECTED OUTPUT: D3
### EXPECTED OUTPUT: 
### EXPECTED OUTPUT: Recommended axes:
### EXPECTED OUTPUT: 3.0 | 0.9999328 x -0.011591071 x -6.123234e-17 || 0.9127528
### EXPECTED OUTPUT: 2.0 | 0.0007012483 x 0.0 x 0.99999976 || 0.998042
### EXPECTED OUTPUT: 
### EXPECTED OUTPUT: All axes:
### EXPECTED OUTPUT: 2.0 | 0.0007012483 x 0.0 x 0.99999976 || 0.998042
### EXPECTED OUTPUT: 3.0 | 0.9999328 x -0.011591071 x -6.123234e-17 || 0.9127528
### EXPECTED OUTPUT: 2.0 | 0.009986875 x 0.8641762 x -0.5030902 || 0.90472996
### EXPECTED OUTPUT: 2.0 | 0.009988657 x 0.8641767 x 0.5030893 || 0.9047275

######################################################
### Release ProSHADE memory
del pStruct
del pSet

######################################################
### Done 
