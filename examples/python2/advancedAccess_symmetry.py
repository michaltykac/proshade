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
#   \version   0.7.4.3
#   \date      SEP 2020
##############################################
##############################################


### System modules
import sys
import numpy

### Import ProSHADE from non-system folder (local installation assumed)
sys.path.append                               ( "/Users/mysak/BioCEV/proshade/experimental/install/python2" )
import proshade

### Create the settings object
pSet                                          = proshade.ProSHADE_settings ()

### Set settings values
pSet.task                                     = proshade.Symmetry
pSet.verbose                                  = 1
pSet.setResolution                            ( 8.0 )
pSet.moveToCOM                                = False
pSet.changeMapResolution                      = True
#pSet.requestedSymmetryType                    = "C"
#pSet.requestedSymmetryFold                    = 12

### Create the structure object
pStruct                                       = proshade.ProSHADE_data ( pSet )

### Read in the structure
pStruct.readInStructure                       ( "./emd_6324.map", 0, pSet ) # This example uses EMD 6324 (PDB 3JA7)

### Process map
pStruct.processInternalMap                    ( pSet )

### Map to spheres
pStruct.mapToSpheres                          ( pSet )

### Compute spherical harmonics
pStruct.computeSphericalHarmonics             ( pSet )

### Compute self-rotation function
pStruct.getRotationFunction                   ( pSet )

### Detect recommended symmetry
pStruct.detectSymmetryInStructurePython       ( pSet )
recSymmetryType                               = pStruct.getRecommendedSymmetryType ( pSet )
recSymmetryFold                               = pStruct.getRecommendedSymmetryFold ( pSet )
recSymmetryAxes                               = proshade.getRecommendedSymmetryAxesPython ( pStruct, pSet )

### Print results
print ( "Detected " + str( recSymmetryType ) + "-" + str( recSymmetryFold ) + " symetry." )
print ( "Fold      x         y         z       Angle     Height" )
for iter in range ( 0, len( recSymmetryAxes ) ):
     print ( "  %s    %+1.3f    %+1.3f    %+1.3f    %+1.3f    %+1.4f" % ( recSymmetryAxes[iter][0], recSymmetryAxes[iter][1], recSymmetryAxes[iter][2], recSymmetryAxes[iter][3], recSymmetryAxes[iter][4], recSymmetryAxes[iter][5] ) )

### Expected output
#   Detected D-12 symetry.
#   Fold      x         y         z       Angle     Height
#     12    -0.004    +0.013    +1.000    +0.524    +0.9552
#      2    -0.190    +0.982    -0.009    +3.142    +0.3469

### Get list of all cyclic axes detected
allCAxes                                      = proshade.getAllDetectedSymmetryAxes ( pStruct, pSet )
print ( "Found a total of " + str( len ( allCAxes ) ) + " cyclic point groups." )

### Expected output
#   Found a total of 9 cyclic point groups.

### Get indices of which C axes form any detected non-C symmetry
allNonCAxesIndices                            = proshade.getNonCSymmetryAxesIndices ( pSet )
print ( "Found a total of " + str( len ( allNonCAxesIndices["D"] ) ) + " dihedral point groups." )

### Expected output
#   Found a total of 22 dihedral point groups.

### Get the group elements for the firs dihedral group - this does not have to be the recommended one, just the first in the list
firstAxisElements                             = proshade.getGroupElementsRotMat ( pStruct, pSet, allNonCAxesIndices["D"][0][0] )
secondAxisElements                            = proshade.getGroupElementsRotMat ( pStruct, pSet, allNonCAxesIndices["D"][0][1] )
allGroupElements                              = firstAxisElements + secondAxisElements
allGroupElements.insert                       ( 0, numpy.identity ( 3, dtype="float32" ) ) ### This is to add the identity element not returned by ProSHADE

### Print the first non-identity element
print ( "The first non-identity element is:" )
print ( "  %+1.3f    %+1.3f    %+1.3f " % ( allGroupElements[1][0][0], allGroupElements[1][0][1], allGroupElements[1][0][2] ) )
print ( "  %+1.3f    %+1.3f    %+1.3f " % ( allGroupElements[1][1][0], allGroupElements[1][1][1], allGroupElements[1][1][2] ) )
print ( "  %+1.3f    %+1.3f    %+1.3f " % ( allGroupElements[1][2][0], allGroupElements[1][2][1], allGroupElements[1][2][2] ) )

### Expected output
#   The first non-identity element is:
#     +0.866    -0.500    +0.006
#     +0.500    +0.866    +0.004
#     -0.007    -0.000    +1.000

### Release C++ pointers
del pStruct
del pSet

### Done
