######################################################
######################################################
#   \file advancedAccess_symmetry.py
#   \brief This code demonstrates the usage of the ProSHADE tool in the advanced mode for the symmetry detection functionality.
#
#   This file shows a fast demonstration of how the advanced access interfacte can be used to compute the
#   symmetry of a particular structure and how the results can be obtained. This file does not contain all
#   the explanations and possible settings, for complete documentation, please see the directAccess.py
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
#   \version   0.7.6.2
#   \date      DEC 2021
######################################################
######################################################

######################################################
### Import system modules
import sys
import numpy

######################################################
### Import ProSHADE from non-system folder (local installation assumed)
sys.path.append                                       ( "/Users/mysak/BioCEV/proshade/experimental/proshade/install/pythonModule" )
import proshade

######################################################
### Create the settings object
pSet                                                  = proshade.ProSHADE_settings ( )

######################################################
### Set basic settings values
pSet.task                                             = proshade.Symmetry
pSet.setResolution                                    ( 12.0 )
pSet.setMapResolutionChange                           ( True )
pSet.setMapCentering                                  ( True )
pSet.verbose                                          = -1

######################################################
### Create the structure object
pStruct                                               = proshade.ProSHADE_data ( )

######################################################
### Read in the structure
pStruct.readInStructure                               ( "/Users/mysak/BioCEV/proshade/playground/emd_0116.map.gz", 0, pSet ) # This example uses EMD 0116

######################################################
### Process map
pStruct.processInternalMap                            ( pSet )

######################################################
### Map to spheres
pStruct.mapToSpheres                                  ( pSet )

######################################################
### Compute spherical harmonics
pStruct.computeSphericalHarmonics                     ( pSet )

######################################################
### Compute self-rotation function
pStruct.computeRotationFunction                       ( pSet )

######################################################
### Symmetry detection algorithm run
pStruct.detectSymmetryInStructure                     ( pSet )
    
######################################################
### Retrieve the results
recSymmetryType                                       = pStruct.getRecommendedSymmetryType ( pSet )
recSymmetryFold                                       = pStruct.getRecommendedSymmetryFold ( pSet )
recSymmetryAxes                                       = pStruct.getRecommendedSymmetryAxes ( pSet )

######################################################
### Print results
print                                                 ( "Detected " + str( recSymmetryType ) + "-" + str( recSymmetryFold ) + " symetry." )
print                                                 ( "Fold              x           y         z       Angle      Height    Average FSC" )
for iter in range ( 0, len( recSymmetryAxes ) ):
     print                                            ( "  %s    \t%+1.3f      %+1.3f    %+1.3f    %+1.3f    %+1.4f      %+1.4f" % ( recSymmetryAxes[iter][0],
                                                                                                                         recSymmetryAxes[iter][1],
                                                                                                                         recSymmetryAxes[iter][2],
                                                                                                                         recSymmetryAxes[iter][3],
                                                                                                                         recSymmetryAxes[iter][4],
                                                                                                                         recSymmetryAxes[iter][5],
                                                                                                                         recSymmetryAxes[iter][6] ) )

######################################################
### Expected output
#   Detected D-6 symetry.
#   Fold              x           y         z       Angle      Height    Average FSC
#     6.0        +0.000      +0.000    +1.000    +1.047    +0.9709      +0.9758
#     2.0        +1.000      -0.005    -0.000    -3.142    +1.0103      +0.8624

######################################################
### Get list of all cyclic axes detected
allCAxes                                              = pStruct.getAllCSyms ( pSet )
print                                                 ( "Found a total of " + str( len ( allCAxes ) ) + " cyclic point groups." )

######################################################
### Expected output
#   Found a total of 9 cyclic point groups.

######################################################
### Get indices of which C axes form any detected non-C symmetry
allNonCAxesIndices                                    = pStruct.getNonCSymmetryAxesIndices ( pSet )
print                                                 ( "Found a total of " + str( len ( allNonCAxesIndices["D"] ) ) + " dihedral point groups." )

######################################################
### Expected output
#   Found a total of 21 dihedral point groups.

######################################################
### Find internal map COM shift
comShift                                              = pStruct.getMapCOMProcessChange ()
print                                                 ( "The internal map has been shifted by " + str( comShift[0] ) + " , " + str( comShift[1] ) + " , " + str( comShift[2] ) )

######################################################
### Expected output
#   The internal map has been shifted by 4.111518 , 4.1201262 , 2.39567

######################################################
#  NOTE: To get all the point group elements, one needs to supply the list of all cyclic point groups which comprise the
#        requested point group. This is relatively simple for T, O and I symmetries, as such a list is already produced by
#        ProSHADE - see the following examples:
#
#        allGroupElements = proshade.getAllGroupElements ( pSet, pStruct, allNonCAxesIndices['T'], "T" )
#        allGroupElements = proshade.getAllGroupElements ( pSet, pStruct, allNonCAxesIndices['O'], "O" )
#        allGroupElements = proshade.getAllGroupElements ( pSet, pStruct, allNonCAxesIndices['I'], "I" )
#
#        For cyclic point groups, this is also simple, as one can select the required >index< from the allCs variable and use
#        NOTE: The [] around index is required, as the function expects an array (list) and not an int!
#
#        allGroupElements = proshade.getAllGroupElements ( pSet, pStruct, [index], "C" )
#
#        The only problem comes when D is to be used, as ProSHADE gives a vector (list) of all combinations (also as vector/list) of cyclic point groups which form the
#        D point groups. Therefore, to select the recommended D point group from this list, a search needs to be done. This is shown in the following code.
######################################################

######################################################
### Define isclose() for comparing floats
def isclose ( a, b, tol = 1e-04 ):
    if abs(a) - tol < abs(b) and abs(a) + tol > abs(b):
        return True
    else:
        return False

######################################################
### Find the indices of the best dihedral combination
bestDCombination                                      = []
for dIt in range ( 0, len ( allNonCAxesIndices['D'] ) ):
    firstMatch                                        = False
    secondMatch                                       = False
    for recIt in range ( 0, len( recSymmetryAxes ) ):
        if ( isclose ( allCAxes[allNonCAxesIndices['D'][dIt][0]][1], recSymmetryAxes[0][1] ) ) and \
           ( isclose ( allCAxes[allNonCAxesIndices['D'][dIt][0]][2], recSymmetryAxes[0][2] ) ) and \
           ( isclose ( allCAxes[allNonCAxesIndices['D'][dIt][0]][3], recSymmetryAxes[0][3] ) ) and \
           ( isclose ( allCAxes[allNonCAxesIndices['D'][dIt][0]][4], recSymmetryAxes[0][4] ) ) and \
           ( isclose ( allCAxes[allNonCAxesIndices['D'][dIt][0]][5], recSymmetryAxes[0][5] ) ):
            firstMatch                                = True

        if ( isclose ( allCAxes[allNonCAxesIndices['D'][dIt][1]][1], recSymmetryAxes[1][1] ) ) and \
           ( isclose ( allCAxes[allNonCAxesIndices['D'][dIt][1]][2], recSymmetryAxes[1][2] ) ) and \
           ( isclose ( allCAxes[allNonCAxesIndices['D'][dIt][1]][3], recSymmetryAxes[1][3] ) ) and \
           ( isclose ( allCAxes[allNonCAxesIndices['D'][dIt][1]][4], recSymmetryAxes[1][4] ) ) and \
           ( isclose ( allCAxes[allNonCAxesIndices['D'][dIt][1]][5], recSymmetryAxes[1][5] ) ):
            secondMatch                               = True

    if (firstMatch and secondMatch) and (len(bestDCombination)==0):
        bestDCombination.append                       ( allNonCAxesIndices['D'][dIt][0] )
        bestDCombination.append                       ( allNonCAxesIndices['D'][dIt][1] )


######################################################
### Get the group elements for the best dihedral group
allGroupElements                                      = pStruct.getAllGroupElements ( pSet, bestDCombination, "D" )

######################################################
### Print the first non-identity element
print                                                 ( "Found a total of " + str( len ( allGroupElements ) ) + " group " + str( bestDCombination ) + " elements." )
print                                                 ( "The first non-identity element is:" )
print                                                 ( "  %+1.3f    %+1.3f    %+1.3f " % ( allGroupElements[1][0][0], allGroupElements[1][0][1], allGroupElements[1][0][2] ) )
print                                                 ( "  %+1.3f    %+1.3f    %+1.3f " % ( allGroupElements[1][1][0], allGroupElements[1][1][1], allGroupElements[1][1][2] ) )
print                                                 ( "  %+1.3f    %+1.3f    %+1.3f " % ( allGroupElements[1][2][0], allGroupElements[1][2][1], allGroupElements[1][2][2] ) )

######################################################
### Expected output
#   Found a total of 12 group [6, 0] elements.
#   The first non-identity element is:
#     +0.500    +0.866    +0.000 
#     -0.866    +0.500    +0.000 
#     +0.000    +0.000    +1.000 

######################################################
### Assuming you have modifiex / created your own C groups
### and want to obtain their group elements:
gr1Elements                                           = proshade.computeGroupElementsForGroup ( 1,      # X-axis element of the group rotation axis
                                                                                                0,      # Y-axis element of the group rotation axis
                                                                                                0,      # Z-axis element of the group rotation axis
                                                                                                4 )     # Fold

gr2Elements                                           = proshade.computeGroupElementsForGroup ( 0,      # X-axis element of the group rotation axis
                                                                                                1,      # Y-axis element of the group rotation axis
                                                                                                0,      # Z-axis element of the group rotation axis
                                                                                                2 )     # Fold

combinedGroupElements                                 = proshade.joinElementsFromDifferentGroups ( gr1Elements, # Elements of the first group
                                                                                                   gr2Elements, # Elements of the second group
                                                                                                   0.0001,      # Matrix tolerance (i.e. if the abs ( trace ( Mat1 * Mat2^(T) ) - 3.0 ) must be below this number for two matrices to be considered identical )
                                                                                                   True )       # Should new group elements resulting from multiplication of elements from group 1 with elements from group 2 be computed?

print ( len ( combinedGroupElements ) )

######################################################
### Expected output
#   8

######################################################
### Release C++ pointers
del pStruct
del pSet
