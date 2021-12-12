######################################################
######################################################
#   \file howto_useGemmiObject.py
#   \brief This code demonstrates how gemmi::Structure object can be used to create ProSHADE_data object from python.
#
#   This code demonstrates how gemmi::Structure object (in this case read from file) can be passed to ProSHADE instead
#   of letting ProSHADE read the data from file. This could be useful if the gemmi::Structure object was modified by the
#   user between it being read in and passed to ProSHADE. 
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
### Import ProSHADE from system folder
import proshade

######################################################
### Import Gemmi
import gemmi

######################################################
### Create the settings object
pSet                                                  = proshade.ProSHADE_settings ( proshade.Distances )

######################################################
### Set basic settings values
pSet.requestedResolution                              = 1.0
pSet.changeMapResolution                              = True
pSet.moveToCOM                                        = False
pSet.verbose                                          = -1

######################################################
### Create Gemmi structure
gStruct                                               = gemmi.read_structure( "/Users/mysak/BioCEV/proshade/xx_Metals/Sandwich_5_3_times10.pdb" )

######################################################
### Find position of a metal atom (this needs to be
### improved if used outside of an example)
metalCoords                                           = numpy.full ( ( 3 ), numpy.inf )
for model in gStruct:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.element.is_metal:
                    metalCoords[0]                    = atom.pos[0]
                    metalCoords[1]                    = atom.pos[1]
                    metalCoords[2]                    = atom.pos[2]
                   
######################################################
### Force the metal co-ordinates to be at the centre
### of the map box.
###
### NOTE: This shifts the map so that it no longer
### overlays with the co-ordinate file in visualisation.
pSet.setBoxCentering                                  ( metalCoords[0], metalCoords[1], metalCoords[2] )

######################################################
### Read in the map and process it for RF calculation
pStruct                                               = proshade.ProSHADE_data ( )
pStruct.readInStructure                               ( gStruct, 0, pSet )

######################################################
### Write internal map and PDB to see if all worked fine
pStruct.writeMap                                      ( "gemmiTest.map" );
pStruct.writeGemmi                                    ( "gemmiTest.pdb", gStruct )

######################################################
### Compute distance to itself to make sure it works
pStruct.processInternalMap                            ( pSet )
pStruct.mapToSpheres                                  ( pSet )
pStruct.computeSphericalHarmonics                     ( pSet )
energyLevelsDescriptor                                = proshade.computeEnergyLevelsDescriptor     ( pStruct, pStruct, pSet )
traceSigmaDescriptor                                  = proshade.computeTraceSigmaDescriptor       ( pStruct, pStruct, pSet )
fullRotationFunctionDescriptor                        = proshade.computeRotationFunctionDescriptor ( pStruct, pStruct, pSet )
print                                                 ( "The energy levels distance is          %+1.3f" % ( energyLevelsDescriptor ) )
print                                                 ( "The trace sigma distance is            %+1.3f" % ( traceSigmaDescriptor ) )
print                                                 ( "The rotation function distance is      %+1.3f" % ( fullRotationFunctionDescriptor ) )

### Expected output: TThe energy levels distance is          +1.000
### Expected output: TThe trace sigma distance is            +1.000
### Expected output: TThe rotation function distance is      +1.000

######################################################
### Release ProSHADE memory
del gStruct
del pStruct
del pSet

######################################################
### Done
