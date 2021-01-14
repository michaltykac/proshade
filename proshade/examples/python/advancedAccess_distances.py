######################################################
######################################################
#   \file advancedAccess_distances.py
#   \brief This code demonstrates the usage of the ProSHADE tool in the advanced mode for the shape distances mode.
#
#   This file shows a fast demonstration of how the advanced access interface can be used to compute the
#   distances between two (or more) structures and how the results can be obtained. This file does not
#   contain all the explanations and possible settings, for complete documentation, please see the
#   directAccess.py file instead.
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
pSet.task                                             = proshade.Distances
pSet.verbose                                          = -1
pSet.setResolution                                    ( 6.0 )

######################################################
### Create the structure objects
pStruct1                                              = proshade.ProSHADE_data ( pSet )
pStruct2                                              = proshade.ProSHADE_data ( pSet )

######################################################
### Read in the structures
pStruct1.readInStructure                              ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/bf/1BFO_A_dom_1.pdb", 0, pSet ) # A BALBES domain 1BFO_A_dom_1
pStruct2.readInStructure                              ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/h8/1H8N_A_dom_1.pdb", 1, pSet ) # A BALBES domain 1H8N_A_dom_1

######################################################
### Process maps
pStruct1.processInternalMap                           ( pSet )
pStruct2.processInternalMap                           ( pSet )

######################################################
### Map to spheres
pStruct1.mapToSpheres                                 ( pSet )
pStruct2.mapToSpheres                                 ( pSet )

######################################################
### Compute spherical harmonics
pStruct1.computeSphericalHarmonics                    ( pSet )
pStruct2.computeSphericalHarmonics                    ( pSet )

######################################################
### Get the distances
energyLevelsDescriptor                                = proshade.computeEnergyLevelsDescriptor    ( pStruct1, pStruct2, pSet )
traceSigmaDescriptor                                  = proshade.computeTraceSigmaDescriptor      ( pStruct1, pStruct2, pSet )
fullRotationFunctionDescriptor                        = proshade.computeRotationunctionDescriptor ( pStruct1, pStruct2, pSet )

######################################################
### Print results
print                                                 ( "The energy levels distance is          %+1.3f" % ( energyLevelsDescriptor ) )
print                                                 ( "The trace sigma distance is            %+1.3f" % ( traceSigmaDescriptor ) )
print                                                 ( "The rotation function distance is      %+1.3f" % ( fullRotationFunctionDescriptor ) )

######################################################
### Expected output
#   The energy levels distance is          +0.895
#   The trace sigma distance is            +0.960
#   The rotation function distance is      +0.756

######################################################
### Release C++ pointers
del pStruct1
del pStruct2
del pSet
