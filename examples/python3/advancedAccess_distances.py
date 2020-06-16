##############################################
##############################################
# \file advancedAccess_distances.py
# \brief This code demonstrates the usage of the ProSHADE tool in the advanced mode.
#
# This file shows a fast demonstration of how the advanced access interfacte can be used to compute the
# distances between two (or more) structures and how the results can be obtained. This file does not
# contain all the explanations and possible settings, for complete documentation, please see the
# advancedAccess.py file instead.
#
# Please note that the input structure paths are hard coded and you will need
# to change these before this file runs.
#
# This file is part of the ProSHADE library for calculating
# shape descriptors and symmetry operators of protein structures.
# This is a prototype code, which is by no means complete or fully
# tested. Its use is at your own risk only. There is no quarantee
# that the results are correct.
#
# \author    Michal Tykac
# \author    Garib N. Murshudov
# \version   0.7.3
# \date      MAY 2020
##############################################
##############################################


### System modules
import sys
import numpy

### Import ProSHADE from non-system folder (local installation assumed)
#sys.path.append                               ( "/Users/mysak/BioCEV/proshade/experimental/install/python3" )
import proshade

### Create the settings object
pSet                                          = proshade.ProSHADE_settings ()

### Set settings values
pSet.task                                     = proshade.Distances
pSet.verbose                                  = 1
pSet.setResolution                            ( 6.0 )

### Create the structure objects
pStruct1                                      = proshade.ProSHADE_data ( pSet )
pStruct2                                      = proshade.ProSHADE_data ( pSet )

### Read in the structure
pStruct1.readInStructure                      ( "/Users/mysak/LMB/proshade/exp/demo/C2.pdb", 0, pSet )
pStruct2.readInStructure                      ( "/Users/mysak/LMB/proshade/exp/demo/testMap.map", 1, pSet )

### Process map
pStruct1.processInternalMap                   ( pSet )
pStruct2.processInternalMap                   ( pSet )

### Map to spheres
pStruct1.mapToSpheres                         ( pSet )
pStruct2.mapToSpheres                         ( pSet )

### Compute spherical harmonics
pStruct1.computeSphericalHarmonics            ( pSet )
pStruct2.computeSphericalHarmonics            ( pSet )

### Get the distances
energyLevelsDescriptor                        = proshade.computeEnergyLevelsDescriptor    ( pStruct1, pStruct2, pSet )
traceSigmaDescriptor                          = proshade.computeTraceSigmaDescriptor      ( pStruct1, pStruct2, pSet )
fullRotationFunctionDescriptor                = proshade.computeRotationunctionDescriptor ( pStruct1, pStruct2, pSet )

### Print results
print ( "The energy levels distance is          %+1.3f" % ( energyLevelsDescriptor ) )
print ( "The trace sigma distance is            %+1.3f" % ( traceSigmaDescriptor ) )
print ( "The rotation function distance is      %+1.3f" % ( fullRotationFunctionDescriptor ) )

### Release C++ pointers
del pStruct1
del pStruct2
del pSet

### Done
