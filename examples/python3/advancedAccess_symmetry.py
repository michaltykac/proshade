##############################################
##############################################
# \file advancedAccess_symmetry.py
# \brief This code demonstrates the usage of the ProSHADE tool in the advanced mode.
#
# This file shows a fast demonstration of how the advanced access interfacte can be used to compute the
# symmetry of a particular structure and how the results can be obtained. This file does not contain all
# the explanations and possible settings, for complete documentation, please see the advancedAccess.py
# file instead.
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
# \version   0.7.1
# \date      NOV 2019
##############################################
##############################################


### System modules
import sys
import numpy

### Import ProSHADE from non-system folder (local installation assumed)
sys.path.append                               ( "/Users/mysak/BioCEV/proshade/master/install/python3" )
import proshade

### Create the settings object
pSet                                          = proshade.ProSHADE_settings ()

### Set settings values
pSet.task                                     = proshade.Symmetry
pSet.verbose                                  = 1
pSet.setResolution                            ( 10.0 )
pSet.moveToCOM                                = True

### Create the structure object
pStruct                                       = proshade.ProSHADE_data ( pSet )

### Read in the structure
pStruct.readInStructure                       ( "/Users/mysak/LMB/proshade/exp/demo/C3.pdb", 0, pSet )

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

### Release C++ pointers
del pStruct
del pSet

### Done
