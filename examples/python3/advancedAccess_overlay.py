##############################################
##############################################
# \file advancedAccess_overlay.py
# \brief This code demonstrates the usage of the ProSHADE tool in the advanced mode for overlay functionality.
#
# ...
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

### Create the object
pSet                                          = proshade.ProSHADE_settings ()

### Set settings values for optimal overlay
pSet.task                                     = proshade.OverlayMap
pSet.verbose                                  = 4
pSet.requestedResolution                      = 8.0;
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
pStruct_static.readInStructure                ( "/Users/mysak/BioCEV/proshade/00_GeneralTests/04_MapOverlay/test1.map", 0, pSet )
pStruct_moving.readInStructure                ( "/Users/mysak/BioCEV/proshade/00_GeneralTests/04_MapOverlay/test1_higherRotTrs.map", 1, pSet )

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
optimalRotationMatrix                         = proshade.getRotationMatrixFromEulerZXZ ( optimalRotationAngles )

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
pStruct_static.readInStructure                ( "/Users/mysak/BioCEV/proshade/00_GeneralTests/04_MapOverlay/test1.map", 0, pSet )
pStruct_moving.readInStructure                ( "/Users/mysak/BioCEV/proshade/00_GeneralTests/04_MapOverlay/test1_higherRotTrs.map", 1, pSet )

### Get spherical harmonics for moving structure only
pStruct_static.processInternalMap             ( pSet )

pStruct_moving.processInternalMap             ( pSet )
pStruct_moving.mapToSpheres                   ( pSet )
pStruct_moving.computeSphericalHarmonics      ( pSet )

### Rotate the moving structure (the signs with the Euler angles need to be like this!)
pStruct_moving.rotateMap                      ( pSet, -optimalRotationAngles[0], optimalRotationAngles[1], -optimalRotationAngles[2] )

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

### Print results
print ( "Optimal Euler angles         :  %+1.3f    %+1.3f    %+1.3f\n" % ( optimalRotationAngles[0], optimalRotationAngles[1], optimalRotationAngles[2] ) )
print ( "Optimal Euler rotation matrix:  %+1.3f    %+1.3f    %+1.3f" % ( optimalRotationMatrix[0][0], optimalRotationMatrix[0][1], optimalRotationMatrix[0][2] ) )
print ( "                             :  %+1.3f    %+1.3f    %+1.3f" % ( optimalRotationMatrix[1][0], optimalRotationMatrix[1][1], optimalRotationMatrix[1][2] ) )
print ( "                             :  %+1.3f    %+1.3f    %+1.3f\n" % ( optimalRotationMatrix[2][0], optimalRotationMatrix[2][1], optimalRotationMatrix[2][2] ) )
print ( "Optimal translation          :  %+1.3f    %+1.3f    %+1.3f" % ( optimalTranslationVector[0], optimalTranslationVector[1], optimalTranslationVector[2] ) )

### Write out the map - it has some artefacts, this is caused by the double interpolation - I recommend applying the rotation and translation in EMDA instead of using this map.
pStruct_moving.translateMap                   ( pSet, optimalTranslationVector[0], optimalTranslationVector[1], optimalTranslationVector[2] );
pStruct_moving.writeMap                       ( "/Users/mysak/Desktop/movedPy.map" )

### Delete the data
del pStruct_static
del pStruct_moving

### Done
