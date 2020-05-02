##############################################
##############################################
# \file simpleAccess_reboxing.py
# \brief This code demonstrates the usage of the ProSHADE tool in the simple mode for the reboxing functionality.
#
# This code demonstrates how the ProSHADE python module can be used in the simple access mode to re-box structures.
# It shows how the module can be loaded from non-system directory  (assuming local installation), how the settings
# object can be created and filled with values, how proshade run can be done and how the results can be programmatically
# accessed, including access to the map data in Numpy array format.
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

##############################################
### Import system modules
import sys
import numpy

### Import ProSHADE from non-system folder (local installation assumed)
#sys.path.append                               ( "/Users/mysak/BioCEV/proshade/experimental/install/python3" )
import proshade

##############################################
### Create the ProSHADE_settings object
pSet                                          = proshade.ProSHADE_settings ()

##############################################
### Set the settings for Symmetry detection defaults

### Required values
pSet.task                                     = proshade.MapManip                      # The task ProSHADE is expected to perform
pSet.verbose                                  = -1                                     # How verbose should the run be? Use -1 for absolute silence.
pSet.setResolution                            ( 10.0 )                                 # The resolution to which computations are to be done. May be lower or higher than the experimentally measured resolution.
pSet.addStructure                             ( "/Users/mysak/LMB/1_ProteinDomains/14_FullRotFunction/14_AdditionalSymmetries/TET1.map" )         # The path to the structure to be processed.
pSet.addStructure                             ( "/Users/mysak/LMB/proshade/exp/demo/testMap.map" )                                                # The path to the structure to be processed.
pSet.addStructure                             ( "/Users/mysak/LMB/proshade/exp/demo/testMap2.map" )                                               # The path to the structure to be processed.


### Useful settings
pSet.setMasking                               ( True )                                 # Should maps be masked by blurring?
pSet.setMapReboxing                           ( True )                                 # Should the structure be re-boxed? Required masking to be done in order to be meaningful.
pSet.setNormalisation                         ( False )                                # Should internal map representation be normalised to mean 0 and standard deviation 1?
pSet.setMapCentering                          ( False )                                # Move structure COM to the centre of map box?
pSet.setOutputFilename                        ( str ( "reBoxed" ) )                    # Filename to where re-boxed structure will be written to.
pSet.setMapResolutionChange                   ( False )                                # Should maps be re-sample to the computation resolution?
pSet.setBoundsSpace                           ( 3.0 )                                  # The extra space in Angs to add to the minimal boundaries when re-boxing.
pSet.setBoundsThreshold                       ( 0 )                                    # If two boundaries are within this threshold, the smaller one will be increased to have the same value as the larger one.
pSet.setSameBoundaries                        ( False )                                # Make multiple structures have the same boundaries. This is useful for half-maps.
pSet.setPDBBFactor                            ( -1.0 )                                 # Should all B-factors in a PDB file changed to this value? If no, set to negative value.
pSet.setMaskBlurFactor                        ( 350.0 )                                # If masking, what blur factor should be used? 350 seems to work for most maps.
pSet.setMaskIQR                               ( 3.0 )                                  # Number of inter-quartile ranges from median to use as the masking threshold.
pSet.setMaskSaving                            ( False )                                # Should map mask be saved?
pSet.setMaskFilename                          ( str ( "maskFile" ) )                   # The filename (no extension) to which the map masks will be saved into.


### All other (possible other tasks related) settings
pSet.setExtraSpace                            ( 10.0 )                                 # Extra space in Angs to be added when creating internap map representation. This helps avoid map effects from other cells.
pSet.setPeakNeighboursNumber                  ( 1 )                                    # Numer of points in each direction which needs to be lower in order for the central point to be considered a peak.
pSet.setPeakNaiveNoIQR                        ( 5.0 )                                  # Peak searching threshold for too low peaks in number of inter-quartile ranges from median of the non-peak point values.
pSet.setMissingPeakThreshold                  ( 0.3 )                                  # Fraction of peaks that can be missing for missing axis search to be initiated.
pSet.setAxisComparisonThreshold               ( 0.05 )                                 # The dot product difference within which two axes are considered the same.
pSet.setRequestedSymmetry                     ( "" )                                   # Which symmetry type (C,D,T,O or I) is requested to be detected? If none, then leave empty
pSet.setRequestedFold                         ( 0 )                                    # For C and D symmetries, which symmetry fold is requested to be detected? If none, leave 0.
pSet.setOverlaySaveFile                       ( "moved" )                              # Filename where the overlayed moving structure should be saved.
pSet.setMapInversion                          ( False )                                # Should all map positions x,y,z be swapped to -x,-y,-z? Use this only if your helices have the wrong hand ...
                                                                                       # ... as a result of first runs of map computation.
pSet.setBandwidth                             ( 0 )                                    # The spherical harmonics bandwidth to which to compute. Set to 0 for automatic determination.
pSet.setAngularResolution                     ( 0 )                                    # The resolution of the sphere mapping. Set to 0 for automatic determination.
pSet.setPhaseUsage                            ( True )                                 # Use full maps, or Patterson-like maps?
pSet.setSphereDistances                       ( 0.0 )                                  # The distance between spheres. Use 0.0 for automatic determination.
pSet.setIntegrationOrder                      ( 0 )                                    # The order of the Gauss-Legendre integration computation. Set to 0 for automatic determination.
pSet.setTaylorSeriesCap                       ( 10 )                                   # Set the Taylor series approximation cap. 10 seems like a fast and accurate value, but feel free to change.
pSet.setProgressiveSphereMapping              ( False )                                # Should smaller spheres be less sampled? It is considerably faster, but may sacrifice some (little) accuracy.
pSet.setEnergyLevelsComputation               ( True )                                 # Should energy levels descriptor be computed, assuming Distances are required (irrelevant otherwise)?
pSet.setTraceSigmaComputation                 ( True )                                 # Should trace sigma descriptor be computed, assuming Distances are required (irrelevant otherwise)?
pSet.setRotationFunctionComputation           ( True )                                 # Should rotation function descriptor be computed, assuming Distances are required (irrelevant otherwise)?
pSet.setEnLevShellWeight                      ( 1.0 )                                  # The weighting of shell distances for energy levels descriptor.


##############################################
### Print the settings for reference
#pSet.printSettings                            ( )

##############################################
### Run the Distances task
pRun                                          = proshade.ProSHADE_run ( pSet )

##############################################
### Get the reboxed information and maps
originalBounds                                = proshade.getOrigBounds  ( pRun )
reboxedBounds                                 = proshade.getReboxBounds ( pRun )
reboxedMaps                                   = proshade.getReboxMap    ( pRun )

##############################################
### Delete the C++ pointers
del pRun
del pSet

##############################################
### Print the results
for iter in range ( 0, len ( originalBounds ) ):
    print ( "Structure " + str( iter ) + " had bounds " + str( originalBounds[iter][0] ) + " - " + str( originalBounds[iter][1] ) + " ; " + str( originalBounds[iter][2] ) + " - " + str( originalBounds[iter][3] ) + " ; " + str( originalBounds[iter][4] ) + " - " + str( originalBounds[iter][5] ) + " and is now re-boxed to " + str( reboxedBounds[iter][0] ) + " - " + str( reboxedBounds[iter][1] ) + " ; " + str( reboxedBounds[iter][2] ) + " - " + str( reboxedBounds[iter][3] ) + " ; " + str( reboxedBounds[iter][4] ) + " - " + str( reboxedBounds[iter][5] ) )
    print ( "  This structure's re-boxed map has been read into Python as " + str( type( reboxedMaps[iter] ) ) + " with " + str( len ( reboxedMaps[iter] ) ) + " points having mean of " + str( numpy.mean ( reboxedMaps[iter] ) ) )

##############################################
### Done
