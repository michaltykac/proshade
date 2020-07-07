##############################################
##############################################
#   \file simpleAccess_reboxing.py
#   \brief This code demonstrates the usage of the ProSHADE tool in the simple mode for the reboxing functionality.
#
#   This code demonstrates how the ProSHADE python module can be loaded, the settings object created and filled (including all
#   optional values), supplied to the run object and how the optimal rotation matrix, Euler angles and optimal translation vector
#   can be obtained using the simple access interface.
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
#   \date      JUL 2020
##############################################
##############################################

##############################################
### Import system modules
import sys
import numpy

### Import ProSHADE from non-system folder (local installation assumed)
#sys.path.append                               ( "/Users/mysak/BioCEV/proshade/master/install/python2" )
import proshade

##############################################
### Create the ProSHADE_settings object
pSet                                          = proshade.ProSHADE_settings ()

##############################################
### Set the settings for Symmetry detection defaults

### Required values
pSet.task                                     = proshade.OverlayMap                    # The task ProSHADE is expected to perform
pSet.verbose                                  = -1                                     # How verbose should the run be? Use -1 for absolute silence.
pSet.setResolution                            ( 4.0 )                                  # The resolution to which computations are to be done. May be lower or higher than the experimentally measured resolution.
pSet.addStructure                             ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/bf/1BFO_A_dom_1.pdb" ) # The path to the structure to be processed. This is a BALBES domain 1BFO_A_dom_1.
pSet.addStructure                             ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/h8/1H8N_A_dom_1.pdb" ) # The path to the structure to be processed. This is a BALBES domain 1H8N_A_dom_1.


### Useful settings
pSet.setProgressiveSphereMapping              ( False )                                # Should smaller spheres be less sampled? It is considerably faster, but may sacrifice some (little) accuracy.
pSet.setOverlaySaveFile                       ( "moved" )                              # Filename where the overlayed moving structure should be saved.
pSet.setMasking                               ( False )                                # Should maps be masked by blurring?
pSet.setMapReboxing                           ( False )                                # Should the structure be re-boxed? Required masking to be done in order to be meaningful.
pSet.setNormalisation                         ( False )                                # Should internal map representation be normalised to mean 0 and standard deviation 1?
pSet.setMapCentering                          ( False )                                # Move structure COM to the centre of map box?
pSet.setMapResolutionChange                   ( False )                                # Should maps be re-sample to the computation resolution?

### All other (possible other tasks related) settings
pSet.setExtraSpace                            ( 10.0 )                                 # Extra space in Angs to be added when creating internap map representation. This helps avoid map effects from other cells.
pSet.setPeakNeighboursNumber                  ( 1 )                                    # Numer of points in each direction which needs to be lower in order for the central point to be considered a peak.
pSet.setPeakNaiveNoIQR                        ( 5.0 )                                  # Peak searching threshold for too low peaks in number of inter-quartile ranges from median of the non-peak point values.
pSet.setMissingPeakThreshold                  ( 0.3 )                                  # Fraction of peaks that can be missing for missing axis search to be initiated.
pSet.setAxisComparisonThreshold               ( 0.05 )                                 # The dot product difference within which two axes are considered the same.
pSet.setRequestedSymmetry                     ( "" )                                   # Which symmetry type (C,D,T,O or I) is requested to be detected? If none, then leave empty
pSet.setRequestedFold                         ( 0 )                                    # For C and D symmetries, which symmetry fold is requested to be detected? If none, leave 0.
pSet.setMapInversion                          ( False )                                # Should all map positions x,y,z be swapped to -x,-y,-z? Use this only if your helices have the wrong hand ...
                                                                                       # ... as a result of first runs of map computation.
pSet.setBandwidth                             ( 0 )                                    # The spherical harmonics bandwidth to which to compute. Set to 0 for automatic determination.
pSet.setPhaseUsage                            ( True )                                 # Use full maps, or Patterson-like maps?
pSet.setSphereDistances                       ( 0.0 )                                  # The distance between spheres. Use 0.0 for automatic determination.
pSet.setIntegrationOrder                      ( 0 )                                    # The order of the Gauss-Legendre integration computation. Set to 0 for automatic determination.
pSet.setTaylorSeriesCap                       ( 10 )                                   # Set the Taylor series approximation cap. 10 seems like a fast and accurate value, but feel free to change.
pSet.setEnergyLevelsComputation               ( True )                                 # Should energy levels descriptor be computed, assuming Distances are required (irrelevant otherwise)?
pSet.setTraceSigmaComputation                 ( True )                                 # Should trace sigma descriptor be computed, assuming Distances are required (irrelevant otherwise)?
pSet.setRotationFunctionComputation           ( True )                                 # Should rotation function descriptor be computed, assuming Distances are required (irrelevant otherwise)?
pSet.setEnLevShellWeight                      ( 1.0 )                                  # The weighting of shell distances for energy levels descriptor.
pSet.setOutputFilename                        ( str ( "reBoxed" ) )                    # Filename to where re-boxed structure will be written to.
pSet.setBoundsSpace                           ( 3.0 )                                  # The extra space in Angs to add to the minimal boundaries when re-boxing.
pSet.setBoundsThreshold                       ( 0 )                                    # If two boundaries are within this threshold, the smaller one will be increased to have the same value as the larger one.
pSet.setSameBoundaries                        ( False )                                # Make multiple structures have the same boundaries. This is useful for half-maps.
pSet.setPDBBFactor                            ( -1.0 )                                 # Should all B-factors in a PDB file changed to this value? If no, set to negative value.
pSet.setMaskBlurFactor                        ( 350.0 )                                # If masking, what blur factor should be used? 350 seems to work for most maps.
pSet.setMaskIQR                               ( 3.0 )                                  # Number of inter-quartile ranges from median to use as the masking threshold.
pSet.setMaskSaving                            ( False )                                # Should map mask be saved?
pSet.setMaskFilename                          ( str ( "maskFile" ) )                   # The filename (no extension) to which the map masks will be saved into.




##############################################
### Print the settings for reference
#pSet.printSettings                            ( )

##############################################
### Run the Distances task
pRun                                          = proshade.ProSHADE_run ( pSet )

##############################################
### Get the reboxed information and maps
eulerAngles                                   = proshade.getEulerAngles ( pRun )
rotMatrix                                     = proshade.getRotationMat ( pRun )
translation                                   = proshade.getTranslation ( pRun )

##############################################
### Delete the C++ pointers
del pRun
del pSet

##############################################
### Print the results
print ( "Optimal overlay rotation (Euler angles) : %+2.4f   %+2.4f   %+2.4f" % ( eulerAngles[0] , eulerAngles[1] , eulerAngles[2] ) )
print ( "Optimal overlay rotation (rot. matrix)  : %+2.4f   %+2.4f   %+2.4f" % ( rotMatrix[0][0], rotMatrix[0][1], rotMatrix[0][2] ) )
print ( "                                        : %+2.4f   %+2.4f   %+2.4f" % ( rotMatrix[1][0], rotMatrix[1][1], rotMatrix[1][2] ) )
print ( "                                        : %+2.4f   %+2.4f   %+2.4f" % ( rotMatrix[2][0], rotMatrix[2][1], rotMatrix[2][2] ) )
print ( "Optimal overlay translation (Angstroms) : %+2.4f   %+2.4f   %+2.4f" % ( translation[0] , translation[1] , translation[2] ) )

### Expected output
#  Optimal overlay rotation (Euler angles) : +5.5084   +0.7639   +3.8732
#  Optimal overlay rotation (rot. matrix)  : -0.8513   -0.1015   +0.5147
#                                          : -0.1758   -0.8692   -0.4621
#                                          : +0.4943   -0.4839   +0.7221
#  Optimal overlay translation (Angstroms) : +8.0000   +6.0000   -6.0000

##############################################
### Done
