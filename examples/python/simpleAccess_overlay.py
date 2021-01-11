######################################################
######################################################
#   \file simpleAccess_overlay.py
#   \brief This code demonstrates the usage of the ProSHADE tool in the simple mode for the structure overlay functionality.
#
#   This code demonstrates how the ProSHADE python module can be loaded,
#   the settings object created and filled (including all optional values),
#   supplied to the run object and how the optimal overlay results (i.e.)
#   the rotation matrix, the rotation centre and the rotation centre to
#   optimal overlay translation can be obained after the run.
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
import numpy

######################################################
### Import ProSHADE (assumes installation through pip)
import proshade

######################################################
### Create the settings object
ps                                                    = proshade.ProSHADE_settings ( )

######################################################
### Set up the run
ps.task                                               = proshade.OverlayMap
ps.verbose                                            = -1                                                               # How verbose should the run be? -1 Means no verbal output at all.
ps.setResolution                                      ( 4.0 )                                                            # The resolution to which the calculations will be done. NOTE: Not necessarily the resolution of the structure!
ps.addStructure                                       ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/bf/1BFO_A_dom_1.pdb" ) # A path to the structure to be processed. This is a BALBES domain 1BFO_A_dom_1.
ps.addStructure                                       ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/h8/1H8N_A_dom_1.pdb" ) # A path to the structure to be processed. This is a BALBES domain 1H8N_A_dom_1.

######################################################
### Set up the run Further useful settings
ps.forceP1                                            = True                                                             # Should PDB files be forced to have P1 spacegroup?
ps.removeWaters                                       = True                                                             # Should PDB files have their water molecules removed?
ps.firstModelOnly                                     = True                                                             # Should PDB files have only their first model used, or should ProSHADE use all models?
ps.setProgressiveSphereMapping                        ( False )                                                          # Should smaller spheres be less sampled? It is considerably faster, but may sacrifice some (little) accuracy.
ps.setOverlaySaveFile                                 ( "overlayResults" )                                               # Filename where the overlayed moving structure should be saved.
ps.setOverlayJsonFile                                 ( "movedStructureOperations.json" )                                # Filename where the overlay operations should be saved.
ps.setMasking                                         ( False )                                                          # Should maps be masked by blurring?
ps.setMapReboxing                                     ( False )                                                          # Should the structure be re-boxed? Required masking to be done in order to be meaningful.
ps.setNormalisation                                   ( False )                                                          # Should internal map representation be normalised to mean 0 and standard deviation 1?
ps.setMapCentering                                    ( False )                                                          # Move structure COM to the centre of map box?
ps.setMapResolutionChange                             ( True )                                                           # Should maps be re-sample to the computation resolution using reciprocal-space re-sampling?
ps.setMapResolutionChangeTriLinear                    ( False )                                                          # Should maps be re-sample to the computation resolution using real-space tri-linear interpolation?
ps.setExtraSpace                                      ( 10.0 )                                                           # Extra space in Angs to be added when creating internap map representation. This helps avoid map effects from other cells.

######################################################
### All other (possibly other tasks related) settings
ps.setSymmetryRotFunPeaks                             ( True )                                                           # Should the new angle-axis space symmetry detection be used?
ps.setBicubicInterpolationSearch                      ( True )                                                           # Should bi-cubic interpolation between peak grid indices be done?
ps.setMaxSymmetryFold                                 ( 30 )                                                             # The maximum prime number fold that will be searched for.
ps.setPeakNeighboursNumber                            ( 1 )                                                              # Numer of points in each direction which needs to be lower in order for the central point to be considered a peak.
ps.setPeakNaiveNoIQR                                  ( -999.9 )                                                         # Peak searching threshold for too low peaks in number of inter-quartile ranges from median of the non-peak point values.
ps.setMissingPeakThreshold                            ( 0.3 )                                                            # Fraction of peaks that can be missing for missing axis search to be initiated.
ps.setAxisComparisonThreshold                         ( 0.1 )                                                            # The dot product difference within which two axes are considered the same.
ps.setMinimumPeakForAxis                              ( 0.3 )                                                            # The minimum peak height for axis to be used.
ps.setRequestedSymmetry                               ( "" )                                                             # Which symmetry type (C,D,T,O or I) is requested to be detected? If none, then leave empty
ps.setRequestedFold                                   ( 0 )                                                              # For C and D symmetries, which symmetry fold is requested to be detected? If none, leave 0.
ps.setMapInversion                                    ( False )                                                          # Should all map positions x,y,z be swapped to -x,-y,-z? Use this only if your helices have the wrong hand as a result of first runs of map computation.
ps.setBandwidth                                       ( 0 )                                                              # The spherical harmonics bandwidth to which to compute. Set to 0 for automatic determination.
ps.setPhaseUsage                                      ( True )                                                           # Use full maps, or Patterson-like maps?
ps.setSphereDistances                                 ( 0.0 )                                                            # The distance between spheres. Use 0.0 for automatic determination.
ps.setIntegrationOrder                                ( 0 )                                                              # The order of the Gauss-Legendre integration computation. Set to 0 for automatic determination.
ps.setTaylorSeriesCap                                 ( 10 )                                                             # Set the Taylor series approximation cap. 10 seems like a fast and accurate value, but feel free to change.
ps.setEnergyLevelsComputation                         ( True )                                                           # Should energy levels descriptor be computed, assuming Distances are required (irrelevant otherwise)?
ps.setTraceSigmaComputation                           ( True )                                                           # Should trace sigma descriptor be computed, assuming Distances are required (irrelevant otherwise)?
ps.setRotationFunctionComputation                     ( True )                                                           # Should rotation function descriptor be computed, assuming Distances are required (irrelevant otherwise)?
ps.setEnLevShellWeight                                ( 1.0 )                                                            # The weighting of shell distances for energy levels descriptor.
ps.setOutputFilename                                  ( "reBoxed" )                                                      # Filename to where re-boxed structure will be written to.
ps.setBoundsSpace                                     ( 3.0 )                                                            # The extra space in Angs to add to the minimal boundaries when re-boxing.
ps.setBoundsThreshold                                 ( 0 )                                                              # If two boundaries are within this threshold, the smaller one will be increased to have the same value as the larger one.
ps.setSameBoundaries                                  ( False )                                                          # Make multiple structures have the same boundaries. This is useful for half-maps.
ps.setPDBBFactor                                      ( -1.0 )                                                           # Should all B-factors in a PDB file changed to this value? If no, set to negative value.
ps.setMaskBlurFactor                                  ( 350.0 )                                                          # If masking, what blur factor should be used? 350 seems to work for most maps.
ps.setMaskIQR                                         ( 3.0 )                                                            # Number of inter-quartile ranges from median to use as the masking threshold.
ps.setMaskSaving                                      ( False )                                                          # Should map mask be saved?
ps.setMaskFilename                                    ( "maskFile" )                                                     # The filename (no extension) to which the map masks will be saved into.

######################################################
### Run ProSHADE
rn                                                    = proshade.ProSHADE_run ( ps )

######################################################
### Retrieve results
eulerAngles                                           = rn.getEulerAngles                ( )
rotationMatrix                                        = rn.getOptimalRotMat              ( )
translationToOrigin                                   = rn.getTranslationToOrigin        ( )
translationToOverlay                                  = rn.getOriginToOverlayTranslation ( )

######################################################
### Print results
print                                                 ( "Optimal rotation Euler angles : " + str( eulerAngles[0] ) + "\t" + str( eulerAngles[1] ) + "\t" + str( eulerAngles[2] ) )
print                                                 ( "Optimal rotation matrix       : " + str( rotationMatrix[0][0] ) + "\t" + str( rotationMatrix[0][1] ) + "\t" + str( rotationMatrix[0][2] ) )
print                                                 ( "                              : " + str( rotationMatrix[1][0] ) + "\t" + str( rotationMatrix[1][1] ) + "\t" + str( rotationMatrix[1][2] ) )
print                                                 ( "                              : " + str( rotationMatrix[2][0] ) + "\t" + str( rotationMatrix[2][1] ) + "\t" + str( rotationMatrix[2][2] ) )
print                                                 ( "Translation to origin         : " + str( translationToOrigin[0] ) + "\t" + str( translationToOrigin[1] ) + "\t" + str( translationToOrigin[2] ) )
print                                                 ( "Translation to overlay        : " + str( translationToOverlay[0] ) + "\t" + str( translationToOverlay[1] ) + "\t" + str( translationToOverlay[2] ) )

######################################################
### Expected outuput
#   Optimal rotation Euler angles : 5.4325056    0.7526409    3.9270065
#   Optimal rotation matrix       : -0.87191415    -0.07835774    0.4833485
#                                 : -0.19118044    -0.85428894    -0.48336366
#                                 : 0.45079455    -0.5138584    0.7298862
#   Translation to origin         : -16.0    -20.0    -24.0
#   Translation to overlay        : 8.0    8.0    -6.0

######################################################
### Release memory
del rn
del ps

######################################################
### Done
