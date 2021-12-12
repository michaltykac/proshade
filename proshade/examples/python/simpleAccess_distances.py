######################################################
######################################################
#   \file simpleAccess_distances.py
#   \brief This code demonstrates the usage of the ProSHADE tool in the simple mode for distances functionality.
#
#   This code demonstrates how the ProSHADE python module can be loaded,
#   the settings object created and filled (including all optional values),
#   supplied to the run object and how to three descriptor values can be
#   obtained programatically.
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
import numpy

######################################################
### Import ProSHADE (assumes installation through pip)
import proshade

######################################################
### Create the settings object
ps                                                    = proshade.ProSHADE_settings ( )

######################################################
### Set up the run
ps.task                                               = proshade.Distances
ps.verbose                                            = -1                                                               # How verbose should the run be? -1 Means no verbal output at all.
ps.setResolution                                      ( 6.0 )                                                            # The resolution to which the calculations will be done. NOTE: Not necessarily the resolution of the structure!
ps.addStructure                                       ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/bf/1BFO_A_dom_1.pdb" ) # A BALBES domain 1BFO_A_dom_1
ps.addStructure                                       ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/h8/1H8N_A_dom_1.pdb" ) # A BALBES domain 1H8N_A_dom_1
ps.addStructure                                       ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/ig/3IGU_A_dom_1.pdb" ) # A BALBES domain 3IGU_A_dom_1

######################################################
### Further useful settings
ps.forceP1                                            = True                                                             # Should PDB files be forced to have P1 spacegroup?
ps.setNegativeDensity                                 ( True )                                                           # Should the negative density be removed from input files?
ps.removeWaters                                       = True                                                             # Should PDB files have their water molecules removed?
ps.firstModelOnly                                     = True                                                             # Should PDB files have only their first model used, or should ProSHADE use all models?
ps.setProgressiveSphereMapping                        ( False )                                                          # Should smaller spheres be less sampled? It is considerably faster, but may sacrifice some (little) accuracy.
ps.setMapResolutionChange                             ( False )                                                          # Should maps be re-sample to the computation resolution using reciprocal-space re-sampling?
ps.setMapResolutionChangeTriLinear                    ( False )                                                          # Should maps be re-sample to the computation resolution using real-space tri-linear interpolation?
ps.setPDBBFactor                                      ( -1.0 )                                                           # Should all B-factors in a PDB file changed to this value? If no, set to negative value.
ps.setBandwidth                                       ( 0 )                                                              # The spherical harmonics bandwidth to which to compute. Set to 0 for automatic determination.
ps.setPhaseUsage                                      ( True )                                                           # Use full maps, or Patterson-like maps?
ps.setSphereDistances                                 ( 0.0 )                                                            # The distance between spheres. Use 0.0 for automatic determination.
ps.setIntegrationOrder                                ( 0 )                                                              # The order of the Gauss-Legendre integration computation. Set to 0 for automatic determination.
ps.setTaylorSeriesCap                                 ( 10 )                                                             # Set the Taylor series approximation cap. 10 seems like a fast and accurate value, but feel free to change.
ps.setNormalisation                                   ( False )                                                          # Should internal map representation be normalised to mean 0 and standard deviation 1?
ps.setMapInversion                                    ( False )                                                          # Should all map positions x,y,z be swapped to -x,-y,-z? Use this only if your helices have the wrong hand as a result of first runs of map computation.
ps.setMasking                                         ( False )                                                          # Should maps be masked by blurring?
ps.setMapReboxing                                     ( False )                                                          # Should the structure be re-boxed? Required masking to be done in order to be meaningful.
ps.setMapCentering                                    ( False )                                                          # Move structure COM to the centre of map box?
ps.setEnergyLevelsComputation                         ( True )                                                           # Should energy levels descriptor be computed, assuming Distances are required (irrelevant otherwise)?
ps.setTraceSigmaComputation                           ( True )                                                           # Should trace sigma descriptor be computed, assuming Distances are required (irrelevant otherwise)?
ps.setRotationFunctionComputation                     ( True )                                                           # Should rotation function descriptor be computed, assuming Distances are required (irrelevant otherwise)?
ps.setEnLevShellWeight                                ( 1.0 )                                                            # The weighting of shell distances for energy levels descriptor.


######################################################
### All other (possibly other tasks related) settings
ps.setSymmetryCentreSearch                            ( False )                                                          # Should symmetry centre be searched for? Takes a lot of time...
ps.setBicubicInterpolationSearch                      ( True )                                                           # Should bi-cubic interpolation between peak grid indices be done?
ps.setMaxSymmetryFold                                 ( 30 )                                                             # The maximum prime number fold that will be searched for.
ps.setFSCThreshold                                    ( 0.75 )                                                           # The threshold for FSC value under which the axis is considered to be likely noise.
ps.setPeakThreshold                                   ( 0.80 )                                                           # The threshold for peak height above which axes are considered possible.
ps.setMaskBlurFactor                                  ( 350.0 )                                                          # If masking, what blur factor should be used? 350 seems to work for most maps.
ps.setMaskIQR                                         ( 3.0 )                                                            # Number of inter-quartile ranges from median to use as the masking threshold.
ps.setMaskSaving                                      ( False )                                                          # Should map mask be saved?
ps.setMaskFilename                                    ( "maskFile" )                                                     # The filename (no extension) to which the map masks will be saved into.
ps.setAppliedMaskFilename                             ( "" )                                                             # The filename from which mask data will be read from.
ps.setBoundsSpace                                     ( 3.0 )                                                            # The extra space in Angs to add to the minimal boundaries when re-boxing.
ps.setBoundsThreshold                                 ( 0 )                                                              # If two boundaries are within this threshold, the smaller one will be increased to have the same value as the larger one.
ps.setSameBoundaries                                  ( False )                                                          # Make multiple structures have the same boundaries. This is useful for half-maps.
ps.setOutputFilename                                  ( "reBoxed" )                                                      # Filename to where re-boxed structure will be written to.
ps.setExtraSpace                                      ( 10.0 )                                                           # Extra space in Angs to be added when creating internap map representation. This helps avoid map effects from other cells.
ps.setPeakNeighboursNumber                            ( 1 )                                                              # Numer of points in each direction which needs to be lower in order for the central point to be considered a peak.
ps.setPeakNaiveNoIQR                                  ( -999.9 )                                                         # Peak searching threshold for too low peaks in number of inter-quartile ranges from median of the non-peak point values.
ps.setMissingPeakThreshold                            ( 0.3 )                                                            # Fraction of peaks that can be missing for missing axis search to be initiated.
ps.setAxisComparisonThreshold                         ( 0.1 )                                                            # The dot product difference within which two axes are considered the same.
ps.setMinimumPeakForAxis                              ( 0.3 )                                                            # The minimum peak height for axis to be used.
ps.setRequestedSymmetry                               ( "" )                                                             # Which symmetry type (C,D,T,O or I) is requested to be detected? If none, then leave empty
ps.setRequestedFold                                   ( 0 )                                                              # For C and D symmetries, which symmetry fold is requested to be detected? If none, leave 0.
ps.setOverlaySaveFile                                 ( "moved" )                                                        # Filename where the overlayed moving structure should be saved.
ps.setOverlayJsonFile                                 ( "movedStructureOperations.json" )                                # Filename where the overlay operations should be saved.

######################################################
### Run ProSHADE
rn                                                    = proshade.ProSHADE_run ( ps )

######################################################
### Retrieve results
energyLevelsDesc                                      = rn.getEnergyLevelsVector()
traceSigmaDesc                                        = rn.getTraceSigmaVector()
rotationFunctionDesc                                  = rn.getRotationFunctionVector()

######################################################
### Write out results
print                                                 ( "Energy levels distance          : " + str ( energyLevelsDesc[0] )     + "\t\t" + str ( energyLevelsDesc[1]     ) )
print                                                 ( "Trace sigma distance            : " + str ( traceSigmaDesc[0] )       + "\t\t" + str ( traceSigmaDesc[1]       ) )
print                                                 ( "Full rotation function distance : " + str ( rotationFunctionDesc[0] ) + "\t\t" + str ( rotationFunctionDesc[1] ) )

######################################################
### Expected output
#   Energy levels distance          : 0.8556063        0.5606434
#   Trace sigma distance            : 0.96467125        0.7530951
#   Full rotation function distance : 0.62452674        0.46518552

######################################################
### Release memory
del rn
del ps

######################################################
### Done
