/*! \file advancedAccess_overlay.cpp
    \brief This file showcases how the ProSHADE structure overlay can be computed using the dynamic library's advanced interface.

    This code is an example of how any C++ project linking the ProSHADE library can access the advanced ProSHADE interface. More specifically, to compute the structure overlay, this example firstly reads in two structures and then removes the phase information
    from both of them. This results in the internal map representation being a Patterson map. Next, these two maps are mapped onto a set of spheres, have their spherical harmonics computed and finally a rotation function map is computed from the two structures'
    spherical harmonics decomposition. From the rottion function, the optimal rotation angles are obtained and the two structures are deleted. Then, the structures are read in again, this time keeping the phase information intact. The moving structure is then processed
    to have its spherical harmonics computed as to allow this structure to be rotated in the spherical harmonics space. Once the moving structure is rotated, the two structures are zero padded to have the same size and have their translation map computed. Finally, from
    this translation map, the optimal translation vector is computed and thus the optimal overlay is obtained.

    Copyright by Michal Tykac and individual contributors. All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.

    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.6.2
    \date      DEC 2021
*/

//==================================================== ProSHADE
#include "../../src/proshade/ProSHADE.hpp"

//==================================================== Main
int main ( int argc, char **argv )
{
    //================================================ Create the settings object and parse the command line arguments
    ProSHADE_Task task                                = OverlayMap;                          // Setting the task ahead sets most of the default settings to best values for the task.
    ProSHADE_settings* settings                       = new ProSHADE_settings ( task );      // Creating the ProSHADE_settings object, which caries all of the settings and where everything can be set.

    //================================================ Required settings
    settings->verbose                                 = -1;                                  // How verbose should the run be? -1 Means no verbal output at all.
    settings->setResolution                           ( 4.0 );                               // The resolution to which the calculations will be done. NOTE: Not necessarily the resolution of the structure!
    settings->setPhaseUsage                           ( false );                             // Use full maps, or Patterson-like maps?
    settings->setMapResolutionChange                  ( true );                              // Should maps be re-sample to the computation resolution using reciprocal-space re-sampling?
    settings->setMapResolutionChangeTriLinear         ( false );                             // Should maps be re-sample to the computation resolution using real-space tri-linear interpolation?
    settings->setMasking                              ( false );                             // Should maps be masked by blurring?
    settings->setMapCentering                         ( false );                             // Move structure COM to the centre of map box?
    settings->setMapReboxing                          ( false );                             // Should the structure be re-boxed? Required masking to be done in order to be meaningful.
    
    //================================================ Further useful settings
    settings->forceP1                                 = true;                                // Should PDB files be forced to have P1 spacegroup?
    settings->setNegativeDensity                      ( true );                              // Should the negative density be removed from input files?
    settings->removeWaters                            = true;                                // Should PDB files have their water molecules removed?
    settings->firstModelOnly                          = true;                                // Should PDB files have only their first model used, or should ProSHADE use all models?
    settings->setProgressiveSphereMapping             ( false );                             // Should smaller spheres be less sampled? It is considerably faster, but may sacrifice some (little) accuracy.
    settings->setOverlaySaveFile                      ( "overlayResuls" );                   // Filename where the overlayed moving structure should be saved.
    settings->setOverlayJsonFile                      ( "movedStructureOperations.json" );   // Filename where the overlay operations should be saved.
    settings->setNormalisation                        ( false );                             // Should internal map representation be normalised to mean 0 and standard deviation 1?
    settings->setExtraSpace                           ( 10.0 );                              // Extra space in Angs to be added when creating internap map representation. This helps avoid map effects from other cells.
    
    //================================================ All other (possibly other tasks related) settings
    settings->setSymmetryCentreSearch                 ( false );                             // Should symmetry centre be searched for? Takes a lot of time...
    settings->setBicubicInterpolationSearch           ( true );                              // Should bi-cubic interpolation between peak grid indices be done?
    settings->setMaxSymmetryFold                      ( 30 );                                // The maximum prime number fold that will be searched for.
    settings->setFSCThreshold                         ( 0.75 );                              // Sets the minimum FSC threshold for axis to be considered detected.
    settings->setPeakThreshold                        ( 0.80 );                              // Sets the minimum peak height threshold for axis to be considered possible.
    settings->setPeakNeighboursNumber                 ( 1 );                                 // Numer of points in each direction which needs to be lower in order for the central point to be considered a peak.
    settings->setPeakNaiveNoIQR                       ( -999.9 );                            // Peak searching threshold for too low peaks in number of inter-quartile ranges from median of the non-peak point values.
    settings->setMissingPeakThreshold                 ( 0.3 );                               // Fraction of peaks that can be missing for missing axis search to be initiated.
    settings->setAxisComparisonThreshold              ( 0.1 );                               // The dot product difference within which two axes are considered the same.
    settings->setMinimumPeakForAxis                   ( 0.3 );                               // The minimum peak height for axis to be used.
    settings->setRequestedSymmetry                    ( "" );                                // Which symmetry type (C,D,T,O or I) is requested to be detected? If none, then leave empty
    settings->setRequestedFold                        ( 0 );                                 // For C and D symmetries, which symmetry fold is requested to be detected? If none, leave 0.
    settings->setMapInversion                         ( false );                             // Should all map positions x,y,z be swapped to -x,-y,-z? Use this only if your helices have the wrong hand as a result of first runs of map computation.
    settings->setBandwidth                            ( 0 );                                 // The spherical harmonics bandwidth to which to compute. Set to 0 for automatic determination.
    settings->setSphereDistances                      ( 0.0 );                               // The distance between spheres. Use 0.0 for automatic determination.
    settings->setIntegrationOrder                     ( 0 );                                 // The order of the Gauss-Legendre integration computation. Set to 0 for automatic determination.
    settings->setTaylorSeriesCap                      ( 10 );                                // Set the Taylor series approximation cap. 10 seems like a fast and accurate value, but feel free to change.
    settings->setEnergyLevelsComputation              ( true );                              // Should energy levels descriptor be computed, assuming Distances are required (irrelevant otherwise)?
    settings->setTraceSigmaComputation                ( true );                              // Should trace sigma descriptor be computed, assuming Distances are required (irrelevant otherwise)?
    settings->setRotationFunctionComputation          ( true );                              // Should rotation function descriptor be computed, assuming Distances are required (irrelevant otherwise)?
    settings->setEnLevShellWeight                     ( 1.0 );                               // The weighting of shell distances for energy levels descriptor.
    settings->setOutputFilename                       ( "reBoxed" );                         // Filename to where re-boxed structure will be written to.
    settings->setBoundsSpace                          ( 3.0 );                               // The extra space in Angs to add to the minimal boundaries when re-boxing.
    settings->setBoundsThreshold                      ( 0 );                                 // If two boundaries are within this threshold, the smaller one will be increased to have the same value as the larger one.
    settings->setSameBoundaries                       ( false );                             // Make multiple structures have the same boundaries. This is useful for half-maps.
    settings->setPDBBFactor                           ( -1.0 );                              // Should all B-factors in a PDB file changed to this value? If no, set to negative value.
    settings->setMaskBlurFactor                       ( 350.0 );                             // If masking, what blur factor should be used? 350 seems to work for most maps.
    settings->setMaskIQR                              ( 3.0 );                               // Number of inter-quartile ranges from median to use as the masking threshold.
    settings->setMaskSaving                           ( false );                             // Should map mask be saved?
    settings->setMaskFilename                         ( "maskFile" );                        // The filename (no extension) to which the map masks will be saved into.
    settings->setAppliedMaskFilename                  ( "" );                                // The filename from which mask data will be read from.

    //================================================ Print all the settings values
//    settings->printSettings                           ( );                                   // Prints all the ProSHADE_settings values. Mostly for debugging purposes.

    //================================================ Create the structure objects
    ProSHADE_internal_data::ProSHADE_data* staticStr  = new ProSHADE_internal_data::ProSHADE_data ( ); // This line initialises the strcture object
    ProSHADE_internal_data::ProSHADE_data* movingStr  = new ProSHADE_internal_data::ProSHADE_data ( ); // This line initialises the strcture object
    
    //================================================ Read in the structures
    staticStr->readInStructure                        ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/bf/1BFO_A_dom_1.pdb", 0, settings ); // This is how a particular structure file is read into the ProSHADE object. This example uses BALBES domain 1BFO_A_dom_1.
    movingStr->readInStructure                        ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/h8/1H8N_A_dom_1.pdb", 1, settings ); // This is how a particular structure file is read into the ProSHADE object. This example uses BALBES domain 1H8N_A_dom_1.

    //================================================ Process maps
    staticStr->processInternalMap                     ( settings );  // This function does the internal map processing such as map centering, masking, invertion, phase removal, etc. for the structure which calls it.
    movingStr->processInternalMap                     ( settings );  // This function does the internal map processing such as map centering, masking, invertion, phase removal, etc. for the structure which calls it.

    //================================================ Map to spheres
    staticStr->mapToSpheres                           ( settings );  // This function maps the processed internal map onto a set of concentric spheres in preparation for spherical harmonics computation for the structure which calls it.
    movingStr->mapToSpheres                           ( settings );  // This function maps the processed internal map onto a set of concentric spheres in preparation for spherical harmonics computation for the structure which calls it.

    //================================================ Compute spherical harmonics
    staticStr->computeSphericalHarmonics              ( settings );  // This function computes the spherical harmonics for this structure.
    movingStr->computeSphericalHarmonics              ( settings );  // This function computes the spherical harmonics for this structure.

    //================================================ Compute the rotation function on the Patterson maps (they are Patterson maps as the phase was removed)
    movingStr->getOverlayRotationFunction             ( settings, staticStr ); // This function computes the rotation function between the two objects (the calling and the argument).
    
    //================================================ Get the optimal rotation Euler angles
    std::vector< proshade_double > optimalEulerRot    = movingStr->getBestRotationMapPeaksEulerAngles ( settings ); // This function returns the optimal rotation Euler angles from the rotation function.
    
    //================================================ Get optimal rotation matrix. This is not required for ProSHADE, but users may be interested in this, so it is showcased here.
    proshade_double* rotMat                           = new proshade_double[9];                 // The rotation matrix is returned as an array of 9 doubles, rows first. This is where the memory is allocated.
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ ); // Function for checking memory allocation, feel free to change to your preferred memory checking approach.
    ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( optimalEulerRot.at(0), optimalEulerRot.at(1), optimalEulerRot.at(2), rotMat ); // This is internal ProSHADE function which computes the rotation matrix from the Euler angles.
    std::cout << "Optimal rotation Euler angles are:      " << optimalEulerRot.at(0) << " ; " << optimalEulerRot.at(1) << " ; " << optimalEulerRot.at(2) << std::endl;
    std::cout << "Optimal rotation matrix is       :      " << rotMat[0] << " ; " << rotMat[1] << " ; " << rotMat[2] << std::endl;
    std::cout << "                                 :      " << rotMat[3] << " ; " << rotMat[4] << " ; " << rotMat[5] << std::endl;
    std::cout << "                                 :      " << rotMat[6] << " ; " << rotMat[7] << " ; " << rotMat[8] << std::endl;

    //================================================ Expected output
//  Optimal rotation Euler angles are:      4.01301 ; 0.698058 ; 5.40991
//  Optimal rotation matrix is       :      -0.90328 ; 0.116827 ; -0.41284
//                                   :      0.113553 ; -0.862809 ; -0.492612
//                                   :      -0.413753 ; -0.491846 ; 0.766092
    
    //================================================ Delete the Patterson maps. They are no longer needed as we will now proceed with phased maps.
    delete staticStr;
    delete movingStr;
    
    //================================================ Change the settings to reflect than now we want to work with phased maps
    settings->setPhaseUsage                           ( true );                              // Use full maps, or Patterson-like maps?
    settings->setMapResolutionChange                  ( true );                              // Should maps be re-sample to the computation resolution?
    
    //================================================ Create new structures
    staticStr                                         = new ProSHADE_internal_data::ProSHADE_data ( );
    movingStr                                         = new ProSHADE_internal_data::ProSHADE_data ( );
    
    //================================================ Read in the structures again
    staticStr->readInStructure                        ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/bf/1BFO_A_dom_1.pdb", 0, settings );  // This is how a particular structure file is read into the ProSHADE object. This example uses BALBES domain 1BFO_A_dom_1.
    movingStr->readInStructure                        ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/h8/1H8N_A_dom_1.pdb", 1, settings );  // This is how a particular structure file is read into the ProSHADE object. This example uses BALBES domain 1H8N_A_dom_1.

    //================================================ Process maps
    staticStr->processInternalMap                     ( settings );  // This function does the internal map processing such as map centering, masking, invertion, phase removal, etc. for the structure which calls it.
    movingStr->processInternalMap                     ( settings );  // This function does the internal map processing such as map centering, masking, invertion, phase removal, etc. for the structure which calls it.

    //================================================ Map to spheres - this is done only for moving structure as it will need to be rotated, but it is not needed for static structure
    movingStr->mapToSpheres                           ( settings );  // This function maps the processed internal map onto a set of concentric spheres in preparation for spherical harmonics computation for the structure which calls it.

    //================================================ Compute spherical harmonics - this is done only for moving structure as it will need to be rotated, but it is not needed for static structure
    movingStr->computeSphericalHarmonics              ( settings );  // This function computes the spherical harmonics for this structure.
    
    //================================================ Rotate the moving structure using the optimal rotation Euler angles computed before.
    movingStr->rotateMapRealSpaceInPlace              ( optimalEulerRot.at(0), optimalEulerRot.at(1), optimalEulerRot.at(2) ); // This function rotates the internal map representation in the real space using tri-linear interpolation.
    
    //================================================ Zero padding for both structures (only really applied to the smaller one, as nothing is added if the dimensions already have the requested size). This is needed to make the structures Fourier coefficients comparable
    staticStr->zeroPaddToDims                         ( int ( std::max ( staticStr->getXDim(), movingStr->getXDim() ) ),
                                                        int ( std::max ( staticStr->getYDim(), movingStr->getYDim() ) ),
                                                        int ( std::max ( staticStr->getZDim(), movingStr->getZDim() ) ) ); // This function increases the map to required sizes, adding zeroes to all new points.
    movingStr->zeroPaddToDims                         ( int ( std::max ( staticStr->getXDim(), movingStr->getXDim() ) ),
                                                        int ( std::max ( staticStr->getYDim(), movingStr->getYDim() ) ),
                                                        int ( std::max ( staticStr->getZDim(), movingStr->getZDim() ) ) ); // This function increases the map to required sizes, adding zeroes to all new points.
    
    //================================================ Compute the translation map
    movingStr->computeTranslationMap                  ( staticStr );  // This function computes the translation map for the two structures, assuming they have the same dimensions.
    
    //================================================ Find the optimal translation vector from the translation map
    std::vector< proshade_double > optimalTranslation = movingStr->getBestTranslationMapPeaksAngstrom ( staticStr ); // This function finds the best translation from the translation map using peak search algorithm and applies it to the internal map.
 
    //================================================ Find the translation vectors
    std::vector< proshade_double > rotationCentre;
    ProSHADE_internal_misc::addToDoubleVector         ( &rotationCentre, movingStr->originalPdbRotCenX );
    ProSHADE_internal_misc::addToDoubleVector         ( &rotationCentre, movingStr->originalPdbRotCenY );
    ProSHADE_internal_misc::addToDoubleVector         ( &rotationCentre, movingStr->originalPdbRotCenZ );
    
    //================================================ Write out the translations
    std::cout << "Rot. Centre to origin translation:           " << rotationCentre.at(0) << " ; " << rotationCentre.at(1) << " ; " << rotationCentre.at(2) << std::endl;
    std::cout << "Rot. Centre to optimal overlay translation:  " << optimalTranslation.at(0) << " ; " << optimalTranslation.at(1) << " ; " << optimalTranslation.at(2) << std::endl;
    
    //================================================ Expected output
//  Rot. Centre to origin translation:           9.14286 ; 17.7778 ; 17.7778
//  Rot. Centre to optimal overlay translation:  16 ; 16 ; -8
    
    //================================================ Write out the output files
    movingStr->writeOutOverlayFiles                   ( settings, optimalEulerRot.at(0), optimalEulerRot.at(1), optimalEulerRot.at(2), &rotationCentre, &optimalTranslation );
    
    //================================================ Release the settings and runProshade objects
    delete[] rotMat;
    delete staticStr;
    delete movingStr;
    delete settings;
    
    //================================================ DONE
    return                                            ( EXIT_SUCCESS );
}
