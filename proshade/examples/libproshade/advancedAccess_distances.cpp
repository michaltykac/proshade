/*! \file advancedAccess_distances.cpp
    \brief This file showcases how the ProSHADE shape distances can be computed using the dynamic library's advanced interface.

    This code is an example of how any C++ project linking the ProSHADE library can access the advanced ProSHADE interface to first set the ProSHADE settings it needs, then create the ProSHADE_data structures, make these
    read in the structure files, process these and if needed and/or also manipulate the map. It then proceeds to show how these internal maps can be mapped onto a set of concentric spheres and have their spherical harmocs decomposition
    computed. Finally, once all of this is done, the code also shows how the three different shape distances can be computed. If you want to learn more about the advanced interface and the functions available, please look into the
    DOxygen documentation provided in the Documentation folder provided with the ProSHADE source codes.

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

//==================================================== Standard header for the setprecision call
#include <iomanip>

//==================================================== Main
int main ( int argc, char **argv )
{
    //================================================ Create the settings object and parse the command line arguments
    ProSHADE_Task task                                = Distances;                           // Setting the task ahead sets most of the default settings to best values for the task.
    ProSHADE_settings* settings                       = new ProSHADE_settings ( task );      // Creating the ProSHADE_settings object, which caries all of the settings and where everything can be set.

    //================================================ Set up the run
    settings->verbose                                 = -1;                                  // How verbose should the run be? -1 Means no verbal output at all.
    settings->setResolution                           ( 6.0 );                               // The resolution to which the calculations will be done. NOTE: Not necessarily the resolution of the structure!
    
    //================================================ Further useful settings
    settings->forceP1                                 = true;                                // Should PDB files be forced to have P1 spacegroup?
    settings->setNegativeDensity                      ( true );                              // Should the negative density be removed from input files?
    settings->removeWaters                            = true;                                // Should PDB files have their water molecules removed?
    settings->firstModelOnly                          = true;                                // Should PDB files have only their first model used, or should ProSHADE use all models?
    settings->setProgressiveSphereMapping             ( false );                             // Should smaller spheres be less sampled? It is considerably faster, but may sacrifice some (little) accuracy.
    settings->setMapResolutionChange                  ( false );                             // Should maps be re-sample to the computation resolution using reciprocal-space re-sampling?
    settings->setMapResolutionChangeTriLinear         ( false );                             // Should maps be re-sample to the computation resolution using real-space tri-linear interpolation?
    settings->setPDBBFactor                           ( -1.0 );                              // Should all B-factors in a PDB file changed to this value? If no, set to negative value.
    settings->setBandwidth                            ( 0 );                                 // The spherical harmonics bandwidth to which to compute. Set to 0 for automatic determination.
    settings->setPhaseUsage                           ( true );                              // Use full maps, or Patterson-like maps?
    settings->setSphereDistances                      ( 0.0 );                               // The distance between spheres. Use 0.0 for automatic determination.
    settings->setIntegrationOrder                     ( 0 );                                 // The order of the Gauss-Legendre integration computation. Set to 0 for automatic determination.
    settings->setTaylorSeriesCap                      ( 10 );                                // Set the Taylor series approximation cap. 10 seems like a fast and accurate value, but feel free to change.
    settings->setNormalisation                        ( false );                             // Should internal map representation be normalised to mean 0 and standard deviation 1?
    settings->setMapInversion                         ( false );                             // Should all map positions x,y,z be swapped to -x,-y,-z? Use this only if your helices have the wrong hand as a result of first runs of map computation.
    settings->setMasking                              ( false );                             // Should maps be masked by blurring?
    settings->setMapReboxing                          ( false );                             // Should the structure be re-boxed? Required masking to be done in order to be meaningful.
    settings->setMapCentering                         ( false );                             // Move structure COM to the centre of map box?
    settings->setEnergyLevelsComputation              ( true );                              // Should energy levels descriptor be computed, assuming Distances are required (irrelevant otherwise)?
    settings->setTraceSigmaComputation                ( true );                              // Should trace sigma descriptor be computed, assuming Distances are required (irrelevant otherwise)?
    settings->setRotationFunctionComputation          ( true );                              // Should rotation function descriptor be computed, assuming Distances are required (irrelevant otherwise)?
    settings->setEnLevShellWeight                     ( 1.0 );                               // The weighting of shell distances for energy levels descriptor.


    //================================================ All other (possibly other tasks related) settings
    settings->setSymmetryCentreSearch                 ( false );                             // Should symmetry centre be searched for? Takes a lot of time...
    settings->setBicubicInterpolationSearch           ( true );                              // Should bi-cubic interpolation between peak grid indices be done?
    settings->setMaxSymmetryFold                      ( 30 );                                // The maximum prime number fold that will be searched for.
    settings->setFSCThreshold                         ( 0.75 );                              // Sets the minimum FSC threshold for axis to be considered detected.
    settings->setPeakThreshold                        ( 0.80 );                              // Sets the minimum peak height threshold for axis to be considered possible.
    settings->setMaskBlurFactor                       ( 350.0 );                             // If masking, what blur factor should be used? 350 seems to work for most maps.
    settings->setMaskIQR                              ( 3.0 );                               // Number of inter-quartile ranges from median to use as the masking threshold.
    settings->setMaskSaving                           ( false );                             // Should map mask be saved?
    settings->setMaskFilename                         ( "maskFile" );                        // The filename (no extension) to which the map masks will be saved into.
    settings->setAppliedMaskFilename                  ( "" );                                // The filename from which mask data will be read from.
    settings->setBoundsSpace                          ( 3.0 );                               // The extra space in Angs to add to the minimal boundaries when re-boxing.
    settings->setBoundsThreshold                      ( 0 );                                 // If two boundaries are within this threshold, the smaller one will be increased to have the same value as the larger one.
    settings->setSameBoundaries                       ( false );                             // Make multiple structures have the same boundaries. This is useful for half-maps.
    settings->setOutputFilename                       ( "reBoxed" );                         // Filename to where re-boxed structure will be written to.
    settings->setExtraSpace                           ( 10.0 );                              // Extra space in Angs to be added when creating internap map representation. This helps avoid map effects from other cells.
    settings->setPeakNeighboursNumber                 ( 1 );                                 // Numer of points in each direction which needs to be lower in order for the central point to be considered a peak.
    settings->setPeakNaiveNoIQR                       ( -999.9 );                            // Peak searching threshold for too low peaks in number of inter-quartile ranges from median of the non-peak point values.
    settings->setMissingPeakThreshold                 ( 0.3 );                               // Fraction of peaks that can be missing for missing axis search to be initiated.
    settings->setAxisComparisonThreshold              ( 0.1 );                               // The dot product difference within which two axes are considered the same.
    settings->setMinimumPeakForAxis                   ( 0.3 );                               // The minimum peak height for axis to be used.
    settings->setRequestedSymmetry                    ( "" );                                // Which symmetry type (C,D,T,O or I) is requested to be detected? If none, then leave empty
    settings->setRequestedFold                        ( 0 );                                 // For C and D symmetries, which symmetry fold is requested to be detected? If none, leave 0.
    settings->setOverlaySaveFile                      ( "moved" );                           // Filename where the overlayed moving structure should be saved.
    settings->setOverlayJsonFile                      ( "movedStructureOperations.json" );   // Filename where the overlay operations should be saved.

    //================================================ Print all the settings values
//    settings->printSettings                           ( );                                   // Prints all the ProSHADE_settings values. Mostly for debugging purposes.
    
    //================================================ Create the structure objects
    ProSHADE_internal_data::ProSHADE_data* structure1 = new ProSHADE_internal_data::ProSHADE_data ( ); // This line initialises the strcture object
    ProSHADE_internal_data::ProSHADE_data* structure2 = new ProSHADE_internal_data::ProSHADE_data ( ); // This line initialises the strcture object
    
    //================================================ Read in the structures
    structure1->readInStructure                       ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/bf/1BFO_A_dom_1.pdb", 0, settings ); // This is how a particular structure file is read into the ProSHADE object.
    structure2->readInStructure                       ( "/Users/mysak/LMB/1_ProteinDomains/0_DOMS/h8/1H8N_A_dom_1.pdb", 1, settings ); // This is how a particular structure file is read into the ProSHADE object.

    //================================================ Process maps
    structure1->processInternalMap                    ( settings );  // This function does the internal map processing such as map centering, masking, invertion, phase removal, etc. for the structure which calls it.
    structure2->processInternalMap                    ( settings );  // This function does the internal map processing such as map centering, masking, invertion, phase removal, etc. for the structure which calls it.

    //================================================ Map to spheres
    structure1->mapToSpheres                          ( settings );  // This function maps the processed internal map onto a set of concentric spheres in preparation for spherical harmonics computation for the structure which calls it.
    structure2->mapToSpheres                          ( settings );  // This function maps the processed internal map onto a set of concentric spheres in preparation for spherical harmonics computation for the structure which calls it.

    //================================================ Compute spherical harmonics
    structure1->computeSphericalHarmonics             ( settings );  // This function computes the spherical harmonics for this structure.
    structure2->computeSphericalHarmonics             ( settings );  // This function computes the spherical harmonics for this structure.

    
    //================================================ Get results - these are already computed by the previous line, they are just copied by these commands.
    proshade_double energyDistance                    = ProSHADE_internal_distances::computeEnergyLevelsDescriptor     ( structure1, structure2, settings ); // This is how the energy levels distance is computed between structure1 and structure2.
    proshade_double traceDistance                     = ProSHADE_internal_distances::computeTraceSigmaDescriptor       ( structure1, structure2, settings ); // This is how the trace sigma distance is computed between structure1 and structure2.
    proshade_double rotFunDistance                    = ProSHADE_internal_distances::computeRotationFunctionDescriptor ( structure1, structure2, settings ); // This is how the rotation function distance is computed between structure1 and structure2.

    //================================================ Print results
    std::cout << std::setprecision(5) << std::flush;
    std::cout << "Energy levels distance       : " << energyDistance << std::endl;
    std::cout << "Trace sigma distance         : " << traceDistance  << std::endl;
    std::cout << "Rotation function distance   : " << rotFunDistance << std::endl;
    
    //================================================ Expected output
//  Energy levels distance       : 0.85561
//  Trace sigma distance         : 0.96467
//  Rotation function distance   : 0.62453

    //================================================ Release the settings and runProshade objects
    delete structure1;
    delete structure2;
    delete settings;
    
    //================================================ DONE
    return                                            ( EXIT_SUCCESS );
}
