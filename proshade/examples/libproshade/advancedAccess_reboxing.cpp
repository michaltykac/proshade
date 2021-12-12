/*! \file advancedAccess_reboxing.cpp
    \brief This file showcases how the ProSHADE structure re-boxing can be computed using the dynamic library's advanced interface.

    This code is an example of how any C++ project linking the ProSHADE library can access the advanced ProSHADE interface. More specifically, in this example, a structure object is created and a structure is read into it. Then, the current object boundaries are copied
    (opitonal) and the internal structure is processed according to the settings. In this stage, map resolution and other changes are applied. Once the internal representation is processed, it is possible to do two main things. Firstly, the code shows how to creat a new ProSHADE
    internal object containing only a portion of the original object's map contained by the new boundaries. Secondaly, it shows how these new boundaries can be obtained. Finally, the new re-boxed map is written out into a file using the new ProSHADE object, which is a full
    ProSHADE object, which can be used for any other purposes, for example symmetry detection or anything else. If you want to learn more about the advanced interface and the functions available, please look into the
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

//==================================================== Main
int main ( int argc, char **argv )
{
    //================================================ Create the settings object and parse the command line arguments
    ProSHADE_Task task                                = MapManip;                            // Setting the task ahead sets most of the default settings to best values for the task.
    ProSHADE_settings* settings                       = new ProSHADE_settings ( task );      // Creating the ProSHADE_settings object, which caries all of the settings and where everything can be set.

    //================================================ Required settings
    settings->verbose                                 = -1;                                  // How verbose should the run be? -1 Means no verbal output at all.
    settings->setMapReboxing                          ( true );                              // Should the structure be re-boxed? Required masking to be done in order to be meaningful.
    
    //================================================ Further useful settings
    settings->forceP1                                 = true;                                // Should PDB files be forced to have P1 spacegroup?
    settings->setNegativeDensity                      ( true );                              // Should the negative density be removed from input files?
    settings->removeWaters                            = true;                                // Should PDB files have their water molecules removed?
    settings->firstModelOnly                          = true;                                // Should PDB files have only their first model used, or should ProSHADE use all models?
    settings->setOutputFilename                       ( "reBoxed" );                         // Filename to where re-boxed structure will be written to.
    settings->setBoundsSpace                          ( 3.0 );                               // The extra space in Angs to add to the minimal boundaries when re-boxing.
    settings->setBoundsThreshold                      ( 0 );                                 // If two boundaries are within this threshold, the smaller one will be increased to have the same value as the larger one.
    settings->setSameBoundaries                       ( false );                             // Make multiple structures have the same boundaries. This is useful for half-maps.
    settings->setMasking                              ( false );                             // Should maps be masked by blurring?
    settings->setMapCentering                         ( false );                             // Move structure COM to the centre of map box?
    settings->setMaskBlurFactor                       ( 350.0 );                             // If masking, what blur factor should be used? 350 seems to work for most maps.
    settings->setMaskIQR                              ( 3.0 );                               // Number of inter-quartile ranges from median to use as the masking threshold.
    settings->setMaskSaving                           ( false );                             // Should map mask be saved?
    settings->setMaskFilename                         ( "maskFile" );                        // The filename (no extension) to which the map masks will be saved into.
    settings->setAppliedMaskFilename                  ( "" );                                // The filename from which mask data will be read from.
    settings->setResolution                           ( 4.0 );                               // The resolution to which the calculations will be done. NOTE: Not necessarily the resolution of the structure!
    settings->setMapResolutionChange                  ( false );                             // Should maps be re-sample to the computation resolution using reciprocal-space re-sampling?
    settings->setMapResolutionChangeTriLinear         ( false );                             // Should maps be re-sample to the computation resolution using real-space tri-linear interpolation?
    
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
    settings->setProgressiveSphereMapping             ( false );                             // Should smaller spheres be less sampled? It is considerably faster, but may sacrifice some (little) accuracy.
    settings->setEnergyLevelsComputation              ( true );                              // Should energy levels descriptor be computed, assuming Distances are required (irrelevant otherwise)?
    settings->setTraceSigmaComputation                ( true );                              // Should trace sigma descriptor be computed, assuming Distances are required (irrelevant otherwise)?
    settings->setRotationFunctionComputation          ( true );                              // Should rotation function descriptor be computed, assuming Distances are required (irrelevant otherwise)?
    settings->setEnLevShellWeight                     ( 1.0 );                               // The weighting of shell distances for energy levels descriptor.
    settings->setPDBBFactor                           ( -1.0 );                              // Should all B-factors in a PDB file changed to this value? If no, set to negative value.
    settings->setPhaseUsage                           ( true );                              // Use full maps, or Patterson-like maps?
    settings->setOverlaySaveFile                      ( "overlayResuls" );                   // Filename where the overlayed moving structure should be saved.
    settings->setOverlayJsonFile                      ( "movedStructureOperations.json" );   // Filename where the overlay operations should be saved.
    settings->setNormalisation                        ( false );                             // Should internal map representation be normalised to mean 0 and standard deviation 1?
    settings->setExtraSpace                           ( 20.0 );                              // Extra space in Angs to be added when creating internap map representation. This helps avoid map effects from other cells.

    //================================================ Print all the settings values
//    settings->printSettings                           ( );                                   // Prints all the ProSHADE_settings values. Mostly for debugging purposes.

    //================================================ Create the structure objects
    ProSHADE_internal_data::ProSHADE_data* reboxStr   = new ProSHADE_internal_data::ProSHADE_data ( ); // This line initialises the structure object
    
    //================================================ Read in the structures
    reboxStr->readInStructure                         ( "./emd_5762.map", 0, settings );     // This is how a particular structure file is read into the ProSHADE object. This example uses EMDB map 5762 (PDB accession code 3J4S)

    //================================================ Get the original structure boundaries
    std::vector < proshade_signed* > oldBounds;
    ProSHADE_internal_misc::deepCopyBoundsSigPtrVector ( &oldBounds,
                                                          reboxStr->getXFromPtr(),
                                                          reboxStr->getXToPtr(),
                                                          reboxStr->getYFromPtr(),
                                                          reboxStr->getYToPtr(),
                                                          reboxStr->getZFromPtr(),
                                                          reboxStr->getZToPtr() ); // This function copies the structure boundaries into a vector of pointers, so that they can be accessible even after the structure is deleted.
    
    //================================================ Process maps
    reboxStr->processInternalMap                     ( settings );  // This function does the internal map processing such as map centering, masking, invertion, phase removal, etc. for the structure which calls it.
    
    //================================================ Create a new structure where the re-boxed data will live
    ProSHADE_internal_data::ProSHADE_data* reBoxedStr = new ProSHADE_internal_data::ProSHADE_data ( ); // This line initialises a new structure object

    //================================================ Find new boundaries (in a temporary variable)
    proshade_signed* newBounds                        = new proshade_signed[6];
    reboxStr->getReBoxBoundaries                      ( settings, newBounds ); // This function sets the new boundaries to the supplied pointer.
    
    //================================================ Create new structure from the new bounds
    reboxStr->createNewMapFromBounds                  ( settings, reBoxedStr, newBounds ); // This function will initialise a new ProSHADE object given the original object and the new boundaries.
    
    //================================================ Release the temporary variable and the no longer needed original structure object
    delete[] newBounds;
    delete reboxStr;
    
    //================================================ Save the re-boxed boundaries in a nice format
    std::vector < proshade_signed* > reboxedBounds;
    ProSHADE_internal_misc::deepCopyBoundsSigPtrVector ( &reboxedBounds,
                                                          reBoxedStr->getXFromPtr(),
                                                          reBoxedStr->getXToPtr(),
                                                          reBoxedStr->getYFromPtr(),
                                                          reBoxedStr->getYToPtr(),
                                                          reBoxedStr->getZFromPtr(),
                                                          reBoxedStr->getZToPtr() ); // This function copies the structure boundaries into a vector of pointers, so that they can be accessible even after the structure is deleted.
                                                          
    
    //================================================ Write out the rotated and translated map
    std::stringstream hlpName;
    hlpName << settings->outName << ".map";
    reBoxedStr->writeMap                              ( hlpName.str() ); // Write out the current state of the internal map representation as map.
    
    //================================================ Write out the change in bounds text (to show that we know them :-)
    std::cout << "The original structure boundaries were: " << oldBounds.at(0)[0] << " to " << oldBounds.at(0)[1] << " ; " << oldBounds.at(0)[2] << " to " << oldBounds.at(0)[3] << " and " << oldBounds.at(0)[4] << " to " << oldBounds.at(0)[5] << std::endl;
    std::cout << "The new structure boundaries are:       " << reboxedBounds.at(0)[0] << " to " << reboxedBounds.at(0)[1] << " ; " << reboxedBounds.at(0)[2] << " to " << reboxedBounds.at(0)[3] << " and " << reboxedBounds.at(0)[4] << " to " << reboxedBounds.at(0)[5] << std::endl;
    
    //================================================ Expected output
//  The original structure boundaries were: -128 to 127 ; -128 to 127 and -128 to 127
//  The new structure boundaries are:       -94 to 93 ; -94 to 93 and -136 to 135
//  // And of course the reBoxed.map file :-)
    
    //================================================ Release the settings and runProshade objects
    delete reBoxedStr;
    delete settings;
    
    //================================================ DONE
    return                                            ( EXIT_SUCCESS );
}
