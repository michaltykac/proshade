/*! \file pyProSHADE.cpp
    \brief This file contains the PyBind11 bindings for the ProSHADE_settings class.
    
    This file provides the bindings for hte ProSHADE_settings class members and functions. It also defines several python specific functions (written as C++ lambda functions) which allow direct access
    to the computed results as numpy.ndarrays, while making sure the memory is released correctly for such cross-language shared variables.
    
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

//==================================================== Include PyBind11 header
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

//==================================================== Add the ProSHADE_settings and ProSHADE_run classes to the PyBind11 module
void add_settingsClass ( pybind11::module& pyProSHADE )
{
    //================================================ Export the ProSHADE_settings class
    pybind11::class_ < ProSHADE_settings >            ( pyProSHADE, "ProSHADE_settings" )
        
        //============================================ Constructors (destructors do not need wrappers???)
        .def                                          ( pybind11::init < > ( ) )
        .def                                          ( pybind11::init < ProSHADE_Task > ( ), pybind11::arg ( "task" ) )
    
        //============================================ Member variables
        .def_readwrite                                ( "task",                                 &ProSHADE_settings::task                                )
                                
        .def_readwrite                                ( "inputFiles",                           &ProSHADE_settings::inputFiles                          )
        .def_readwrite                                ( "forceP1",                              &ProSHADE_settings::forceP1                             )
        .def_readwrite                                ( "removeNegativeDensity",                &ProSHADE_settings::removeNegativeDensity               )
        .def_readwrite                                ( "removeWaters",                         &ProSHADE_settings::removeWaters                        )
        .def_readwrite                                ( "firstModelOnly",                       &ProSHADE_settings::firstModelOnly                      )
                
        .def_readwrite                                ( "requestedResolution",                  &ProSHADE_settings::requestedResolution                 )
        .def_readwrite                                ( "changeMapResolution",                  &ProSHADE_settings::changeMapResolution                 )
        .def_readwrite                                ( "changeMapResolutionTriLinear",         &ProSHADE_settings::changeMapResolutionTriLinear        )
    
        .def_readwrite                                ( "pdbBFactorNewVal",                     &ProSHADE_settings::pdbBFactorNewVal                    )
    
        .def_readwrite                                ( "maxBandwidth",                         &ProSHADE_settings::maxBandwidth                        )
        .def_readwrite                                ( "rotationUncertainty",                  &ProSHADE_settings::rotationUncertainty                 )
    
        .def_readwrite                                ( "usePhase",                             &ProSHADE_settings::usePhase                            )
    
        .def_readwrite                                ( "maxSphereDists",                       &ProSHADE_settings::maxSphereDists                      )
    
        .def_readwrite                                ( "integOrder",                           &ProSHADE_settings::integOrder                          )
        .def_readwrite                                ( "taylorSeriesCap",                      &ProSHADE_settings::taylorSeriesCap                     )
    
        .def_readwrite                                ( "normaliseMap",                         &ProSHADE_settings::normaliseMap                        )
    
        .def_readwrite                                ( "invertMap",                            &ProSHADE_settings::invertMap                           )
    
        .def_readwrite                                ( "blurFactor",                           &ProSHADE_settings::blurFactor                          )
        .def_readwrite                                ( "maskingThresholdIQRs",                 &ProSHADE_settings::maskingThresholdIQRs                )
        .def_readwrite                                ( "maskMap",                              &ProSHADE_settings::maskMap                             )
        .def_readwrite                                ( "useCorrelationMasking",                &ProSHADE_settings::useCorrelationMasking               )
        .def_readwrite                                ( "halfMapKernel",                        &ProSHADE_settings::halfMapKernel                       )
        .def_readwrite                                ( "correlationKernel",                    &ProSHADE_settings::correlationKernel                   )
        .def_readwrite                                ( "saveMask",                             &ProSHADE_settings::saveMask                            )
        .def_readwrite                                ( "maskFileName",                         &ProSHADE_settings::maskFileName                        )
        .def_readwrite                                ( "appliedMaskFileName",                  &ProSHADE_settings::appliedMaskFileName                 )
    
        .def_readwrite                                ( "fourierWeightsFileName",               &ProSHADE_settings::fourierWeightsFileName              )
    
        .def_readwrite                                ( "reBoxMap",                             &ProSHADE_settings::reBoxMap                            )
        .def_readwrite                                ( "boundsExtraSpace",                     &ProSHADE_settings::boundsExtraSpace                    )
        .def_readwrite                                ( "boundsSimilarityThreshold",            &ProSHADE_settings::boundsSimilarityThreshold           )
        .def_readwrite                                ( "useSameBounds",                        &ProSHADE_settings::useSameBounds                       )
        .def_readwrite                                ( "forceBounds",                          &ProSHADE_settings::forceBounds                         )
    
        .def_readwrite                                ( "moveToCOM",                            &ProSHADE_settings::moveToCOM                           )
    
        .def_readwrite                                ( "addExtraSpace",                        &ProSHADE_settings::addExtraSpace                       )
        .def_readwrite                                ( "coOrdsExtraSpace",                     &ProSHADE_settings::coOrdsExtraSpace                    )
    
        .def_readwrite                                ( "progressiveSphereMapping",             &ProSHADE_settings::progressiveSphereMapping            )
    
        .def_readwrite                                ( "outName",                              &ProSHADE_settings::outName                             )
    
        .def_readwrite                                ( "computeEnergyLevelsDesc",              &ProSHADE_settings::computeEnergyLevelsDesc             )
        .def_readwrite                                ( "enLevMatrixPowerWeight",               &ProSHADE_settings::enLevMatrixPowerWeight              )
        .def_readwrite                                ( "computeTraceSigmaDesc",                &ProSHADE_settings::computeTraceSigmaDesc               )
        .def_readwrite                                ( "computeRotationFuncDesc",              &ProSHADE_settings::computeRotationFuncDesc             )
    
        .def_readwrite                                ( "peakNeighbours",                       &ProSHADE_settings::peakNeighbours                      )
        .def_readwrite                                ( "noIQRsFromMedianNaivePeak",            &ProSHADE_settings::noIQRsFromMedianNaivePeak           )
    
        .def_readwrite                                ( "smoothingFactor",                      &ProSHADE_settings::smoothingFactor                     )
    
        .def_readwrite                                ( "symMissPeakThres",                     &ProSHADE_settings::symMissPeakThres                    )
        .def_readwrite                                ( "axisErrTolerance",                     &ProSHADE_settings::axisErrTolerance                    )
        .def_readwrite                                ( "axisErrToleranceDefault",              &ProSHADE_settings::axisErrToleranceDefault             )
        .def_readwrite                                ( "minSymPeak",                           &ProSHADE_settings::minSymPeak                          )
        .def_readwrite                                ( "recommendedSymmetryType",              &ProSHADE_settings::recommendedSymmetryType             )
        .def_readwrite                                ( "recommendedSymmetryFold",              &ProSHADE_settings::recommendedSymmetryFold             )
        .def_readwrite                                ( "requestedSymmetryType",                &ProSHADE_settings::requestedSymmetryType               )
        .def_readwrite                                ( "requestedSymmetryFold",                &ProSHADE_settings::requestedSymmetryFold               )
        .def_readwrite                                ( "useBiCubicInterpolationOnPeaks",       &ProSHADE_settings::useBiCubicInterpolationOnPeaks      )
        .def_readwrite                                ( "maxSymmetryFold",                      &ProSHADE_settings::maxSymmetryFold                     )
        .def_readwrite                                ( "fscThreshold",                         &ProSHADE_settings::fscThreshold                        )
        .def_readwrite                                ( "peakThresholdMin",                     &ProSHADE_settings::peakThresholdMin                    )
    
        .def_readwrite                                ( "overlayStructureName",                 &ProSHADE_settings::overlayStructureName                )
        .def_readwrite                                ( "rotTrsJSONFile",                       &ProSHADE_settings::rotTrsJSONFile                      )
    
        .def_readwrite                                ( "verbose",                              &ProSHADE_settings::verbose                             )
        .def_readwrite                                ( "messageShift",                         &ProSHADE_settings::messageShift                        )
    
        //============================================ Mutators
        .def                                          ( "addStructure",                         &ProSHADE_settings::addStructure,                           "Adds a structure file name to the appropriate variable.",                                                                  pybind11::arg ( "structure"     ) )
        .def                                          ( "setResolution",                        &ProSHADE_settings::setResolution,                          "This function sets the resolution in the appropriate variable.",                                                           pybind11::arg ( "resolution"    ) )
        .def                                          ( "setPDBBFactor",                        &ProSHADE_settings::setPDBBFactor,                          "Sets the requested B-factor value for PDB files in the appropriate variable.",                                             pybind11::arg ( "newBF"         ) )
        .def                                          ( "setNormalisation",                     &ProSHADE_settings::setNormalisation,                       "Sets the requested map normalisation value in the appropriate variable.",                                                  pybind11::arg ( "normalise"     ) )
        .def                                          ( "setMapInversion",                      &ProSHADE_settings::setMapInversion,                        "Sets the requested map inversion value in the appropriate variable.",                                                      pybind11::arg ( "mInv"          ) )
        .def                                          ( "setVerbosity",                         &ProSHADE_settings::setVerbosity,                           "Sets the requested verbosity in the appropriate variable.",                                                                pybind11::arg ( "verbosity"     ) )
        .def                                          ( "setMaskBlurFactor",                    &ProSHADE_settings::setMaskBlurFactor,                      "Sets the requested map blurring factor in the appropriate variable.",                                                      pybind11::arg ( "blurFac"       ) )
        .def                                          ( "setMaskIQR",                           &ProSHADE_settings::setMaskIQR,                             "Sets the requested number of IQRs for masking threshold in the appropriate variable.",                                     pybind11::arg ( "noIQRs"        ) )
        .def                                          ( "setMasking",                           &ProSHADE_settings::setMasking,                             "Sets the requested map masking decision value in the appropriate variable.",                                               pybind11::arg ( "mask"          ) )
        .def                                          ( "setCorrelationMasking",                &ProSHADE_settings::setCorrelationMasking,                  "Sets the requested map masking type in the appropriate variable.",                                                         pybind11::arg ( "corMask"       ) )
        .def                                          ( "setTypicalNoiseSize",                  &ProSHADE_settings::setTypicalNoiseSize,                    "Sets the requested \"fake\" half-map kernel in the appropriate variable.",                                                 pybind11::arg ( "typNoi"        ) )
        .def                                          ( "setMinimumMaskSize",                   &ProSHADE_settings::setMinimumMaskSize,                     "Sets the requested minimum mask size.",                                                                                    pybind11::arg ( "minMS"         ) )
        .def                                          ( "setMaskSaving",                        &ProSHADE_settings::setMaskSaving,                          "Sets whether the mask should be saved.",                                                                                   pybind11::arg ( "savMsk"        ) )
        .def                                          ( "setMaskFilename",                      &ProSHADE_settings::setMaskFilename,                        "Sets where the mask should be saved.",                                                                                     pybind11::arg ( "mskFln"        ) )
        .def                                          ( "setAppliedMaskFilename",               &ProSHADE_settings::setAppliedMaskFilename,                 "Sets the filename of the mask data that should be applied to the input map.",                                             pybind11::arg ( "mskFln"        ) )
        .def                                          ( "setFourierWeightsFilename",            &ProSHADE_settings::setFourierWeightsFilename,              "Sets the filename of the Fourier weights data that should be applied to the input map.",                                  pybind11::arg ( "fWgFln"        ) )
        .def                                          ( "setMapReboxing",                       &ProSHADE_settings::setMapReboxing,                         "Sets whether re-boxing needs to be done in the appropriate variable.",                                                     pybind11::arg ( "reBx"          ) )
        .def                                          ( "setBoundsSpace",                       &ProSHADE_settings::setBoundsSpace,                         "Sets the requested number of angstroms for extra space in re-boxing in the appropriate variable.",                         pybind11::arg ( "boundsExSp"    ) )
        .def                                          ( "setBoundsThreshold",                   &ProSHADE_settings::setBoundsThreshold,                     "Sets the threshold for number of indices difference acceptable to make index sizes same in the appropriate variable.",     pybind11::arg ( "boundsThres"   ) )
        .def                                          ( "setSameBoundaries",                    &ProSHADE_settings::setSameBoundaries,                      "Sets whether same boundaries should be used in the appropriate variable.",                                                 pybind11::arg ( "sameB"         ) )
        .def                                          ( "setOutputFilename",                    &ProSHADE_settings::setOutputFilename,                      "Sets the requested output file name in the appropriate variable.",                                                         pybind11::arg ( "oFileName"     ) )
        .def                                          ( "setMapResolutionChange",               &ProSHADE_settings::setMapResolutionChange,                 "Sets the requested map resolution change decision in the appropriate variable.",                                           pybind11::arg ( "mrChange"      ) )
        .def                                          ( "setMapResolutionChangeTriLinear",      &ProSHADE_settings::setMapResolutionChangeTriLinear,        "Sets the requested map resolution change decision using tri-linear interpolation in the appropriate variable.",            pybind11::arg ( "mrChange"      ) )
        .def                                          ( "setMapCentering",                      &ProSHADE_settings::setMapCentering,                        "Sets the requested map centering decision value in the appropriate variable.",                                             pybind11::arg ( "com"           ) )
        .def                                          ( "setExtraSpace",                        &ProSHADE_settings::setExtraSpace,                          "Sets the requested map extra space value in the appropriate variable.",                                                    pybind11::arg ( "exSpace"       ) )
        .def                                          ( "setCoordExtraSpace",                   &ProSHADE_settings::setCoordExtraSpace,                     "Sets the requested co-ordinates extra space value in the appropriate variable.",                                                    pybind11::arg ( "exSpace"       ) )
        .def                                          ( "setBoxCentering",                      &ProSHADE_settings::setBoxCentering,                        "Sets the requested centre of box to be at the real space co-ordinates supplied.",                                           pybind11::arg ( "xPos"          ), pybind11::arg ( "yPos"          ), pybind11::arg ( "zPos"          ) )
        .def                                          ( "setBandwidth",                         &ProSHADE_settings::setBandwidth,                           "Sets the requested spherical harmonics bandwidth in the appropriate variable.",                                            pybind11::arg ( "band"          ) )
        .def                                          ( "setProgressiveSphereMapping",          &ProSHADE_settings::setProgressiveSphereMapping,            "Sets the requested sphere mapping value settings approach in the appropriate variable.",                                   pybind11::arg ( "progSphMap"    ) )
        .def                                          ( "setSphereDistances",                   &ProSHADE_settings::setSphereDistances,                     "Sets the requested distance between spheres in the appropriate variable.",                                                 pybind11::arg ( "sphDist"       ) )
        .def                                          ( "setIntegrationOrder",                  &ProSHADE_settings::setIntegrationOrder,                    "Sets the requested order for the Gauss-Legendre integration in the appropriate variable.",                                 pybind11::arg ( "intOrd"        ) )
        .def                                          ( "setTaylorSeriesCap",                   &ProSHADE_settings::setTaylorSeriesCap,                     "Sets the requested Taylor series cap for the Gauss-Legendre integration in the appropriate variable.",                     pybind11::arg ( "tayCap"        ) )
        .def                                          ( "setEnergyLevelsComputation",           &ProSHADE_settings::setEnergyLevelsComputation,             "Sets whether the energy level distance descriptor should be computed.",                                                    pybind11::arg ( "enLevDesc"     ) )
        .def                                          ( "setTraceSigmaComputation",             &ProSHADE_settings::setTraceSigmaComputation,               "Sets whether the trace sigma distance descriptor should be computed.",                                                     pybind11::arg ( "trSigVal"      ) )
        .def                                          ( "setRotationFunctionComputation",       &ProSHADE_settings::setRotationFunctionComputation,         "Sets whether the rotation function distance descriptor should be computed.",                                               pybind11::arg ( "rotfVal"       ) )
        .def                                          ( "setPeakNeighboursNumber",              &ProSHADE_settings::setPeakNeighboursNumber,                "Sets the number of neighbour values that have to be smaller for an index to be considered a peak.",                        pybind11::arg ( "pkS"           ) )
        .def                                          ( "setPeakNaiveNoIQR",                    &ProSHADE_settings::setPeakNaiveNoIQR,                      "Sets the number of IQRs from the median for threshold height a peak needs to be considered a peak.",                       pybind11::arg ( "noIQRs"        ) )
        .def                                          ( "setPhaseUsage",                        &ProSHADE_settings::setPhaseUsage,                          "Sets whether the phase information will be used.",                                                                         pybind11::arg ( "phaseUsage"    ) )
        .def                                          ( "setEnLevShellWeight",                  &ProSHADE_settings::setEnLevShellWeight,                    "Sets the weight of shell position for the energy levels computation.",                                                  pybind11::arg ( "mPower"        ) )
        .def                                          ( "setGroupingSmoothingFactor",           &ProSHADE_settings::setGroupingSmoothingFactor,             "Sets the grouping smoothing factor into the proper variable.",                                                     pybind11::arg ( "smFact"        ) )
        .def                                          ( "setMissingPeakThreshold",              &ProSHADE_settings::setMissingPeakThreshold,                "Sets the threshold for starting the missing peaks procedure.",                                                    pybind11::arg ( "mpThres"       ) )
        .def                                          ( "setSymmetryCentreSearch",              &ProSHADE_settings::setSymmetryCentreSearch,                "Sets the symmetry centre search on or off.",                                                                               pybind11::arg ( "sCen"          ) )
        .def                                          ( "setAxisComparisonThreshold",           &ProSHADE_settings::setAxisComparisonThreshold,             "Sets the threshold for matching symmetry axes.",                                                                           pybind11::arg ( "axThres"       ) )
        .def                                          ( "setAxisComparisonThresholdBehaviour",  &ProSHADE_settings::setAxisComparisonThresholdBehaviour,    "Sets the automatic symmetry axis tolerance decreasing.",                                                                   pybind11::arg ( "behav"         ) )
        .def                                          ( "setMinimumPeakForAxis",                &ProSHADE_settings::setMinimumPeakForAxis,                  "Sets the minimum peak height for symmetry axis to be considered.",                                                         pybind11::arg ( "minSP"         ) )
        .def                                          ( "setRecommendedSymmetry",               &ProSHADE_settings::setRecommendedSymmetry,                 "Sets the ProSHADE detected symmetry type.",                                                                                pybind11::arg ( "val"           ) )
        .def                                          ( "setRecommendedFold",                   &ProSHADE_settings::setRecommendedFold,                     "Sets the ProSHADE detected symmetry fold.",                                                                                pybind11::arg ( "val"           ) )
        .def                                          ( "setRequestedSymmetry",                 &ProSHADE_settings::setRequestedSymmetry,                   "Sets the user requested symmetry type.",                                                                                   pybind11::arg ( "val"           ) )
        .def                                          ( "setRequestedFold",                     &ProSHADE_settings::setRequestedFold,                       "Sets the user requested symmetry fold.",                                                                                   pybind11::arg ( "val"           ) )
        .def                                          ( "setDetectedSymmetry",                  &ProSHADE_settings::setDetectedSymmetry,                    "Sets the final detected symmetry axes information.",                                                                       pybind11::arg ( "sym"           ) )
        .def                                          ( "setOverlaySaveFile",                   &ProSHADE_settings::setOverlaySaveFile,                     "Sets the filename to which the overlay structure is to be save into.",                                                    pybind11::arg ( "filename"      ) )
        .def                                          ( "setOverlayJsonFile",                   &ProSHADE_settings::setOverlayJsonFile,                     "Sets the filename to which the overlay operations are to be save into.",                                                 pybind11::arg ( "filename"      ) )
        .def                                          ( "setBicubicInterpolationSearch",        &ProSHADE_settings::setBicubicInterpolationSearch,          "Sets the bicubic interpolation on peaks.",                                                                                 pybind11::arg ( "bicubPeaks"    ) )
        .def                                          ( "setMaxSymmetryFold",                   &ProSHADE_settings::setMaxSymmetryFold,                     "Sets the maximum symmetry fold (well, the maximum prime symmetry fold).",                                               pybind11::arg ( "maxFold"       ) )
        .def                                          ( "setFSCThreshold",                      &ProSHADE_settings::setFSCThreshold,                        "Sets the minimum FSC threshold for axis to be considered detected.",                                                     pybind11::arg ( "fscThr"        ) )
        .def                                          ( "setPeakThreshold",                     &ProSHADE_settings::setPeakThreshold,                       "Sets the minimum peak height threshold for axis to be considered possible.",                                          pybind11::arg ( "peakThr"       ) )
        .def                                          ( "setNegativeDensity",                   &ProSHADE_settings::setNegativeDensity,                     "Sets the internal variable deciding whether input files negative density should be removed.",                           pybind11::arg ( "nDens"         ) )
    
        //============================================ Command line parsing
        .def                                          ( "getCommandLineParams",
                                                        [] ( ProSHADE_settings &self, std::vector < std::string > args )
                                                        {
                                                            std::vector < char * > cstrs; cstrs.reserve ( args.size() );
            
                                                            for ( auto &s : args )
                                                                cstrs.push_back ( const_cast < char * > ( s.c_str ( ) ) );
                                                            
                                                            return self.getCommandLineParams ( static_cast< int > ( cstrs.size ( ) ), cstrs.data ( ) );
                                                        }, "This function takes a VectorOfStrings and parses it as if it were command line arguments, filling in the calling ProSHADE_settings class with the values." )
        
        //============================================ Debugging
        .def                                          ( "printSettings", &ProSHADE_settings::printSettings, "This function prints the current values in the settings object." )
    
        //============================================ Description
        .def                                          ( "__repr__", [] ( ) { return "<ProSHADE_settings class object> (Settings class is used to set all settings values in a single place)"; } );
    
    //================================================ Export the ProSHADE_run class
    pybind11::class_ < ProSHADE_run >                 ( pyProSHADE, "ProSHADE_run" )
    
        //============================================ Constructors (destructors do not need wrappers???)
        .def                                          ( pybind11::init < ProSHADE_settings* > ( ) )
    
        //============================================ General accessors
        .def                                          ( "getNoStructures", &ProSHADE_run::getNoStructures, "This function returns the number of structures used." )
        .def                                          ( "getVerbose", &ProSHADE_run::getVerbose, "This function returns the verbose value." )
    
        //============================================ Distances results accessor functions wrapped as lambda functions for numpy return types
        .def                                          ( "getEnergyLevelsVector",
                                                        [] ( ProSHADE_run &self ) -> pybind11::array_t < float >
                                                        {
                                                            //== Get the values
                                                            std::vector< proshade_double > vals = self.getEnergyLevelsVector ();
            
                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[static_cast<unsigned int> (vals.size())];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );
            
                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( vals.size() ); iter++ ) { npVals[iter] = static_cast< float > ( vals.at(iter) ); }
            
                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleEnLevs ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );
            
                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { static_cast<unsigned int> (vals.size()) },  // Shape
                                                                                                                            { sizeof(float) },                            // C-stype strides
                                                                                                                            npVals,                                       // Data
                                                                                                                            pyCapsuleEnLevs );                            // Capsule (C++ destructor, basically)
            
                                                            //== Done
                                                            return ( retArr );
                                                        }, "This function returns the energy level distances vector from the first to all other structures." )
    
        .def                                          ( "getTraceSigmaVector",
                                                        [] ( ProSHADE_run &self ) -> pybind11::array_t < float >
                                                        {
                                                            //== Get the values
                                                            std::vector< proshade_double > vals = self.getTraceSigmaVector ();
        
                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[static_cast<unsigned int> (vals.size())];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );
        
                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( vals.size() ); iter++ ) { npVals[iter] = static_cast< float > ( vals.at(iter) ); }
        
                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleTrSigs ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );
        
                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { static_cast<unsigned int> (vals.size()) },  // Shape
                                                                                                                            { sizeof(float) },                            // C-stype strides
                                                                                                                            npVals,                                       // Data
                                                                                                                            pyCapsuleTrSigs );                            // Capsule (C++ destructor, basically)
        
                                                            //== Done
                                                            return ( retArr );
                                                        }, "This function returns the trace sigma distances vector from the first to all other structures." )
    
        .def                                          ( "getRotationFunctionVector",
                                                        [] ( ProSHADE_run &self ) -> pybind11::array_t < float >
                                                        {
                                                            //== Get the values
                                                            std::vector< proshade_double > vals = self.getRotationFunctionVector ();
        
                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[static_cast<unsigned int> (vals.size())];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );
        
                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( vals.size() ); iter++ ) { npVals[iter] = static_cast< float > ( vals.at(iter) ); }
        
                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleRotFun ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );
        
                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { static_cast<unsigned int> (vals.size()) },  // Shape
                                                                                                                            { sizeof(float) },                            // C-stype strides
                                                                                                                            npVals,                                       // Data
                                                                                                                            pyCapsuleRotFun );                            // Capsule (C++ destructor, basically)
        
                                                            //== Done
                                                            return ( retArr );
                                                        }, "This function returns the full rotation function distances vector from the first to all other structures." )
    
        //============================================ Symmetry results accessor functions
        .def                                          ( "getSymmetryType", &ProSHADE_run::getSymmetryType, "This is the main accessor function for the user to get to know what symmetry type ProSHADE has detected and recommends." )
        .def                                          ( "getSymmetryFold", &ProSHADE_run::getSymmetryFold, "This is the main accessor function for the user to get to know what symmetry fold ProSHADE has detected and recommends." )
        .def                                          ( "getSymmetryAxis", &ProSHADE_run::getSymmetryAxis, "This function returns a single symmetry axis as a vector of strings from the recommended symmetry axes list.", pybind11::arg ( "axisNo" ) )
        .def                                          ( "getAllCSyms",
                                                        [] ( ProSHADE_run &self ) -> pybind11::array_t < float >
                                                        {
                                                            //== Get the values
                                                            std::vector< std::vector< proshade_double > > vals = self.getAllCSyms ();
    
                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[static_cast<unsigned int> ( vals.size() * 7 )];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );
            
                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( vals.size() ); iter++ ) { for ( proshade_unsign it = 0; it < 7; it++ ) { npVals[(iter*7)+it] = static_cast< float > ( vals.at(iter).at(it) ); } }
            
                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleSymList ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );
    
                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { static_cast<int> ( vals.size() ), static_cast<int> ( 7 ) },  // Shape
                                                                                                                            { 7 * sizeof(float), sizeof(float) },                          // C-stype strides
                                                                                                                            npVals,                                                        // Data
                                                                                                                            pyCapsuleSymList );                                            // Capsule
    
                                                            //== Done
                                                            return ( retArr );
                                                        }, "This function returns a all symmetry axes as a 2D numpy array." )
        .def                                          ( "getMapCOMProcessChange",
                                                       [] ( ProSHADE_run &self ) -> pybind11::array_t < float >
                                                       {
                                                            //== Get the values
                                                            std::vector< proshade_double > vals = self.getMapCOMProcessChange ();

                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[static_cast<unsigned int> ( 3 )];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );
        
                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < 3; iter++ ) { npVals[iter] = static_cast< float > ( vals.at(iter) ); }
        
                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleSymShift ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );

                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { static_cast<int> ( vals.size() ) },      // Shape
                                                                                                                          { sizeof(float) },                           // C-stype strides
                                                                                                                          npVals,                                      // Data
                                                                                                                          pyCapsuleSymShift );                         // Capsule

                                                            //== Done
                                                            return ( retArr );
                                                        }, "This function returns the shift in Angstrom applied to the internal map representation in order to align its COM with the centre of box." )
    
        //============================================ Reboxing results accessor functions as lambda functions directly returning numpy arrays
        .def                                          ( "getOriginalBounds",
                                                        [] ( ProSHADE_run &self, proshade_unsign strNo ) -> pybind11::array_t < float >
                                                        {
                                                            //== Get values
                                                            std::vector< proshade_signed > vals = self.getOriginalBounds ( strNo );
            
                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[static_cast<proshade_unsign> ( vals.size() )];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );
        
                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( vals.size() ); iter++ ) { npVals[iter] = vals.at(iter); }
        
                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleOrigBnds ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );

                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { static_cast<proshade_unsign> ( vals.size() ) },  // Shape
                                                                                                                            { sizeof(float) },                                 // C-stype strides
                                                                                                                            npVals,                                            // Data
                                                                                                                            pyCapsuleOrigBnds );                               // Capsule (C++ destructor, basically)

                                                            //== Done
                                                            return ( retArr );
                                                        }, "This function returns the original structure boundaries as numpy array." )
    
        .def                                          ( "getReBoxedBounds",
                                                        [] ( ProSHADE_run &self, proshade_unsign strNo ) -> pybind11::array_t < float >
                                                        {
                                                            //== Get values
                                                            std::vector< proshade_signed > vals = self.getReBoxedBounds ( strNo );
            
                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[static_cast<proshade_unsign> ( vals.size() )];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );

                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( vals.size() ); iter++ ) { npVals[iter] = vals.at(iter); }

                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleReBoBnds ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );

                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { static_cast<proshade_unsign> ( vals.size() ) },  // Shape
                                                                                                                            { sizeof(float) },                                 // C-stype strides
                                                                                                                            npVals,                                            // Data
                                                                                                                            pyCapsuleReBoBnds );                               // Capsule (C++ destructor, basically)

                                                            //== Done
                                                            return ( retArr );
                                                        }, "This function returns the re-boxed structure boundaries as numpy array." )
    
    
        .def                                          ( "getReBoxedMap",
                                                        [] ( ProSHADE_run &self, proshade_unsign strNo ) -> pybind11::array_t < float >
                                                        {
                                                            //== Get the values
                                                            std::vector< proshade_signed > vals = self.getReBoxedBounds ( strNo );

                                                            //== Determine dimensions
                                                            proshade_unsign xDim = static_cast< proshade_unsign > ( vals.at(1) ) - static_cast< proshade_unsign > ( vals.at(0) ) + 1;
                                                            proshade_unsign yDim = static_cast< proshade_unsign > ( vals.at(3) ) - static_cast< proshade_unsign > ( vals.at(2) ) + 1;
                                                            proshade_unsign zDim = static_cast< proshade_unsign > ( vals.at(5) ) - static_cast< proshade_unsign > ( vals.at(4) ) + 1;
            
                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[xDim * yDim * zDim];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );

                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < (xDim * yDim * zDim); iter++ ) { npVals[iter] = static_cast< float > ( self.getMapValue ( strNo, iter ) ); }

                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleRebMap ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );

                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { xDim, yDim, zDim },                                                 // Shape
                                                                                                                            { yDim * zDim * sizeof(float), zDim * sizeof(float), sizeof(float) }, // C-stype strides
                                                                                                                            npVals,                                                               // Data
                                                                                                                            pyCapsuleRebMap );                                                    // Capsule

                                                            //== Done
                                                            return ( retArr );
                                                    }, "This function returns the re-boxed structure map as a numpy 3D array." )
    
        //============================================ Overlay results accessor functions as lambda functions directly returning numpy arrays
        .def                                          ( "getEulerAngles",
                                                        [] ( ProSHADE_run &self ) -> pybind11::array_t < float >
                                                        {
                                                            //== Get the values
                                                            std::vector< proshade_double > vals = self.getEulerAngles ( );
        
                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[static_cast<proshade_unsign> ( vals.size() )];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );

                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( vals.size() ); iter++ ) { npVals[iter] = static_cast< float > ( vals.at(iter) ); }

                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleEulAngs ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );

                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { static_cast<proshade_unsign> ( vals.size() ) },  // Shape
                                                                                                                            { sizeof(float) },                                 // C-stype strides
                                                                                                                            npVals,                                            // Data
                                                                                                                            pyCapsuleEulAngs );                                // Capsule

                                                            //== Done
                                                            return ( retArr );
                                                        }, "This function returns the vector of Euler angles with best overlay correlation." )
    
        .def                                          ( "getOptimalRotMat",
                                                        [] ( ProSHADE_run &self ) -> pybind11::array_t < float >
                                                        {
                                                            //== Get the values
                                                            std::vector< proshade_double > vals = self.getOptimalRotMat ( );
    
                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[static_cast<proshade_unsign> ( vals.size() )];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );

                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( vals.size() ); iter++ ) { npVals[iter] = static_cast< float > ( vals.at(iter) ); }

                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleRotMat ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );

                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { 3, 3 },                              // Shape
                                                                                                                            { 3 * sizeof(float), sizeof(float) },  // C-stype strides
                                                                                                                            npVals,                                // Data
                                                                                                                            pyCapsuleRotMat );                     // Capsule

                                                            //== Done
                                                            return ( retArr );
                                                        }, "This function returns the vector of Euler angles with best overlay correlation." )
    
        .def                                          ( "getTranslationToOrigin",
                                                        [] ( ProSHADE_run &self ) -> pybind11::array_t < float >
                                                        {
                                                            //== Get the values
                                                            std::vector< proshade_double > vals = self.getTranslationToOrigin ( );

                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[static_cast<proshade_unsign> ( vals.size() )];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );

                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( vals.size() ); iter++ ) { npVals[iter] = static_cast< float > ( vals.at(iter) ); }

                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleTTO ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );

                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { static_cast<proshade_unsign> ( vals.size() ) },  // Shape
                                                                                                                            { sizeof(float) },                                 // C-stype strides
                                                                                                                            npVals,                                            // Data
                                                                                                                            pyCapsuleTTO );                                    // Capsule

                                                            //== Done
                                                            return ( retArr );
                                                        }, "This function returns the negative values of the position of the rotation centre (the point about which the rotation should be done)." )
        .def                                          ( "getOriginToOverlayTranslation",
                                                        [] ( ProSHADE_run &self ) -> pybind11::array_t < float >
                                                        {
                                                            //== Get the values
                                                            std::vector< proshade_double > vals = self.getOriginToOverlayTranslation ( );

                                                            //== Allocate memory for the numpy values
                                                            float* npVals = new float[static_cast<proshade_unsign> ( vals.size() )];
                                                            ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );

                                                            //== Copy values
                                                            for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( vals.size() ); iter++ ) { npVals[iter] = static_cast< float > ( vals.at(iter) ); }

                                                            //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
                                                            pybind11::capsule pyCapsuleOTOT ( npVals, []( void *f ) { float* foo = reinterpret_cast< float* > ( f ); delete foo; } );

                                                            //== Copy the value
                                                            pybind11::array_t < float > retArr = pybind11::array_t<float> ( { static_cast<proshade_unsign> ( vals.size() ) },  // Shape
                                                                                                                            { sizeof(float) },                                 // C-stype strides
                                                                                                                            npVals,                                            // Data
                                                                                                                            pyCapsuleOTOT );                                   // Capsule

                                                            //== Done
                                                            return ( retArr );
                                                        }, "This function returns the translation required to move the structure from origin to optimal overlay." )
    
    
        //============================================ Description
        .def                                          ( "__repr__", [] ( ) { return "<ProSHADE_run class object> (Run class constructor takes a ProSHADE_settings object and completes a single run according to the settings object information)"; } );
}
