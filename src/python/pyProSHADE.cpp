/*! \file pyProSHADE.cpp
    \brief This file contains the PyBind11 bindings for the ProSHADE_settings class.
    
    ...
    
    Copyright by Michal Tykac and individual contributors. All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
 
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.5.0
    \date      DEC 2020
 */

//==================================================== Include header
#include "ProSHADE_settings.hpp"

//==================================================== Include full ProSHADE (including cpp files looks horrible, but it is the only way I can find to stop PyBind11 from complaining)
#include "ProSHADE_misc.cpp"
#include "ProSHADE_maths.cpp"
#include "ProSHADE_tasks.cpp"
#include "ProSHADE_io.cpp"
#include "ProSHADE_data.cpp"
#include "ProSHADE_symmetry.cpp"
#include "ProSHADE_overlay.cpp"
#include "ProSHADE_wignerMatrices.cpp"
#include "ProSHADE_spheres.cpp"
#include "ProSHADE_mapManip.cpp"
#include "ProSHADE_messages.cpp"
#include "ProSHADE_distances.cpp"
#include "ProSHADE_peakSearch.cpp"
#include "ProSHADE_sphericalHarmonics.cpp"
#include "ProSHADE.cpp"

//==================================================== Include PyBind11 header
#include <pybind11/pybind11.h>

void add_settingsClass ( pybind11::module& pyProSHADE )
{
    //================================================ Export the ProSHADE_settings class
    pybind11::class_ < ProSHADE_settings >            ( pyProSHADE, "ProSHADE_settings" )
        
        // Constructors / Destructors
        .def                                          ( pybind11::init<>() )
    
        // Member variables
        .def_readwrite                                ( "task",                                 &ProSHADE_settings::task                                )
                                
        .def_readwrite                                ( "inputFiles",                           &ProSHADE_settings::inputFiles                          )
        .def_readwrite                                ( "forceP1",                              &ProSHADE_settings::forceP1                             )
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
    
        .def_readwrite                                ( "reBoxMap",                             &ProSHADE_settings::reBoxMap                            )
        .def_readwrite                                ( "boundsExtraSpace",                     &ProSHADE_settings::boundsExtraSpace                    )
        .def_readwrite                                ( "boundsSimilarityThreshold",            &ProSHADE_settings::boundsSimilarityThreshold           )
        .def_readwrite                                ( "useSameBounds",                        &ProSHADE_settings::useSameBounds                       )
        .def_readwrite                                ( "forceBounds",                          &ProSHADE_settings::forceBounds                         )
    
        .def_readwrite                                ( "moveToCOM",                            &ProSHADE_settings::moveToCOM                           )
    
        .def_readwrite                                ( "addExtraSpace",                        &ProSHADE_settings::addExtraSpace                       )
    
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
        .def_readwrite                                ( "usePeakSearchInRotationFunctionSpace", &ProSHADE_settings::usePeakSearchInRotationFunctionSpace)
        .def_readwrite                                ( "useBiCubicInterpolationOnPeaks",       &ProSHADE_settings::useBiCubicInterpolationOnPeaks      )
        .def_readwrite                                ( "maxSymmetryFold",                      &ProSHADE_settings::maxSymmetryFold                     )
    
        .def_readwrite                                ( "overlayStructureName",                 &ProSHADE_settings::overlayStructureName                )
        .def_readwrite                                ( "rotTrsJSONFile",                       &ProSHADE_settings::rotTrsJSONFile                      )
    
        .def_readwrite                                ( "verbose",                              &ProSHADE_settings::verbose                             )
    
        // Description
        .def                                          ( "__repr__", [] ( const ProSHADE_settings &a ) { return "<ProSHADE_settings class object> (Settings class is used to set all settings values in a single place)"; } );
}



