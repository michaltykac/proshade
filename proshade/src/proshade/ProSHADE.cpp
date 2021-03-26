/*! \file ProSHADE.cpp
    \brief This is the main source file providing the main access class and its functions.
    
    This file contains the definitions for the main access class (ProSHADE_run), which the user can use to
    start a run of the ProSHADE tool. This is generally done by firstly creating and setting a ProSHADE_settings class
    instance, which is then supplied to the constructor of the ProSHADE_run class defined here. Once the class constructor
    is run, ProSHADE run with the supplied settings will also be complete and the resulting instance of ProSHADE_run class
    will contain the results. To access these results, the ProSHADE_run class provides accessor functions.
    
    Copyright by Michal Tykac and individual contributors. All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
 
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.5.4
    \date      MAR 2021
 */

//==================================================== ProSHADE
#include "ProSHADE.hpp"

/*! \brief Contructor for the ProSHADE_settings class.
 
    This is the generic constructor used in cases where the settings object will be filled based on run-time determined values. If you know the specific
    task to be done, it is recommended to use the constructor which takes the task as argument, so that the default values are set specifically for the
    task at hand.
 */
ProSHADE_settings::ProSHADE_settings ( )
{
    //================================================ Settings regarding the task at hand
    this->task                                        = NA;
    
    //================================================ Settings regarding input files
    this->forceP1                                     = true;
    this->removeWaters                                = true;
    this->firstModelOnly                              = true;
    
    //================================================ Settings regarding the resolution of calculations
    this->requestedResolution                         = -1.0;
    this->changeMapResolution                         = false;
    this->changeMapResolutionTriLinear                = false;
    
    //================================================ Settings regarding the PDB B-factor change
    this->pdbBFactorNewVal                            = -1.0;
    
    //================================================ Settings regarding the bandwidth of calculations
    this->maxBandwidth                                = 0;
    this->rotationUncertainty                         = 0;
    
    //================================================ Settings regarding the phase
    this->usePhase                                    = true;
    
    //================================================ Settings regarding the spheres
    this->maxSphereDists                              = 0.0;
    
    //================================================ Settings regarding the Gauss-Legendre integration
    this->integOrder                                  = 0;
    this->taylorSeriesCap                             = 10;
    
    //================================================ Settings regarding map normalisation
    this->normaliseMap                                = false;
    
    //================================================ Settings regarding map inversion
    this->invertMap                                   = false;
    
    //================================================ Settings regarding map masking
    this->blurFactor                                  = 350.0;
    this->maskingThresholdIQRs                        = 3.0;
    this->maskMap                                     = false;
    this->useCorrelationMasking                       = false;
    this->halfMapKernel                               = 0.0;
    this->correlationKernel                           = 0.0;
    this->saveMask                                    = false;
    this->maskFileName                                = "maskFile";
    
    //================================================ Settings regarding re-boxing
    this->reBoxMap                                    = false;
    this->boundsExtraSpace                            = 3.0;
    this->boundsSimilarityThreshold                   = 0;
    this->useSameBounds                               = false;
    this->forceBounds                                 = new proshade_signed [6];
    
    //================================================ Settings regarding COM
    this->moveToCOM                                   = false;
    
    //================================================ Settings regarding extra cell space
    this->addExtraSpace                               = 10.0;
    
    //================================================ Settings regarding shell settings
    this->progressiveSphereMapping                    = false;
    
    //================================================ Settings regarding output file name
    this->outName                                     = "reBoxed";
    
    //================================================ Settings regarding distances computation
    this->computeEnergyLevelsDesc                     = true;
    this->computeTraceSigmaDesc                       = true;
    this->computeRotationFuncDesc                     = true;
    this->enLevMatrixPowerWeight                      = 1.0;
    
    //================================================ Settings regarding peak searching
    this->peakNeighbours                              = 1;
    this->noIQRsFromMedianNaivePeak                   = -999.9;
    
    //================================================ Settings regarding 1D grouping
    this->smoothingFactor                             =  15.0;
    
    //================================================ Settings regarding the symmetry detection
    this->usePeakSearchInRotationFunctionSpace        = true;
    this->useBiCubicInterpolationOnPeaks              = true;
    this->maxSymmetryFold                             = 30;
    this->symMissPeakThres                            = 0.3;
    this->axisErrTolerance                            = 0.01;
    this->axisErrToleranceDefault                     = true;
    this->minSymPeak                                  = 0.3;
    this->recommendedSymmetryType                     = "";
    this->recommendedSymmetryFold                     = 0;
    this->requestedSymmetryType                       = "";
    this->requestedSymmetryFold                       = 0;
    this->detectedSymmetry.clear                      ( );
    
    //================================================ Settings regarding the structure overlay
    this->overlayStructureName                        = "movedStructure";
    this->rotTrsJSONFile                              = "movedStructureOperations.json";
    
    //================================================ Settings regarding verbosity of the program
    this->verbose                                     = 1;
    
    //================================================ Done
    
}

/*! \brief Contructor for the ProSHADE_settings class for particular task.
 
    This is the generic constructor used in cases where the settings object will be filled based on run-time determined values. If you know the specific
    task to be done, it is recommended to use the constructor which takes the task as argument, so that the default values are set specifically for the
    task at hand.
 
    \param[in] taskToPerform The task that should be performed by ProSHADE.
 */
ProSHADE_settings::ProSHADE_settings ( ProSHADE_Task taskToPerform )
{
    //================================================ Settings regarding the task at hand
    this->task                                        = taskToPerform;
    
    //================================================ Settings regarding input files
    this->forceP1                                     = true;
    this->removeWaters                                = true;
    this->firstModelOnly                              = true;
    
    //================================================ Settings regarding the resolution of calculations
    this->requestedResolution                         = -1.0;
    this->changeMapResolution                         = false;
    this->changeMapResolutionTriLinear                = false;
    
    //================================================ Settings regarding the PDB B-factor change
    this->pdbBFactorNewVal                            = -1.0;
    
    //================================================ Settings regarding the bandwidth of calculations
    this->maxBandwidth                                = 0;
    this->rotationUncertainty                         = 0;
    
    //================================================ Settings regarding the phase
    this->usePhase                                    = true;
    
    //================================================ Settings regarding the spheres
    this->maxSphereDists                              = 0.0;
    
    //================================================ Settings regarding the Gauss-Legendre integration
    this->integOrder                                  = 0;
    this->taylorSeriesCap                             = 10;
    
    //================================================ Settings regarding map normalisation
    this->normaliseMap                                = false;
    
    //================================================ Settings regarding map inversion
    this->invertMap                                   = false;
    
    //================================================ Settings regarding map masking
    this->blurFactor                                  = 350.0;
    this->maskingThresholdIQRs                        = 3.0;
    this->maskMap                                     = false;
    this->useCorrelationMasking                       = false;
    this->halfMapKernel                               = 0.0;
    this->correlationKernel                           = 0.0;
    this->saveMask                                    = false;
    this->maskFileName                                = "maskFile";
    this->detectedSymmetry.clear                      ( );
    
    //================================================ Settings regarding re-boxing
    this->reBoxMap                                    = false;
    this->boundsExtraSpace                            = 3.0;
    this->boundsSimilarityThreshold                   = 0;
    this->useSameBounds                               = false;
    this->forceBounds                                 = new proshade_signed [6];
    
    //================================================ Settings regarding extra cell space
    this->addExtraSpace                               = 10.0;
    
    //================================================ Settings regarding shell settings
    this->progressiveSphereMapping                    = false;
    
    //================================================ Settings regarding output file name
    this->outName                                     = "reBoxed";
    
    //================================================ Settings regarding distances computation
    this->computeEnergyLevelsDesc                     = true;
    this->computeTraceSigmaDesc                       = true;
    this->computeRotationFuncDesc                     = true;
    this->enLevMatrixPowerWeight                      = 1.0;
    
    //================================================ Settings regarding peak searching
    this->peakNeighbours                              = 1;
    this->noIQRsFromMedianNaivePeak                   = -999.9;
    
    //================================================ Settings regarding 1D grouping
    this->smoothingFactor                             =  15.0;
    
    //================================================ Settings regarding the symmetry detection
    this->usePeakSearchInRotationFunctionSpace        = true;
    this->useBiCubicInterpolationOnPeaks              = true;
    this->maxSymmetryFold                             = 30;
    this->symMissPeakThres                            = 0.3;
    this->axisErrTolerance                            = 0.01;
    this->axisErrToleranceDefault                     = true;
    this->minSymPeak                                  = 0.3;
    this->recommendedSymmetryType                     = "";
    this->recommendedSymmetryFold                     = 0;
    this->requestedSymmetryType                       = "";
    this->requestedSymmetryFold                       = 0;
    
    //================================================ Settings regarding the structure overlay
    this->overlayStructureName                        = "movedStructure";
    this->rotTrsJSONFile                              = "movedStructureOperations.json";
    
    //================================================ Settings regarding verbosity of the program
    this->verbose                                     = 1;
    
    //================================================ Task specific settings
    switch ( this->task )
    {
        case NA:
            std::cerr << std::endl << "=====================" << std::endl << "!! ProSHADE ERROR !!" << std::endl << "=====================" << std::endl << std::flush;
            std::cerr << "Error Code          : " << "E000014" << std::endl << std::flush;
            std::cerr << "ProSHADE version    : " << __PROSHADE_VERSION__ << std::endl << std::flush;
            std::cerr << "File                : " << "ProSHADE.cpp" << std::endl << std::flush;
            std::cerr << "Line                : " << 97 << std::endl << std::flush;
            std::cerr << "Function            : " << "ProSHADE_settings (Task) constructor" << std::endl << std::flush;
            std::cerr << "Message             : " << "No task has been specified for task specific constructor." << std::endl << std::flush;
            std::cerr << "Further information : " << "This ProSHADE_settings class constructor is intended to\n                    : set the internal variables to default value given a\n                    : particular taks. By supplying this task as NA, this beats\n                    : the purpose of the constructor. Please use the\n                    : non-argumental constructor if task is not yet known." << std::endl << std::endl << std::flush;
            ProSHADE_internal_messages::printTerminateMessage ( this->verbose );
            exit                                      ( EXIT_FAILURE );
            break;
            
        case Symmetry:
            this->requestedResolution                 = 6.0;
            this->pdbBFactorNewVal                    = 80.0;
            this->changeMapResolution                 = true;
            this->maskMap                             = false;
            this->moveToCOM                           = true;
            this->normaliseMap                        = true;
            this->reBoxMap                            = false;
            break;
            
        case Distances:
            this->changeMapResolution                 = false;
            this->maskMap                             = false;
            this->moveToCOM                           = true;
            this->reBoxMap                            = false;
            break;
                    
        case OverlayMap:
            this->requestedResolution                 = 8.0;
            this->changeMapResolution                 = true;
            this->maskMap                             = false;
            this->moveToCOM                           = false;
            this->normaliseMap                        = false;
            this->reBoxMap                            = false;
            break;
                    
        case MapManip:
            this->changeMapResolution                 = false;
            this->maskMap                             = true;
            this->moveToCOM                           = false;
            break;
    }
    
    //================================================ Done
    
}

/*! \brief Destructor for the ProSHADE_settings class.
 
    This destructor is responsible for releasing all memory used by the settings object
 */
ProSHADE_settings::~ProSHADE_settings ( void )
{
    //================================================ Release boundaries variable
    delete[] this->forceBounds;
    
    //================================================ Release symmetry axes
    if ( this->detectedSymmetry.size() > 0 ) { for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( this->detectedSymmetry.size() ); it++ ) { if ( this->detectedSymmetry.at(it) != NULL ) { delete[] this->detectedSymmetry.at(it); } } }
    
    //================================================ Done
    
}

/*! \brief Function to determine general values that the user left on auto-determination.
 */
void ProSHADE_settings::setVariablesLeftOnAuto ( void  )
{
    //================================================ Determine the peak IQR from median threshold, unless given by user
    if ( this->noIQRsFromMedianNaivePeak == -999.9 )
    {
        //============================================ If using the old symmetry detection algorithm or distances computation, this will be used on many small peaks with few outliers. Use value of 5.0
        if (   this->task == Distances )                                                      { this->noIQRsFromMedianNaivePeak = 5.0; }
        if ( ( this->task == Symmetry  ) && ( !this->usePeakSearchInRotationFunctionSpace ) ) { this->noIQRsFromMedianNaivePeak = 5.0; }
        
        //============================================ If using the new symmetry detection algorithm, this needs to be decreasing with resolution. How much, that is a bit arbitrary...
        if ( ( this->task == Symmetry  ) && (  this->usePeakSearchInRotationFunctionSpace ) ) { this->noIQRsFromMedianNaivePeak = std::max ( 0.0, 1.0 - ( this->requestedResolution * 0.05 ) ); }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Adds a structure file name to the appropriate variable.
 
    This function takes a string defining the filename of a structure to be processed and adds it to
    the list of structures to be processed.
 
    \param[in] structure String file name to be added to the list of structures to process.
 */
void ProSHADE_settings::addStructure ( std::string structure )
{
    //================================================ Use C++ version independent vector processing
    ProSHADE_internal_misc::addToStringVector         ( &( this->inputFiles ), structure );
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested resolution in the appropriate variable.
 
    This function sets the resolution in the appropriate variable.
 
    \param[in] resolution The requested value for the resolution to which the computations are to be done.
 */
void ProSHADE_settings::setResolution ( proshade_single resolution )
{
    //================================================ Set the value
    this->requestedResolution                         = resolution;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested B-factor value for PDB files in the appropriate variable.
 
    This function sets the B-factor value for PDB files in the appropriate variable.
 
    \param[in] newBF The requested value for the B-factor value for PDB files for smooth and processible maps.
 */
void ProSHADE_settings::setPDBBFactor ( proshade_double newBF )
{
    //================================================ Set the value
    this->pdbBFactorNewVal                            = newBF;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested map normalisation value in the appropriate variable.
 
    This function sets the map normalisation between on and off.
 
    \param[in] normalise The requested value for the map normalisation (on = true, off = false).
 */
void ProSHADE_settings::setNormalisation ( bool normalise )
{
    //================================================ Set the value
    this->normaliseMap                                = normalise;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested map inversion value in the appropriate variable.
 
    This function sets the map inversion between on and off.
 
    \param[in] mInv Should the map be inverted? (on = true, off = false).
 */
void ProSHADE_settings::setMapInversion ( bool mInv )
{
    //================================================ Set the value
    this->invertMap                                   = mInv;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested verbosity in the appropriate variable.
 
    This function sets the varbosity of the ProSHADE run in the appropriate variable.
 
    \param[in] verbose The requested value for verbosity. -1 means no output, while 4 means loud output
 */
void ProSHADE_settings::setVerbosity ( proshade_signed verbosity )
{
    //================================================ Set the value
    this->verbose                                     = verbosity;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested map blurring factor in the appropriate variable.
 
    This function sets the blurring / sharpening factor for map masking in the appropriate variable.
 
    \param[in] blurFac The requested value for the blurring factor.
 */
void ProSHADE_settings::setMaskBlurFactor ( proshade_single blurFac )
{
    //================================================ Set the value
    this->blurFactor                                  = blurFac;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested number of IQRs for masking threshold in the appropriate variable.
 
    This function sets the number of interquartile ranges from the median to be used for map masking in the correct
    variable.
 
    \param[in] noIQRs The requested value for the number of IQRs from the median to be used for masking threshold.
 */
void ProSHADE_settings::setMaskIQR ( proshade_single noIQRs )
{
    //================================================ Set the value
    this->maskingThresholdIQRs                        = noIQRs;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested map masking decision value in the appropriate variable.
 
    This function sets the map masking between on and off.
 
    \param[in] mask The requested value for the map masking (on = true, off = false).
 */
void ProSHADE_settings::setMasking ( bool mask )
{
    //================================================ Set the value
    this->maskMap                                     = mask;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested map masking type in the appropriate variable.
 
    This function sets the map masking type. If false, the standard map blurring masking will be used, while
    if true, the new "fake" half-map correlation mask will be used.
 
    \param[in] corMask The requested value for the map masking type.
 */
void ProSHADE_settings::setCorrelationMasking ( bool corMask )
{
    //================================================ Set the value
    this->useCorrelationMasking                       = corMask;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested "fake" half-map kernel in the appropriate variable.
 
    This function sets the kernel for creating the "fake" half-map. What is meant here is that a new map is
    created as the average of neighbours from the original map - this is useful in masking. This value then
    sets how many neighbours.
 
    \param[in] typNoi The requested value for the typical noise size in Angstrom.
 */
void ProSHADE_settings::setTypicalNoiseSize ( proshade_single typNoi )
{
    //================================================ Set the value
    this->halfMapKernel                               = typNoi;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested minimum mask size.
 
    This function sets the kernel for the local correlation computation between the "fake half-map" and the original map.
 
    \param[in] minMS The requested value for the minimum mask size in Angstrom.
 */
void ProSHADE_settings::setMinimumMaskSize ( proshade_single minMS )
{
    //================================================ Set the value
    this->correlationKernel                           = minMS;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets whether the mask should be saved.
 
    This function sets the switch variable to whether mask should be saved.
 
    \param[in] savMsk If true, mask will be saved, otherwise it will not be.
 */
void ProSHADE_settings::setMaskSaving ( bool savMsk )
{
    //================================================ Set the value
    this->saveMask                                    = savMsk;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets where the mask should be saved.
 
    This function sets the the filename to which mask should be saved.
 
    \param[in] mskFln The filename where the mask should be saved.
 */
void ProSHADE_settings::setMaskFilename ( std::string mskFln )
{
    //================================================ Set the value
    this->maskFileName                                = mskFln;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets whether re-boxing needs to be done in the appropriate variable.
 
    This function sets the switch as to whether re-boxing needs to be done to the correct variable.
 
    \param[in] reBx The requested value for the re-boxing switch variable.
 */
void ProSHADE_settings::setMapReboxing ( bool reBx )
{
    //================================================ Set the value
    this->reBoxMap                                    = reBx;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested number of angstroms for extra space in re-boxing in the appropriate variable.
 
    This function sets the number of angstroms to be added both before and after the absolute bounds for re-boxing to
    the correct variable.
 
    \param[in] boundsExSp The requested value for the extra re-boxing space in angstroms.
 */
void ProSHADE_settings::setBoundsSpace ( proshade_single boundsExSp )
{
    //================================================ Set the value
    this->boundsExtraSpace                            = boundsExSp;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the threshold for number of indices difference acceptable to make index sizes same in the appropriate variable.
 
    This function sets the number of indices by which two dimensions can differ for them to be still made the same size.
 
    \param[in] boundsThres The requested value for the bouds difference threhshold.
 */
void ProSHADE_settings::setBoundsThreshold ( proshade_signed boundsThres )
{
    //================================================ Set the value
    this->boundsSimilarityThreshold                   = boundsThres;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets whether same boundaries should be used in the appropriate variable.
 
    This function sets the switch as to whether the same boundaries as for the first map should be forced upon the rest
    if the input maps.
 
    \param[in] sameB The requested value for the same boundaries as first structure switch variable.
 */
void ProSHADE_settings::setSameBoundaries ( bool sameB )
{
    //================================================ Set the value
    this->useSameBounds                               = sameB;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested output file name in the appropriate variable.
 
    This function sets the filename to which the output structure(s) should be saved. This variable is used by multiple tasks
    and therefore cannot be more specifically described here.
 
    \param[in] oFileName The requested value for the output file name variable.
 */
void ProSHADE_settings::setOutputFilename ( std::string oFileName )
{
    //================================================ Set the value
    this->outName                                     = oFileName;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested map resolution change decision in the appropriate variable.
 
    This function sets the map resolution change between on and off.
 
    \param[in] mrChange The requested value for the map resolution change (on = true, off = false).
 */
void ProSHADE_settings::setMapResolutionChange ( bool mrChange )
{
    //================================================ Set the value
    this->changeMapResolution                         = mrChange;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested map resolution change decision using tri-linear interpolation in the appropriate variable.
 
    This function sets the tri-linear interpolation map resolution change between on and off.
 
    \param[in] mrChange The requested value for the map resolution change (on = true, off = false).
 */
void ProSHADE_settings::setMapResolutionChangeTriLinear ( bool mrChange )
{
    //================================================ Set the value
    this->changeMapResolutionTriLinear                = mrChange;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested map centering decision value in the appropriate variable.
 
    This function sets the map centering using COM between on and off.
 
    \param[in] com The requested value for the map centering (on = true, off = false).
 */
void ProSHADE_settings::setMapCentering ( bool com )
{
    //================================================ Set the value
    this->moveToCOM                                   = com;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested map extra space value in the appropriate variable.
 
    This function sets the amount of extra space to be added to internal maps in the appropriate variable.
 
    \param[in] exSpace The requested amount of extra space.
 */
void ProSHADE_settings::setExtraSpace ( proshade_single exSpace )
{
    //================================================ Set the value
    this->addExtraSpace                               = exSpace;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested sphere mapping value settings approach in the appropriate variable.
 
    This function sets the progressive sphere mapping approach between on and off.
 
    \param[in] com The requested value for the progressive sphere mapping (on = true, off = false).
 */
void ProSHADE_settings::setProgressiveSphereMapping ( bool progSphMap )
{
    //================================================ Set the value
    this->progressiveSphereMapping                    = progSphMap;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested spherical harmonics bandwidth in the appropriate variable.
 
    This function sets the bandwidth limit for the spherical harmonics computations in the appropriate variable.
 
    \param[in] band The requested value for spherical harmonics bandwidth (0 = AUTOMATIC DETERMINATION).
 */
void ProSHADE_settings::setBandwidth ( proshade_unsign band )
{
    //================================================ Set the value
    this->maxBandwidth                                = band;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested distance between spheres in the appropriate variable.
 
    This function sets the distance between any two consecutive spheres in the sphere mapping of a map in the appropriate variable.
 
    \param[in] sphDist The requested value for distance between spheres (0 = AUTOMATIC DETERMINATION).
 */
void ProSHADE_settings::setSphereDistances ( proshade_single sphDist )
{
    //================================================ Set the value
    this->maxSphereDists                              = sphDist;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested order for the Gauss-Legendre integration in the appropriate variable.
 
    This function sets the Gauss-Legendre integration between the spheres order value in the appropriate variable.
 
    \param[in] intOrd The requested value for the Gauss-Legendre integration order (0 = AUTOMATIC DETERMINATION).
 */
void ProSHADE_settings::setIntegrationOrder ( proshade_unsign intOrd )
{
    //================================================ Set the value
    this->integOrder                                  = intOrd;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the requested Taylor series cap for the Gauss-Legendre integration in the appropriate variable.
 
    This function sets the Taylor series maximum limit for the Gauss-Legendre integration between the spheres order
    value in the appropriate variable.
 
    \param[in] tayCap The requested value for the Taylor series cap. (0 = AUTOMATIC DETERMINATION).
 */
void ProSHADE_settings::setTaylorSeriesCap ( proshade_unsign tayCap )
{
    //================================================ Set the value
    this->taylorSeriesCap                             = tayCap;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets whether the energy level distance descriptor should be computed.
 
    This function sets the boolean variable deciding whether the RRP matrices and the energy levels descriptor should
    be computed or not.
 
    \param[in] enLevDesc The requested value for the energy levels descriptor computation switch.
 */
void ProSHADE_settings::setEnergyLevelsComputation ( bool enLevDesc )
{
    //======================================== Set the value
    this->computeEnergyLevelsDesc             = enLevDesc;
    
    //======================================== Done
    return ;
    
}

/*! \brief Sets whether the trace sigma distance descriptor should be computed.
 
    This function sets the boolean variable deciding whether the E matrices and the trace sigma descriptor should
    be computed or not.
 
    \param[in] trSigVal The requested value for the trace sigma descriptor computation switch.
 */
void ProSHADE_settings::setTraceSigmaComputation ( bool trSigVal )
{
    //================================================ Set the value
    this->computeTraceSigmaDesc                       = trSigVal;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets whether the rotation function distance descriptor should be computed.
 
    This function sets the boolean variable deciding whether the inverse SO(3) transform and the rotation function descriptor should
    be computed or not.
 
    \param[in] rotfVal The requested value for the rotation function descriptor computation switch.
 */
void ProSHADE_settings::setRotationFunctionComputation  ( bool rotfVal )
{
    //================================================ Set the value
    this->computeRotationFuncDesc                     = rotfVal;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the number of neighbour values that have to be smaller for an index to be considered a peak.
 
    This function sets the number of neighbouring points (in all three dimensions and both positive and negative direction) that
    have to have lower value than the currently considered index in order for this index to be considered as a peak.
 
    \param[in] pkS The requested value for the number of neighbours being lower for a peak.
 */
void ProSHADE_settings::setPeakNeighboursNumber ( proshade_unsign pkS )
{
    //================================================ Set the value
    this->peakNeighbours                              = pkS;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the number of IQRs from the median for threshold height a peak needs to be considered a peak.
 
    This function sets the number of IQRs from the median that determine the threshold used to determine if a 'naive' peak
    is a peak, or just a random local maxim in the background. The set from which median and IQR is computed is the non-peak
    values.
 
    \param[in] noIQRs The requested number of IQRs from the median.
 */
void ProSHADE_settings::setPeakNaiveNoIQR ( proshade_double noIQRs )
{
    //================================================ Set the value
    this->noIQRsFromMedianNaivePeak                   = noIQRs;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets whether the phase information will be used.
 
    This function sets the boolean variable deciding whether the phase information should be used. If not, Patterson maps
    will be used instead of density maps and the 3D data will be converted to them. Also, only even bands of the spherical
    harmonics decomposition will be computed as the odd bands must be 0.
 
    \param[in] phaseUsage The requested value for the phase usage switch.
 */
void ProSHADE_settings::setPhaseUsage ( bool phaseUsage )
{
    //================================================ Set the value
    this->usePhase                                    = phaseUsage;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the weight of shell position for the energy levels computation.
 
    During the computation of the energy levels descriptor, Pearson's correlation coefficient is computed between different shells
    with the same band. The shell index can by expanded to its mPower exponential to give higher shells more weight, or vice versa.
    To do this, set the mPower value as you see fit.
 
    \param[in] mPower The requested value for the shell position exponential.
 */
void ProSHADE_settings::setEnLevShellWeight ( proshade_double mPower )
{
    //================================================ Set the value
    this->enLevMatrixPowerWeight                      = mPower;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the grouping smoothing factor into the proper variable.
 
    When detecting symmetry, it is worth grouping the possible rotations by their self-rotation function peak heights. In this
    process, the distribution of peak heights needs to be smoothen over and this factor decides how smooth it should be. Small
    value leads to all peaks being in the same group, while large number means each peak will be in its own group.
 
    \param[in] smFact The requested value for the grouping smoothing factor.
 */
void ProSHADE_settings::setGroupingSmoothingFactor ( proshade_double smFact )
{
    //================================================ Set the value
    this->smoothingFactor                             = smFact;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the threshold for starting the missing peaks procedure.
 
    When only mpThres percentage of peaks are missing during symmetry detection, the full missing peak detection procedure will
    be started. Otherwise, the symmetry will not be detected at all.
 
    \param[in] mpThres The requested value for the missing peaks procedure starting threshold.
 */
void ProSHADE_settings::setMissingPeakThreshold ( proshade_double mpThres )
{
    //================================================ Set the value
    this->symMissPeakThres                            = mpThres;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the threshold for matching symmetry axes.
 
    When comparing symmetry axes, there needs to be a threshold allowing for some small error comming from the numberical
    inaccuracies. This is where you set this threshold.
 
    \param[in] axThres The requested value for the axes comparison threshold.
 */
void ProSHADE_settings::setAxisComparisonThreshold ( proshade_double axThres )
{
    //================================================ Set the value
    this->axisErrTolerance                            = axThres;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the automatic symmetry axis tolerance decreasing.

    When comparing symmetry axes, there needs to be a threshold allowing for some small error comming from the numberical
    inaccuracies. It turns out that this threshold should take into account the ratio to the next symmetry angles, otherwise it would
    strongly prefer larger symmetries. This variable decides whether the threshold should be decreased based on the fold of sought
    åsymmetry or not.

    \param[in] behav The requested value for the axes comparison threshold decreasing.
*/
void ProSHADE_settings::setAxisComparisonThresholdBehaviour ( bool behav )
{
    //================================================ Set the value
    this->axisErrToleranceDefault                     = behav;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the minimum peak height for symmetry axis to be considered.
 
    When considering if a symmetry axis is "real" and should be acted upon, its average peak height will need to
    be higher than this value.
 
    \param[in] minSP The requested value for the minimum peak height.
 */
void ProSHADE_settings::setMinimumPeakForAxis ( proshade_double minSP )
{
    //================================================ Set the value
    this->minSymPeak                                  = minSP;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the ProSHADE detected symmetry type.
 
    When symmetry detection is done, the resulting recommended symmetry type will be saved in the settings object by this function.
 
    \param[in] val The recommended symmetry type for the structure.
 
    \warning This is an internal function and it should not be used by the user.
 */
void ProSHADE_settings::setRecommendedSymmetry ( std::string val )
{
    //================================================ Set the value
    this->recommendedSymmetryType                     = val;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the ProSHADE detected symmetry fold.
 
    When symmetry detection is done, the resulting recommended symmetry fold  (valid only for C and D symmetry types) will be saved in the
    settings object by this function.
 
    \param[in] val The recommended symmetry fold for the structure.
 
    \warning This is an internal function and it should not be used by the user.
 */
void ProSHADE_settings::setRecommendedFold ( proshade_unsign val )
{
    //================================================ Set the value
    this->recommendedSymmetryFold                     = val;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the user requested symmetry type.
 
    When symmetry detection is started, this symmetry type will be exclusively sought.
 
    \param[in] val The requested symmetry type for the structure.
 */
void ProSHADE_settings::setRequestedSymmetry ( std::string val )
{
    //================================================ Set the value
    this->requestedSymmetryType                       = val;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the user requested symmetry fold.
 
    When symmetry detection is started, this symmetry fold will be exclusively sought.
 
    \param[in] val The requested symmetry fold for the structure.
 */
void ProSHADE_settings::setRequestedFold ( proshade_unsign val )
{
    //================================================ Set the value
    this->requestedSymmetryFold                       = val;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the final detected symmetry axes information.
    
    This function copies (deep copy) the detected and recommended (or requested) symmetry axis information into the settings
    object variable for further processing. For multiple axes, call this function multiple times - the addition is cumulative.
 
    \param[in] sym A pointer to single symmetry axis constituting the detected symmetry.
 */
void ProSHADE_settings::setDetectedSymmetry ( proshade_double* sym )
{
    //================================================ Allocate memory
    proshade_double* hlpAxis                          = new proshade_double [6];
    ProSHADE_internal_misc::checkMemoryAllocation     ( hlpAxis, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy (deep) data
    hlpAxis[0]                                        = sym[0];
    hlpAxis[1]                                        = sym[1];
    hlpAxis[2]                                        = sym[2];
    hlpAxis[3]                                        = sym[3];
    hlpAxis[4]                                        = sym[4];
    hlpAxis[5]                                        = sym[5];
    
    //================================================ Save
    ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( &this->detectedSymmetry, hlpAxis );
    
    //================================================ Release memory
    delete[] hlpAxis;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the filename to which the overlay structure is to be save into.
 
    \param[in] filename The filename to which the overlay structure is to be saved to.
 */
void ProSHADE_settings::setOverlaySaveFile ( std::string filename )
{
    //================================================ Set the value
    this->overlayStructureName                        = filename;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the filename to which the overlay operations are to be save into.
 
    \param[in] filename The filename to which the overlay operations are to be saved to.
 */
void ProSHADE_settings::setOverlayJsonFile ( std::string filename )
{
    //================================================ Set the value
    this->rotTrsJSONFile                              = filename;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the symmetry detection algorithm type.
 
    \param[in] rotFunPeaks Should the original peak detection in rotation function space be used (FALSE), or should the new angle-axis space search be used (DEFAULT - TRUE)?
 */
void ProSHADE_settings::setSymmetryRotFunPeaks ( bool rotFunPeaks )
{
    //================================================ Set the value
    this->usePeakSearchInRotationFunctionSpace        = rotFunPeaks;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the bicubic interpolation on peaks.
 
    \param[in] bicubPeaks Should bicubic interpolation be done to search for improved axis in between peak index values (DEFAULT - TRUE)?
 */
void ProSHADE_settings::setBicubicInterpolationSearch ( bool bicubPeaks )
{
    //================================================ Set the value
    this->useBiCubicInterpolationOnPeaks              = bicubPeaks;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the maximum symmetry fold (well, the maximum prime symmetry fold).
 
    \param[in] maxFold Maximum prime number fold that will be searched for. Still its multiples may also be found.
 */
void ProSHADE_settings::setMaxSymmetryFold ( proshade_unsign maxFold )
{
    //================================================ Set the value
    this->maxSymmetryFold                             = maxFold;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function determines the bandwidth for the spherical harmonics computation.
 
    This function is here to automstically determine the bandwidth to which the spherical harmonics computations should be done.
    It accomplishes this by checking if value is already set, and if not (value is 0), then it sets it to half of the maximum
    circumference of the map, in indices as recommended by Kostelec and Rockmore (2007).
 
    \param[in] circumference The maximum circumference of the map.
 */
void ProSHADE_settings::determineBandwidth ( proshade_unsign circumference )
{
    //================================================ Check the current settings value is set to auto
    if ( this->maxBandwidth != 0 )
    {
        std::stringstream hlpSS;
        hlpSS << "The bandwidth was determined as: " << this->maxBandwidth;
        ProSHADE_internal_messages::printProgressMessage ( this->verbose, 3, hlpSS.str() );
        return ;
    }
    
    //================================================ Determine automatically
    this->maxBandwidth                                = ProSHADE_internal_spheres::autoDetermineBandwidth ( circumference );
    
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "The bandwidth was determined as: " << this->maxBandwidth;
    ProSHADE_internal_messages::printProgressMessage ( this->verbose, 3, hlpSS.str() );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function determines the bandwidth for the spherical harmonics computation from the allowed rotation function angle uncertainty.
 
    This function makes use of the fact that the rotation function dimensions will be 2 * bandwidth and that the dimensions will be covering full
    360 degrees rotation space. Therefore, by saying what is the maximum allowed angle uncertainty, the minimum required bandwidth value can be
    determined.
 
 \param[in] uncertainty The maximum allowed uncertainty on the rotation function.
 */
void ProSHADE_settings::determineBandwidthFromAngle ( proshade_double uncertainty )
{
    //================================================ Determine bandwidth
    if ( static_cast<proshade_unsign> ( std::ceil ( ( 360.0 / uncertainty ) / 2 ) ) % 2 == 0 )
    {
        this->maxBandwidth                            = static_cast<proshade_unsign> ( std::ceil ( ( 360.0 / uncertainty ) / 2.0 ) );
    }
    else
    {
        this->maxBandwidth                            = static_cast<proshade_unsign> ( std::ceil ( ( 360.0 / uncertainty ) / 2.0 ) ) + 1;
    }
    
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "The bandwidth was determined from uncertainty " << uncertainty << " degrees as: " << this->maxBandwidth;
    ProSHADE_internal_messages::printProgressMessage  ( this->verbose, 3, hlpSS.str() );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function determines the sphere distances for sphere mapping.
 
    This function determines the distance between two consecutive spheres in the sphere mappin galgorithm. It checks
    if this values has not been already set and if not, it sets it as the sampling rate (distance between any two map
    points). It then checks that there will be at least 10 spheres and if not, it changes the sphere distance until at
    least 10 spheres are to be produced.
 
 \param[in] maxMapRange The maximum diagonal distance of the map in Angstroms.
 */
void ProSHADE_settings::determineSphereDistances ( proshade_single maxMapRange )
{
    //================================================ Check the current settings value is set to auto
    if ( this->maxSphereDists != 0.0 )
    {
        std::stringstream hlpSS;
        hlpSS << "The sphere distances were determined as " << this->maxSphereDists << " Angstroms.";
        ProSHADE_internal_messages::printProgressMessage ( this->verbose, 3, hlpSS.str() );
        return ;
    }
    
    //================================================ Determine automatically
    this->maxSphereDists                              = ProSHADE_internal_spheres::autoDetermineSphereDistances ( maxMapRange, this->requestedResolution );
    
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "The sphere distances were determined as " << this->maxSphereDists << " Angstroms.";
    ProSHADE_internal_messages::printProgressMessage  ( this->verbose, 3, hlpSS.str() );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function determines the integration order for the between spheres integration.
 
    This function determines the order of the Gauss-Legendre integration which needs to be done between the spheres. To do
    this, it uses the pre-coputed values of maxium distance between integration points for each order and the maxium distance
    between spheres expressed as a fraction of the total.
 
 \param[in] maxMapRange The maximum diagonal distance of the map in Angstroms.
 */
void ProSHADE_settings::determineIntegrationOrder ( proshade_single maxMapRange )
{
    //================================================ Check the current settings value is set to auto
    if ( this->integOrder != 0 )
    {
        std::stringstream hlpSS;
        hlpSS << "The integration order was determined as " << this->integOrder;
        ProSHADE_internal_messages::printProgressMessage ( this->verbose, 3, hlpSS.str() );
        return ;
    }
    
    //================================================ Determine automatically
    this->integOrder                                  = ProSHADE_internal_spheres::autoDetermineIntegrationOrder ( maxMapRange, this->maxSphereDists );
    
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "The integration order was determined as " << this->integOrder;
    ProSHADE_internal_messages::printProgressMessage  ( this->verbose, 3, hlpSS.str() );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function determines all the required values for spherical harmonics computation.
 
    This function takes the maximum dimension size (in indices) and uses the settings pre-set by the user to set up the
    sphherical harmonics bandwidth, sphere sampling, sphere placement and spacing as well as the Gauss-Legendre integration
    order. This is either done using the user set values (if given), or using automated algorithm which only requires the
    resolution and max dimension.
 
    Note that this function will use the resolution value to modify the values to be appropriate for the resolution supplied and not
    necessarily for the map sampling given.
 
    \param[in] xDim The size of the x axis dimension in indices.
    \param[in] yDim The size of the y axis dimension in indices.
    \param[in] zDim The size of the z axis dimension in indices.
    \param[in] xDimAngs The size of the x-axis in Angstroms.
    \param[in] yDimAngs The size of the y-axis in Angstroms.
    \param[in] zDimAngs The size of the z-axis in Angstroms.
 
    \warning Because the automated algorithm decides the values based on the first structure size, by using it one gives up on
    the idea that DIST(A,B) == DIST(B,A). If this is important, then the user should set all of these values manually to the
    settings object to avoid this issue.
 */
void ProSHADE_settings::determineAllSHValues ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_double xDimAngs, proshade_double yDimAngs, proshade_double zDimAngs )
{
    //================================================ Print progress message
    ProSHADE_internal_messages::printProgressMessage  ( this->verbose, 1, "Preparing spherical harmonics environment." );
    
    //================================================ Modify dims by resolution
    proshade_unsign theoXDim                          = std::ceil ( xDimAngs / ( this->requestedResolution / 2.0 ) );
    proshade_unsign theoYDim                          = std::ceil ( yDimAngs / ( this->requestedResolution / 2.0 ) );
    proshade_unsign theoZDim                          = std::ceil ( zDimAngs / ( this->requestedResolution / 2.0 ) );
    
    //================================================ Find maximum circumference
    proshade_unsign maxDim                            = std::max ( theoXDim, std::max ( theoYDim, theoZDim ) );
    proshade_unsign minDim                            = std::min ( theoXDim, std::min ( theoYDim, theoZDim ) );
    proshade_unsign midDim                            = 0;
    if      ( ( xDim < maxDim ) && ( xDim > minDim ) ) { midDim = theoXDim; }
    else if ( ( yDim < maxDim ) && ( yDim > minDim ) ) { midDim = theoYDim; }
    else                                               { midDim = theoZDim; }
    
    proshade_unsign circ                              = ( maxDim ) + ( midDim );
    
    //================================================ Bandwidth
    if ( this->rotationUncertainty > 0.0 ) { this->determineBandwidthFromAngle ( this->rotationUncertainty ); }
    else { this->determineBandwidth ( circ ); }
    
    //================================================ Find maximum diagonal in Angstroms
    proshade_single maxDiag                           = std::sqrt ( std::pow ( static_cast<proshade_single> ( maxDim ) * ( this->requestedResolution / 2.0 ), 2.0 ) +
                                                                    std::pow ( static_cast<proshade_single> ( midDim ) * ( this->requestedResolution / 2.0 ), 2.0 ) );
    
    //================================================ Sphere distances
    this->determineSphereDistances                    ( maxDiag );
    
    //================================================ Integration order
    this->determineIntegrationOrder                   ( maxDiag );
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( this->verbose, 2, "Spherical harmonics environment prepared." );
    
    //================================================ Done
    return ;
    
}

/*! \brief Contructor for the ProSHADE_run class.
 
    This is where all the decisions regarding what should be done are made. It takes the settings and based on them, it decides what to do
    and how to report the results.
 
    \param[in] settings ProSHADE_settings object specifying what should be done.
 */
ProSHADE_run::ProSHADE_run ( ProSHADE_settings* settings )
{
    //================================================ Wellcome message if required
    ProSHADE_internal_messages::printWellcomeMessage  ( settings->verbose );
    
    //================================================ Save the general information
    this->noStructures                                = static_cast<proshade_unsign> ( settings->inputFiles.size() );
    this->verbose                                     = static_cast<proshade_signed> ( settings->verbose );
    
    //================================================ Try to run ProSHADE
    try
    {
        //============================================ Depending on task, switch to correct function to call
        switch ( settings->task )
        {
            case NA:
                throw ProSHADE_exception ( "No task has been specified.", "E000001", __FILE__, __LINE__, __func__, "ProSHADE requires to be told which particular functiona-\n                    : lity (task) is requested from it. In order to do so, the\n                    : command line arguments specifying task need to be used\n                    : (if used from command line), or the ProSHADE_settings\n                    : object needs to have the member variable \'Task\' set to\n                    : one of the following values: Distances, Symmetry,\n                    : OverlayMap or MapManip." );
                break;
                
            case Symmetry:
                ProSHADE_internal_tasks::SymmetryDetectionTask ( settings, &this->RecomSymAxes, &this->allCSymAxes );
                this->setSymmetryResults              ( settings );
                break;
                
            case Distances:
                ProSHADE_internal_tasks::DistancesComputationTask ( settings, &this->enLevs, &this->trSigm, &this->rotFun );
                break;
                
            case OverlayMap:
                ProSHADE_internal_tasks::MapOverlayTask ( settings, &this->coordRotationCentre, &this->eulerAngles, &this->overlayTranslation );
                break;
                
            case MapManip:
                ProSHADE_internal_tasks::MapManipulationTask ( settings, &this->originalBounds, &this->reboxedBounds, &this->manipulatedMaps );
                break;
        }
    }
    
    //================================================ If this is ProSHADE exception, give all available info and terminate gracefully :-)
    catch ( ProSHADE_exception& err )
    {
        std::cerr << std::endl << "=====================" << std::endl << "!! ProSHADE ERROR !!" << std::endl << "=====================" << std::endl << std::flush;
        std::cerr << "Error Code          : " << err.get_errc() << std::endl << std::flush;
        std::cerr << "ProSHADE version    : " << __PROSHADE_VERSION__ << std::endl << std::flush;
        std::cerr << "File                : " << err.get_file() << std::endl << std::flush;
        std::cerr << "Line                : " << err.get_line() << std::endl << std::flush;
        std::cerr << "Function            : " << err.get_func() << std::endl << std::flush;
        std::cerr << "Message             : " << err.what() << std::endl << std::flush;
        std::cerr << "Further information : " << err.get_info() << std::endl << std::endl << std::flush;
        
        //============================================ Done
        ProSHADE_internal_messages::printTerminateMessage ( settings->verbose );
        exit                                          ( EXIT_FAILURE );
    }
    
    //================================================ Well, give all there is and just end
    catch ( ... )
    {
        std::cerr << std::endl << "=====================" << std::endl << "!! ProSHADE ERROR !!" << std::endl << "=====================" << std::endl << std::flush;
        
        //============================================ Try to find out more
#if __cplusplus >= 201103L
            std::exception_ptr exc                    = std::current_exception();
            try
            {
                if (exc)
                {
                    std::rethrow_exception            ( exc );
                }
            }
            catch ( const std::exception& e )
            {
                std::cerr << "Caught unknown exception with following information: " << e.what() << std::endl << std::flush;
            }
#else
            std::cerr << "Unknown error with no further explanation available. Please contact the author for help." << std::endl << std::flush;
#endif
        std::cerr << "Terminating..." << std::endl << std::endl << std::flush;
        
        //============================================ Done
        ProSHADE_internal_messages::printTerminateMessage ( settings->verbose );
        exit                                          ( EXIT_FAILURE );
    }
    
    //================================================ Terminating message
    ProSHADE_internal_messages::printTerminateMessage ( settings->verbose );
    
    //================================================ Done
    
}

/*! \brief Destructor for the ProSHADE class.
 
    This destructor is responsible for releasing all memory used by the executing object
 */
ProSHADE_run::~ProSHADE_run ( )
{
    //================================================ Release reboxing pointers
    if ( this->originalBounds.size() > 0 ) { for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( this->originalBounds.size() ); iter++ ) { delete[] this->originalBounds.at(iter); } }
    if ( this->reboxedBounds.size()  > 0 ) { for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( this->reboxedBounds.size()  ); iter++ ) { delete[] this->reboxedBounds.at(iter); } }
    if ( this->manipulatedMaps.size() > 0 ) { for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( this->manipulatedMaps.size() ); iter++ ) { delete[] this->manipulatedMaps.at(iter); } }
    
    //================================================ Clear vectors
    this->enLevs.clear                                ( );
    this->trSigm.clear                                ( );
    this->rotFun.clear                                ( );
    
    //================================================ Delete symmetry axes memory
    if ( this->RecomSymAxes.size() > 0 )
    {
        for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( this->RecomSymAxes.size() ); iter++ )
        {
            delete[] this->RecomSymAxes.at(iter);
        }
        this->RecomSymAxes.clear                      ( );
    }
    
    //================================================ Done
    
}

/*! \brief This is the main accessor function for the user to get to know what symmetry type ProSHADE has detected and recommends.
 
    \param[out] symRecommType This is the value ( ""=None, C=cyclic, D=Dihedral, T=Tetrahedral, O=Octahedral or I=Icosahedral) of ProSHADE detected and recommended symmetry.
 */
std::string ProSHADE_run::getSymmetryType ( )
{
    //================================================ Return the value
    return                                            ( this->symRecommType );
}

/*! \brief This is the main accessor function for the user to get to know what symmetry fold ProSHADE has detected and recommends.
 
    \param[out] symRecommFold This is the fold of ProSHADE detected and recommended symmetry (C and D symmetry types only).
 */
proshade_unsign ProSHADE_run::getSymmetryFold ( )
{
    //================================================ Return the value
    return                                            ( this->symRecommFold );
}

/*! \brief Sets the ProSHADE detected symmetry type.
 
    When symmetry detection is done, the resulting recommended symmetry type will be saved in the ProSHADE object by this function.
 
    \param[in] val The recommended symmetry type for the structure.
 */
void ProSHADE_run::setRecommendedSymmetry ( std::string val )
{
    //================================================ Set the value
    this->symRecommType                               = val;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the ProSHADE detected symmetry fold.
 
    When symmetry detection is done, the resulting recommended symmetry fold  (valid only for C and D symmetry types) will be saved in the
    ProSHADE object by this function.
 
    \param[in] val The recommended symmetry fold for the structure.
 */
void ProSHADE_run::setRecommendedFold ( proshade_unsign val )
{
    //================================================ Set the value
    this->symRecommFold                               = val;
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the ProSHADE detected symmetry axis.
 
    When symmetry detection is done, the resulting recommended symmetry axis will be saved in the
    ProSHADE object by this function.
 
    \param[in] sym The recommended symmetry axis for the structure.
 */
void ProSHADE_run::setRecommendedAxis ( proshade_double* sym )
{
    //================================================ Set the value
    ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( &this->RecomSymAxes, sym );
    
    //================================================ Done
    return ;
    
}

/*! \brief Sets the ProSHADE detected symmetry information for easy programmatical output.
 
    When symmetry detection is done, the resulting recommended symmetry information will be saved in the
    ProSHADE object by this function.
 
    \param[in] settings ProSHADE_settings object where the results are passed through.
 */
void ProSHADE_run::setSymmetryResults ( ProSHADE_settings* settings )
{
    //================================================ Save type and fold
    this->setRecommendedSymmetry                      ( settings->recommendedSymmetryType );
    this->setRecommendedFold                          ( settings->recommendedSymmetryFold );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function parses the command line arguments into the settings object.

    \param[in] argc The count of the command line arguments (as passed to main function by the system).
    \param[in] argv The string containing the command line arguments (as passed to main function by the system).
 */
void ProSHADE_settings::getCommandLineParams ( int argc, char** argv )
{
    //================================================ If no command line arguments, print help
    if ( argc == 1 ) { ProSHADE_internal_messages::printHelp ( ); }
    
    //================================================ Long options struct
    const struct option_port longopts[] =
    {
        { "version",         no_argument,        NULL, 'v' },
        { "help",            no_argument,        NULL, 'h' },
        { "verbose",         required_argument,  NULL, '!' },
        { "distances",       no_argument,        NULL, 'D' },
        { "mapManip",        no_argument,        NULL, 'M' },
        { "symmetry",        no_argument,        NULL, 'S' },
        { "overlay",         no_argument,        NULL, 'O' },
        { "file",            required_argument,  NULL, 'f' },
        { "forceSpgP1",      no_argument,        NULL, 'u' },
        { "removeWaters",    no_argument,        NULL, 'w' },
        { "firstModel",      no_argument,        NULL, 'x' },
        { "resolution",      required_argument,  NULL, 'r' },
        { "bandwidth",       required_argument,  NULL, 'b' },
        { "sphereDists",     required_argument,  NULL, 's' },
        { "extraSpace",      required_argument,  NULL, 'e' },
        { "integOrder",      required_argument,  NULL, 'i' },
        { "taylorCap",       required_argument,  NULL, 't' },
        { "invertMap",       no_argument,        NULL, '@' },
        { "normalise",       no_argument,        NULL, '#' },
        { "mask",            no_argument,        NULL, '$' },
        { "saveMask",        no_argument,        NULL, '%' },
        { "maskFile",        required_argument,  NULL, '^' },
        { "maskBlurring",    required_argument,  NULL, '&' },
        { "maskThreshold",   required_argument,  NULL, '*' },
        { "mapReboxing",     no_argument,        NULL, 'R' },
        { "boundsSpace",     required_argument,  NULL, '(' },
        { "boundsThreshold", required_argument,  NULL, ')' },
        { "sameBoundaries",  no_argument,        NULL, '-' },
        { "reBoxedFilename", required_argument,  NULL, 'g' },
        { "pdbTempFact",     required_argument,  NULL, 'd' },
        { "center",          no_argument,        NULL, 'c' },
        { "changeMapResol",  no_argument,        NULL, 'j' },
        { "changeMapTriLin", no_argument,        NULL, 'a' },
        { "noPhase",         no_argument,        NULL, 'p' },
        { "progressive",     no_argument,        NULL, 'k' },
        { "noEnL",           no_argument,        NULL, 'l' },
        { "noTrS",           no_argument,        NULL, 'm' },
        { "noFRF",           no_argument,        NULL, 'n' },
        { "EnLWeight",       required_argument,  NULL, '_' },
        { "peakNeigh",       required_argument,  NULL, '=' },
        { "peakThres",       required_argument,  NULL, '+' },
        { "missAxThres",     required_argument,  NULL, '[' },
        { "sameAxComp",      required_argument,  NULL, ']' },
        { "axisComBeh",      no_argument,        NULL, 'q' },
        { "bicubSearch",     no_argument,        NULL, 'A' },
        { "maxSymPrime",     required_argument,  NULL, 'B' },
        { "minPeakHeight",   required_argument,  NULL, 'o' },
        { "reqSym",          required_argument,  NULL, '{' },
        { "overlayFile",     required_argument,  NULL, '}' },
        { "overlayJSONFile", required_argument,  NULL, 'y' },
        { "angUncertain",    required_argument,  NULL, ';' },
        { "usePeaksInRotFun",no_argument,        NULL, 'z' },
        { NULL,              0,                  NULL,  0  }
    };
    
    //================================================ Short options string
    const char* const shortopts                       = "AaB:b:cd:De:f:g:hi:jklmMno:Opqr:Rs:St:uvwxy:z!:@#$%^:&:*:(:):-_:=:+:[:]:{:}:;:";
    
    //================================================ Parsing the options
    while ( true )
    {
        //============================================ Read the next option
        int opt                                       = getopt_long_port ( argc, argv, shortopts, longopts, NULL );
        
        //============================================ Done parsing
        if ( opt == -1 )
        {
            break;
        }
        
        //============================================ For each option, set the internal values appropriately
         switch ( opt )
         {
             //======================================= Print version info
             case 'v':
             {
                 ProSHADE_internal_messages::printWellcomeMessage ( 0 );
                 exit                                 ( EXIT_SUCCESS );
             }
                 
             //======================================= User needs help
             case 'h':
             {
                 ProSHADE_internal_messages::printHelp ( );
                 exit                                 ( EXIT_SUCCESS );
             }
                 
             //======================================= Save the argument as the verbosity value, or if no value given, just set to 3
             case '!':
             {
                 this->setVerbosity                   ( static_cast<proshade_single> ( atoi ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Set task to distances
             case 'D':
             {
                 this->task                           = Distances;
                 continue;
             }
                 
             //======================================= Set task to map manipulation
             case 'M':
             {
                 this->task                           = MapManip;
                 continue;
             }
                 
             //======================================= Set task to symmetry detection
             case 'S':
             {
                 this->task                           = Symmetry;
                 
                 //=================================== Force default unless changed already by the user
                 if (  this->requestedResolution == -1 ) { this->requestedResolution = 6.0;  }
                 if (  this->pdbBFactorNewVal    == -1 ) { this->pdbBFactorNewVal    = 80.0; }
                 this->changeMapResolution            = !this->changeMapResolution;  // Switch value. This can be over-ridden by the user by using -j
                 this->moveToCOM                      = !this->moveToCOM;            // Switch value. This can be over-ridden by the user by using -c.
                 
                 continue;
             }
                 
             //======================================= Set task to map overlay
             case 'O':
             {
                 this->task                           = OverlayMap;
                 continue;
             }
                 
             //======================================= Save the argument as a file to read in
             case 'f':
             {
                 this->addStructure                   ( std::string ( optarg ) );
                 continue;
             }
                 
             //======================================= Force the input PDB files to have P1 spacegroup
             case 'u':
             {
                 this->forceP1                        = !this->forceP1;
                 continue;
             }
                 
             //======================================= Remove waters from PDB input files?
             case 'w':
             {
                 this->removeWaters                   = !this->removeWaters;
                 continue;
             }
                 
             //======================================= Use all models, or just the first one?
             case 'x':
             {
                 this->firstModelOnly                 = !this->firstModelOnly;
                 continue;
             }
                 
             //======================================= Save the argument as the resolution value
             case 'r':
             {
                 this->setResolution                  ( static_cast<proshade_single> ( atof ( optarg ) ) );
                 continue;
             }

             //======================================= Save the argument as the bandwidth value
             case 'b':
             {
                 this->setBandwidth                   ( static_cast<proshade_unsign> ( atoi ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Save the argument as the extra space value
             case 'e':
             {
                 this->setExtraSpace                  ( static_cast<proshade_single> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Save the argument as the intaggration order value
             case 'i':
             {
                 this->setIntegrationOrder            ( static_cast<proshade_unsign> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Save the argument as the sphere distance value
             case 's':
             {
                 this->setSphereDistances             ( static_cast<proshade_single> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Save the argument as the taylor series cap value
             case 't':
             {
                 this->setTaylorSeriesCap             ( static_cast<proshade_unsign> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Set map inversion to true
             case '@':
             {
                 this->setMapInversion                ( true );
                 continue;
             }
                 
             //======================================= Set map normalisation to true
             case '#':
             {
                 this->setNormalisation               ( true );
                 continue;
             }
                 
             //======================================= Set map masking to true
             case '$':
             {
                 this->setMasking                     ( true );
                 continue;
             }
                 
             //======================================= Set map masking to true and mask map saving to true as well
             case '%':
             {
                 this->setMasking                     ( true );
                 this->setMaskSaving                  ( true );
                 continue;
             }
                 
             //======================================= Save the argument as the mask filename value
             case '^':
             {
                 this->setMaskFilename                ( static_cast<std::string> ( optarg ) );
                 continue;
             }
                 
             //======================================= Save the argument as the mask blurring factor value
             case '&':
             {
                 this->setMaskBlurFactor              ( static_cast<proshade_single> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Save the argument as the mask threshold (IQR) value
             case '*':
             {
                 this->setMaskIQR                     ( static_cast<proshade_single> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Set map reboxing to true
             case 'R':
             {
                 this->setMasking                     ( true );
                 this->setMapReboxing                 ( true );
                 continue;
             }
                 
             //======================================= Save the argument as the bounds extra space value
             case '(':
             {
                 this->setBoundsSpace                 ( static_cast<proshade_single> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Save the argument as the bounds similarity threshold value
             case ')':
             {
                 this->setBoundsThreshold             ( static_cast<proshade_signed> ( atoi ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Set same boundaries to true
             case '-':
             {
                 this->setSameBoundaries              ( true );
                 continue;
             }
                 
             //======================================= Save the argument as the re-boxed structure filename value
             case 'g':
             {
                 this->setOutputFilename              ( static_cast<std::string> ( optarg ) );
                 continue;
             }
                 
             //======================================= Save the argument as the PDB B-factor new constant value
             case 'd':
             {
                 this->setPDBBFactor                  ( static_cast<proshade_single> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Set map centering to true
             case 'c':
             {
                 this->moveToCOM                      = !this->moveToCOM;
                 continue;
             }
                 
             //======================================= Set map resolution change using Fourier transforms to true
             case 'j':
             {
                 this->changeMapResolution            = !this->changeMapResolution;
                 continue;
             }
                 
             //======================================= Set map resolution change using real-space tri-linear interpolation to true
             case 'a':
             {
                 this->setMapResolutionChangeTriLinear ( true );
                 continue;
             }
                 
             //======================================= Set map phase removal to true
             case 'p':
             {
                 this->setPhaseUsage                  ( false );
                 continue;
             }
                 
             //======================================= Set progressive shell mapping to true
             case 'k':
             {
                 this->setProgressiveSphereMapping    ( true );
                 continue;
             }
                 
             //======================================= Set energy level descriptor computation to false
             case 'l':
             {
                 this->setEnergyLevelsComputation     ( false );
                 continue;
             }
                 
             //======================================= Set trace sigma descriptor computation to false
             case 'm':
             {
                 this->setTraceSigmaComputation       ( false );
                 continue;
             }
                 
             //======================================= Set full rotation function descriptor computation to false
             case 'n':
             {
                 this->setRotationFunctionComputation ( false );
                 continue;
             }
                 
             //======================================= Save the argument as the energy levels descriptor weight value
             case '_':
             {
                 this->setEnLevShellWeight            ( static_cast<proshade_double> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Save the argument as the peak neighbours minimum value
             case '=':
             {
                 this->setPeakNeighboursNumber        ( static_cast<proshade_unsign> ( atoi ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Save the argument as the peak IQR from median naive small peaks removal value
             case '+':
             {
                 this->setPeakNaiveNoIQR              ( static_cast<proshade_double> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Save the argument as the missing axis threshold value
             case '[':
             {
                 this->setMissingPeakThreshold        ( static_cast<proshade_double> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Save the argument as the missing axis threshold value
             case ']':
             {
                 setAxisComparisonThreshold           ( static_cast<proshade_double> ( atof ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Save the argument as the missing axis threshold value
             case 'q':
             {
                 setAxisComparisonThresholdBehaviour  ( !this->axisErrToleranceDefault );
                 continue;
             }
                 
             //======================================= Save the argument as the bicubic interpolation search requirement value
             case 'A':
             {
                 setBicubicInterpolationSearch        ( !this->useBiCubicInterpolationOnPeaks );
                 continue;
             }
                 
             //======================================= Save the argument as the bicubic interpolation search requirement value
             case 'B':
             {
                 setMaxSymmetryFold                   ( static_cast<proshade_unsign> ( atoi ( optarg ) ) );
                 continue;
             }
                 
             //======================================= Minimum peak height for axis
             case 'o':
             {
                 this->minSymPeak                     = static_cast<proshade_double> ( atof ( optarg ) );
                 continue;
             }
                 
             //======================================= Save the argument as the requested symmetry and potentially fold value
             case '{':
             {
                 std::string input                    = static_cast<std::string> ( optarg );
                 
                 if ( input.at(0) == 'C' )
                 {
                     this->setRequestedSymmetry ( "C" );
                     
                     std::string numHlp ( input.begin()+1, input.end() );
                     if ( numHlp.length() > 0 ) { this->setRequestedFold ( atoi ( numHlp.c_str() ) ); }
                     else { std::cerr << "!!! ProSHADE ERROR !!! The input argument requests search for Cyclic/Dihedral symmetry, but does not specify the requested fold." << std::endl;  exit ( EXIT_FAILURE ); }
                 }
                 else
                 {
                     if ( input.at(0) == 'D' )
                     {
                         this->setRequestedSymmetry ( "D" );
                         
                         std::string numHlp ( input.begin()+1, input.end() );
                         if ( numHlp.length() > 0 ) { this->setRequestedFold ( atoi ( numHlp.c_str() ) ); }
                         else { std::cerr << "!!! ProSHADE ERROR !!! The input argument requests search for Cyclic/Dihedral symmetry, but does not specify the requested fold." << std::endl;  exit ( EXIT_FAILURE ); }
                     }
                     else
                     {
                         if ( input.at(0) == 'T' )
                         {
                             this->setRequestedSymmetry ( "T" );
                         }
                         else
                         {
                             if ( input.at(0) == 'O' )
                             {
                                 this->setRequestedSymmetry ( "O" );
                             }
                             else
                             {
                                 if ( input.at(0) == 'I' )
                                 {
                                     this->setRequestedSymmetry ( "I" );
                                 }
                                 else
                                 {
                                     std::cerr << "!!! ProSHADE ERROR !!! Failed to parse the requested symmetry type. Allowed types are C, D, T, O and I, with C and D requiring to be followed by a number specifying the fold." << std::endl;  exit ( EXIT_FAILURE );
                                 }
                             }
                         }
                     }
                 }
                     
                 continue;
             }
                 
             //======================================= Save the argument as filename to save the overlay moved structure to value
             case '}':
             {
                 this->setOverlaySaveFile             ( static_cast<std::string> ( optarg ) );
                 continue;
             }
                 
             //======================================= Save the argument as filename to save the overlay operations to value
             case 'y':
             {
                 this->setOverlayJsonFile             ( static_cast<std::string> ( optarg ) );
                 continue;
             }
                 
             //======================================= Save the argument as angular uncertainty for bandwidth determination
             case ';':
             {
                 this->rotationUncertainty            = static_cast<proshade_double> ( atof ( optarg ) );
                 continue;
             }
                 
             //======================================= Save the argument as angular uncertainty for bandwidth determination
             case 'z':
             {
                 this->setSymmetryRotFunPeaks         ( false );
                 continue;
             }
                 
             //======================================= Unknown option
             case '?':
             {
                 //=================================== Save the argument as angular uncertainty for bandwidth determination
                 if ( optopt )
                 {
                     std::cout << "!!! ProSHADE ERROR !!! Unrecognised short option -" << static_cast<char> ( optopt ) << " . Please use -h for help on the command line options." << std::endl;
                 }
                 else
                 {
                     std::cout << "!!! ProSHADE ERROR !!! Unrecognised long option " << argv[static_cast<int> (optind)-1] << " . Please use -h for help on the command line options." << std::endl;
                 }
                 
                 //=================================== This case is handled by getopt_long, nothing more needed.
                 exit ( EXIT_SUCCESS );
             }
                 
             //======================================= Fallback option
             default:
             {
                 ProSHADE_internal_messages::printHelp ( );
                 exit ( EXIT_SUCCESS );
             }
         }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function prints the current values in the settings object.
 
    \warning This is a debugging function of no real utility to the user.
 */
void ProSHADE_settings::printSettings ( )
{
    //================================================ Print the currest values in the settings object
    std::stringstream strstr;
    strstr.str(std::string());
    if ( this->task == NA ) { strstr << "NA"; }
    if ( this->task == Distances ) { strstr << "DISTANCES COMPUTATION"; }
    if ( this->task == MapManip ) { strstr << "MAP MANIPULATION"; }
    if ( this->task == Symmetry ) { strstr << "SYMMETRY DETECTION"; }
    if ( this->task == OverlayMap ) { strstr << "MAP OVERLAY"; }
    printf ( "Task to perform     : %37s\n", strstr.str().c_str() );
    
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( this->inputFiles.size() ); iter++ )
    {
        strstr.str(std::string());
        strstr << this->inputFiles.at(iter);
        printf ( "File(s) to process  : %37s\n", strstr.str().c_str() );
    }
    
    strstr.str(std::string());
    strstr << this->verbose;
    printf ( "Verbosity           : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->requestedResolution;
    printf ( "Resolution (comp)   : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->maxBandwidth;
    printf ( "Bandwidth           : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->maxSphereDists;
    printf ( "Sphere distances    : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->addExtraSpace;
    printf ( "Extra space         : %37s\n", strstr.str().c_str() );

    strstr.str(std::string());
    if ( this->forceP1 ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Force P1 spacegroup : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->removeWaters ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Waters removed      : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->firstModelOnly ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Only 1st model      : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->integOrder;
    printf ( "Integration order   : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->taylorSeriesCap;
    printf ( "Taylor series cap   : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->pdbBFactorNewVal;
    printf ( "PDB B-factor const  : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->reBoxMap ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Map re-boxing       : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->invertMap ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Map inversion       : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->normaliseMap ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Map normalisation   : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->maskMap ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Map masking         : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->saveMask ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Saving mask         : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->maskFileName;
    printf ( "Map mask filename   : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->blurFactor;
    printf ( "Map blurring        : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->maskingThresholdIQRs;
    printf ( "Masking threshold   : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->boundsExtraSpace;
    printf ( "Bounds extra space  : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->boundsSimilarityThreshold;
    printf ( "Bounds similarity   : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->useSameBounds ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Same boundaries     : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->outName;
    printf ( "Re-boxed filename   : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->moveToCOM ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Map COM centering   : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->changeMapResolution ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Change map resol    : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->changeMapResolutionTriLinear ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Change map tri-lin  : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->usePhase ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Use phase info      : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->progressiveSphereMapping ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Progressive spheres : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->computeEnergyLevelsDesc ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Energy lvl desc     : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->enLevMatrixPowerWeight;
    printf ( "Energy lvl weight   : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->computeTraceSigmaDesc ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Tr sigma desc       : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->computeRotationFuncDesc ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Full RF desc        : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->peakNeighbours;
    printf ( "Neightbours to peak : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->noIQRsFromMedianNaivePeak;
    printf ( "Peak IQR threshold  : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->symMissPeakThres;
    printf ( "Missing ax. thres   : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->minSymPeak;
    printf ( "Min. sym. peak size : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->axisErrTolerance;
    printf ( "Same ax. threshold  : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    if ( this->axisErrToleranceDefault ) { strstr << "TRUE"; } else { strstr << "FALSE"; }
    printf ( "Same ax. thre. decr.: %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->requestedSymmetryType << "-" << this->requestedSymmetryFold;
    printf ( "Requested symm.     : %37s\n", strstr.str().c_str() );

    strstr.str(std::string());
    strstr << this->overlayStructureName;
    printf ( "Overlay file        : %37s\n", strstr.str().c_str() );
    
    strstr.str(std::string());
    strstr << this->rotTrsJSONFile;
    printf ( "JSON overlay file   : %37s\n", strstr.str().c_str() );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function returns the energy level distances vector from the first to all other structures.
 
    \param[out] enLevs Vector of doubles of the distances.
 */
std::vector< proshade_double > ProSHADE_run::getEnergyLevelsVector ( )
{
    //================================================ Return the value
    return                                            ( this->enLevs );
}

/*! \brief This function returns the trace sigma distances vector from the first to all other structures.
 
    \param[out] trSigm Vector of doubles of the distances.
 */
std::vector< proshade_double > ProSHADE_run::getTraceSigmaVector ( )
{
    //================================================ Return the value
    return                                            ( this->trSigm );
}

/*! \brief This function returns the full rotation function distances vector from the first to all other structures.
 
    \param[out] rotFun Vector of doubles of the distances.
 */
std::vector< proshade_double > ProSHADE_run::getRotationFunctionVector ( )
{
    //================================================ Return the value
    return                                            ( this->rotFun );
}

/*! \brief This function returns the number of structures used.

    \param[in] noStructures Number of structures supplied to the settings object.
*/
proshade_unsign ProSHADE_run::getNoStructures ( )
{
    //================================================ Return the value
    return                                            ( this->noStructures );
}

/*! \brief This function returns the verbose value.

    \param[in] verbose How loud the run should be?
*/
proshade_signed ProSHADE_run::getVerbose ( )
{
    //================================================ Return the value
    return                                            ( this->verbose );
}

/*! \brief This function returns the number of detected recommended symmetry axes.

    \param[out] val The length of the recommended symmetry axes vector.
*/
proshade_unsign ProSHADE_run::getNoSymmetryAxes ( )
{
    //================================================ Return the value
    return                                            ( static_cast<proshade_unsign> ( this->RecomSymAxes.size() ) );
}

/*! \brief This function returns the number of detected recommended symmetry axes.

    \param[out] val The length of the recommended symmetry axes vector.
*/
proshade_unsign ProSHADE_run::getNoRecommendedSymmetryAxes ( )
{
    //================================================ Return the value
    return                                            ( static_cast<proshade_unsign> ( this->RecomSymAxes.size() ) );
}

/*! \brief This function returns a single symmetry axis as a vector of strings from the recommended symmetry axes list.

    \param[in] axisNo The index of the axis to be returned.
    \param[out] val A vector of strings containing the symmetry axis fold, x, y, z axis element, angle and peak height in this order.
*/
std::vector< std::string > ProSHADE_run::getSymmetryAxis ( proshade_unsign axisNo )
{
    //================================================ Sanity checks
    if ( static_cast<proshade_unsign> ( this->RecomSymAxes.size() ) <= axisNo )
    {
        ProSHADE_internal_messages::printWarningMessage ( this->verbose, "!!! ProSHADE WARNING !!! Requested symmetry index does not exist. Returning empty vector.", "WS00039" );
        return                                        ( std::vector< std::string > ( ) );
    }
    
    //================================================ Initialise local variables
    std::vector< std::string > ret;
    
    //================================================ Input the axis data as strings
    std::stringstream ssHlp;
    ssHlp << this->RecomSymAxes.at(axisNo)[0];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
    
    ssHlp << this->RecomSymAxes.at(axisNo)[1];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
            
    ssHlp << this->RecomSymAxes.at(axisNo)[2];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
            
    ssHlp << this->RecomSymAxes.at(axisNo)[3];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
            
    ssHlp << this->RecomSymAxes.at(axisNo)[4];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
            
    ssHlp << this->RecomSymAxes.at(axisNo)[5];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function returns a all symmetry axes as a vector of vectors of doubles.

    \param[out] val A vector of vectors of doubles containing all the symmetries axis fold, x, y, z axis element, angle and peak height in this order.
*/
std::vector < std::vector< proshade_double > > ProSHADE_run::getAllCSyms ( )
{
    //================================================ Done
    return                                            ( this->allCSymAxes );
    
}

/*! \brief This function returns a specific structure original bounds.

    \param[in] strNo The index of the structure for which the bounds are to be returned.
*/
std::vector< proshade_signed > ProSHADE_run::getOriginalBounds ( proshade_unsign strNo )
{
    //================================================ Sanity checks
    if ( noStructures <= strNo )
    {
        ProSHADE_internal_messages::printWarningMessage ( this->verbose, "!!! ProSHADE WARNING !!! Requested bounds for structure index which does not exist. Returning empty vector.", "WB00041" );
        return                                        ( std::vector< proshade_signed > ( ) );
    }
    
    //================================================ Initialise local variables
    std::vector< proshade_signed > ret;
    
    //================================================ Input the axis data as strings
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->originalBounds.at( strNo )[0] );
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->originalBounds.at( strNo )[1] );
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->originalBounds.at( strNo )[2] );
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->originalBounds.at( strNo )[3] );
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->originalBounds.at( strNo )[4] );
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->originalBounds.at( strNo )[5] );
    
    //================================================ Done
    return                                            ( ret );
}

/*! \brief This function returns a specific structure re-boxed bounds.

    \param[in] strNo The index of the structure for which the bounds are to be returned.
*/
std::vector< proshade_signed > ProSHADE_run::getReBoxedBounds ( proshade_unsign strNo )
{
    //================================================ Sanity checks
    if ( noStructures <= strNo )
    {
        ProSHADE_internal_messages::printWarningMessage ( this->verbose, "!!! ProSHADE WARNING !!! Requested bounds for structure index which does not exist. Returning empty vector.", "WB00041" );
        return                                        ( std::vector< proshade_signed > ( ) );
    }
    
    //================================================ Initialise local variables
    std::vector< proshade_signed > ret;
    
    //================================================ Input the axis data as strings
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->reboxedBounds.at( strNo )[0] );
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->reboxedBounds.at( strNo )[1] );
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->reboxedBounds.at( strNo )[2] );
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->reboxedBounds.at( strNo )[3] );
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->reboxedBounds.at( strNo )[4] );
    ProSHADE_internal_misc::addToSignedVector         ( &ret, this->reboxedBounds.at( strNo )[5] );
    
    //================================================ Done
    return                                            ( ret );
}

/*! \brief This function returns a single, specific structure map value.

    \param[in] strNo The index of the structure for which the map value is to be returned.
    \param[in] mapIndex The map array index of which the value is returned.
    \param[out] val The map density value for the particular mapIndex position.
*/
proshade_double ProSHADE_run::getMapValue ( proshade_unsign strNo, proshade_unsign mapIndex )
{
    //================================================ Return the value
    return                                            ( this->manipulatedMaps.at(strNo)[mapIndex] );
}

/*! \brief This function returns the re-boxed structure map 1D array for the processed structure.
 
    \param[in] run The ProSHADE_run object from which the values will be drawn.
    \param[in] strNo Index of the structure for which the bounds are to be returned.
    \param[in] reboxMap Array to which the values are to be loaded into.
    \param[in] len The length of the array.
 */

void getReBoxedMap ( ProSHADE_run* run, proshade_unsign strNo, double *reboxMap, int len )
{
    //================================================ Sanity checks
    if ( run->getNoStructures() <= strNo )
    {
        ProSHADE_internal_messages::printWarningMessage ( run->getVerbose(), "!!! ProSHADE WARNING !!! Requested bounds for structure index which does not exist. Returning empty vector.", "WB00041" );
        return ;
    }
    
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++)
    {
        reboxMap[iter]                                = static_cast<double> ( run->getMapValue ( strNo, iter ) );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function returns the vector of Euler angles with best overlay correlation.

    \param[out] ret Vector of Euler angles (ZXZ convention) which lead to the globally best overlay correlation.
*/
std::vector< proshade_double > ProSHADE_run::getEulerAngles ( )
{
    //================================================ Sanity check
    if ( this->eulerAngles.size() != 3 )
    {
        ProSHADE_internal_messages::printWarningMessage ( this->verbose, "!!! ProSHADE WARNING !!! Requested rotation/translation values for Overlay functionality without having successfully computed it. Please check the correct task was used and no other warnings/errors were obtained.", "WO00042" );
        return                                        ( std::vector< proshade_double > ( ) );
    }
    
    //================================================ Return required value
    return                                            ( this->eulerAngles );
    
}

/*! \brief This function returns the vector forming rotation matrix (rows first) with best overlay correlation.

    \param[out] ret Vector forming rotation matrix (rows first) which lead to the globally best overlay correlation.
*/
std::vector< proshade_double > ProSHADE_run::getOptimalRotMat ( )
{
    //================================================ Sanity check
    if ( this->eulerAngles.size() != 3 )
    {
        ProSHADE_internal_messages::printWarningMessage ( this->verbose, "!!! ProSHADE WARNING !!! Requested rotation/translation values for Overlay functionality without having successfully computed it. Please check the correct task was used and no other warnings/errors were obtained.", "WO00042" );
        return                                        ( std::vector< proshade_double > ( ) );
    }
    
    //================================================ Obtain the optimal rotation matrix
    proshade_double* rotMat                           = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( this->eulerAngles.at(0), this->eulerAngles.at(1), this->eulerAngles.at(2), rotMat );
    
    //================================================ Copy to the output variable
    std::vector< proshade_double > ret;
    for ( proshade_unsign iter = 0; iter < 9; iter++ ) { ProSHADE_internal_misc::addToDoubleVector ( &ret, rotMat[iter] ); }
    
    //================================================ Release the memory
    delete[] rotMat;
    
    //================================================ Return required value
    return                                            ( ret );
    
}

/*! \brief This function returns the negative values of the position of the rotation centre (the point about which the rotation should be done).

    \param[out] ret Vector specifying the negative values of the rotation centre - i.e. the translation of the rotation centre to the origin.
*/
std::vector< proshade_double > ProSHADE_run::getTranslationToOrigin ( )
{
    //================================================ Sanity check
    if ( this->coordRotationCentre.size() != 3 )
    {
        ProSHADE_internal_messages::printWarningMessage ( this->verbose, "!!! ProSHADE WARNING !!! Requested rotation/translation values for Overlay functionality without having successfully computed it. Please check the correct task was used and no other warnings/errors were obtained.", "WO00042" );
        return                                        ( std::vector< proshade_double > ( ) );
    }
    
    //================================================ Create return variable with negative values of the internal varariable
    std::vector < proshade_double > ret;
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, -this->coordRotationCentre.at(0) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, -this->coordRotationCentre.at(1) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, -this->coordRotationCentre.at(2) );
    
    //================================================ Return required value
    return                                            ( ret );
    
}

/*! \brief This function returns the translation required to move the structure from origin to optimal overlay.

    \param[out] ret Translation required to move structure from origin to optimal overlay.
*/
std::vector< proshade_double > ProSHADE_run::getOriginToOverlayTranslation ( )
{
    //================================================ Sanity check
    if ( this->overlayTranslation.size() != 3 )
    {
        ProSHADE_internal_messages::printWarningMessage ( this->verbose, "!!! ProSHADE WARNING !!! Requested rotation/translation values for Overlay functionality without having successfully computed it. Please check the correct task was used and no other warnings/errors were obtained.", "WO00042" );
        return                                        ( std::vector< proshade_double > ( ) );
    }
    
    //================================================ Return required value
    return                                            ( this->overlayTranslation );
    
}
