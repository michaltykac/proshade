/*! \file ProSHADE_settings.hpp
    \brief This header file declares the ProSHADE_settings class, which is the main information carrying object in ProSHADE.
 
    The ProSHADE_settings class declared in this header file is the main information carrying object, which is passed to many functions in order to avoid passing plentitude of arguments and having to
    change the function signatures every time a new variable is added. It is also the object which the user needs to fill in in order to start any ProSHADE run.
 
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
#include "ProSHADE_maths.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wall"

//==================================================== Gemmi
#ifndef __PROSHADE_GEMMI_INCLUDE__
    #define __PROSHADE_GEMMI_INCLUDE__
    #include <gemmi/mmread.hpp>
    #include <gemmi/ccp4.hpp>
    #include <gemmi/it92.hpp>
    #include <gemmi/dencalc.hpp>
    #include <gemmi/fprime.hpp>
    #include <gemmi/gz.hpp>
#endif

//==================================================== FFTW3
#ifdef __cplusplus
extern "C" {
#endif
    
#include <fftw3.h>
    
#ifdef __cplusplus
}
#endif

//==================================================== SOFT
#ifdef __cplusplus
extern "C" {
#endif
    
#include <wrap_fftw.h>
#include <makeweights.h>
#include <s2_primitive.h>
#include <s2_cospmls.h>
#include <s2_legendreTransforms.h>
#include <s2_semi_fly.h>
#include <rotate_so3_utils.h>
#include <utils_so3.h>
#include <soft_fftw.h>
#include <rotate_so3_fftw.h>
    
#ifdef __cplusplus
}
#endif

#pragma GCC diagnostic pop

//==================================================== GetOpt port (BSD License, works on Windows as well as linux)
#include <getopt_port/getopt_port.h>

//==================================================== Overinclusion protection
#ifndef __PROSHADE_SETTINGS__
#define __PROSHADE_SETTINGS__

/*! \class ProSHADE_settings
    \brief This class stores all the settings and is passed to the executive classes instead of a multitude of parameters.
 
    The ProSHADE_settings class is a simple way of keeping all the settings together and easy to set by the user. Its
    constructor sets it to the default settings, so that if the user does not want to change these, he just needs to
    pass the object to the executing class and all is done.
 */
class ProSHADE_settings
{
public:
    //================================================ Settings regarding the task at hand
    ProSHADE_Task task;                               //!< This custom type variable determines which task to perfom (i.e. symmetry detection, distances computation, etc.).
    
    //================================================ Settings regarding the input files
    std::vector < std::string > inputFiles;           //!< This vector contains the filenames of all input structure files.
    bool forceP1;                                     //!< Should the P1 spacegroup be forced on the input PDB files?
    bool removeWaters;                                //!< Should all waters be removed from input PDB files?
    bool firstModelOnly;                              //!< Shoud only the first PDB model be used, or should all models be used?
    
    //================================================ Settings regarding the resolution of calculations
    proshade_single requestedResolution;              //!< The resolution to which the calculations are to be done.
    bool changeMapResolution;                         //!< Should maps be re-sampled to obtain the required resolution?
    bool changeMapResolutionTriLinear;                //!< Should maps be re-sampled to obtain the required resolution?
    
    //================================================ Settings regarding the PDB B-factor change
    proshade_double pdbBFactorNewVal;                 //!< Change all PDB B-factors to this value (for smooth maps).
    
    //================================================ Settings regarding the bandwidth of calculations
    proshade_unsign maxBandwidth;                     //!< The bandwidth of spherical harmonics decomposition for the largest sphere.
    proshade_double rotationUncertainty;              //!< Alternative to bandwidth - the angle in degrees to which the rotation function accuracy should be computed.

    //================================================ Settings regarding the phase
    bool usePhase;                                    //!< If true, the full data will be used, if false, Patterson maps will be used instead and phased data will be converted to them. Also, only half of the spherical harmonics bands will be necessary as odd bands have to be 0 for Patterson maps.
    
    //================================================ Settings regarding the spheres
    proshade_single maxSphereDists;                   //!< The distance between spheres in spherical mapping for the largest sphere.
    
    //================================================ Settings regarding the Gauss-Legendre integration
    proshade_unsign integOrder;                       //!< The order required for full Gauss-Legendre integration between the spheres.
    proshade_unsign taylorSeriesCap;                  //!< The max limit on the Taylor series expansion done for the abscissas of the Gauss-Legendre integration.
    
    //================================================ Settings regarding map normalisation
    bool normaliseMap;                                //!< Should the map be normalised to mean 0 sd 1?
    
    //================================================ Settings regarding map inversion
    bool invertMap;                                   //!< Should the map be inverted? Only use this if you think you have the wrong hand in your map.
    
    //================================================ Settings regarding map masking
    proshade_single blurFactor;                       //!< This is the amount by which B-factors should be increased to create the blurred map for masking.
    proshade_single maskingThresholdIQRs;             //!< Number of inter-quartile ranges from the median to be used for thresholding the blurred map for masking.
    bool maskMap;                                     //!< Should the map be masked from noise?
    bool useCorrelationMasking;                       //!< Should the blurring masking (false) or the correlation masking (true) be used?
    proshade_single halfMapKernel;                    //!< This value in Angstrom will be used as the kernel for the "fake half-map" computation.
    proshade_single correlationKernel;                //!< This value in Angstrom will be used as the kernel for the map-FHM correlation computation.
    bool saveMask;                                    //!< Should the mask be saved?
    std::string maskFileName;                         //!< The filename to which mask should be saved.
    
    //================================================ Settings regarding re-boxing
    bool reBoxMap;                                    //!< This switch decides whether re-boxing is needed.
    proshade_single boundsExtraSpace;                 //!< The number of extra angstroms to be added to all re-boxing bounds just for safety.
    proshade_signed boundsSimilarityThreshold;        //!< Number of indices which can be added just to make sure same size in indices is achieved.
    bool useSameBounds;                               //!< Switch to say that the same boundaries as  used for the first should be used for all input maps.
    proshade_signed* forceBounds;                     //!< These will be the boundaries to be forced upon the map.
    
    //================================================ Settings regarding COM
    bool moveToCOM;                                   //!< Logical value stating whether the structure should be moved to have its Centre Of Mass (COM) in the middle.
    
    //================================================ Settings regarding extra cell space
    proshade_single addExtraSpace;                    //!< If this value is non-zero, this many angstroms of empty space will be added to the internal map.
    
    //================================================ Settings regarding shell settings
    bool progressiveSphereMapping;                    //!< If true, each shell will have its own angular resolution dependent on the actual number of map points which are available to it. If false, all shells will have the same settings.
    
    //================================================ Settings regarding output file name
    std::string outName;                              //!< The file name where the output structure(s) should be saved.
    
    //================================================ Settings regarding distances computation
    bool computeEnergyLevelsDesc;                     //!< If true, the energy levels descriptor will be computed, otherwise all its computations will be omitted.
    proshade_double enLevMatrixPowerWeight;           //!< If RRP matrices shell position is to be weighted by putting the position as an exponent, this variable sets the exponent. Set to 0 for no weighting.
    bool computeTraceSigmaDesc;                       //!< If true, the trace sigma descriptor will be computed, otherwise all its computations will be omitted.
    bool computeRotationFuncDesc;                     //!< If true, the rotation function descriptor will be computed, otherwise all its computations will be omitted.
    
    //================================================ Settings regarding peak searching
    proshade_unsign peakNeighbours;                   //!< Number of points in any direction that have to be lower than the considered index in order to consider this index a peak.
    proshade_double noIQRsFromMedianNaivePeak;        //!< When doing peak searching, how many IQRs from the median the threshold for peak height should be (in terms of median of non-peak values).
    
    //================================================ Settings regarding 1D grouping
    proshade_double smoothingFactor;                  //!< This factor decides how small the group sizes should be - larger factor means more smaller groups.
    
    //================================================ Settings regarding the symmetry detection
    proshade_double symMissPeakThres;                 //!< Percentage of peaks that could be missing that would warrant starting the missing peaks search procedure.
    proshade_double axisErrTolerance;                 //!< Allowed error on vector axis in in dot product ( acos ( 1 - axErr ) is the allowed difference in radians ).
    bool axisErrToleranceDefault;
    proshade_double minSymPeak;                       //!< Minimum average peak for symmetry axis to be considered as "real".
    std::string recommendedSymmetryType;              //!< The symmetry type that ProSHADE finds the best fitting for the structure. Possible values are "" for none, "C" for cyclic, "D" for Dihedral, "T" for Tetrahedral, "O" for Octahedral and "I" for Icosahedral. C and D types also have fold value associated.
    proshade_unsign recommendedSymmetryFold;          //!< The fold of the recommended symmetry C or D type, 0 otherwise.
    std::string requestedSymmetryType;                //!< The symmetry  type requested by the user. Allowed values are C, D, T, O and I.
    proshade_unsign requestedSymmetryFold;            //!< The fold of the requested symmetry (only applicable to C and D symmetry types).
    bool usePeakSearchInRotationFunctionSpace;        //!< This variable switch decides whether symmetry detection will be done using peak search in rotation function or using the angle-axis sperical space.
    bool useBiCubicInterpolationOnPeaks;              //!< This variable switch decides whether best symmetry is detected from peak indices, or whether bicubic interpolation is done to seatch for better axis between indices.
    proshade_unsign maxSymmetryFold;                  //!< The highest symmetry fold to search for.
    
    //================================================ Settings regarding the structure overlay
    std::string overlayStructureName;                 //!< The filename to which the rotated and translated moving structure is to be saved.
    std::string rotTrsJSONFile;                       //!< The filename to which the rotation and translation operations are to be saved into.
    
    //================================================ Settings regarding verbosity of the program
    proshade_signed verbose;                          //!< Should the software report on the progress, or just be quiet? Value between -1 (nothing) and 4 (loud)
        
public:
    //================================================ Symmetry results holding values. This is required for Python being able to access the results without having the ProSHADE_run object.
    std::vector < proshade_double* > detectedSymmetry; //!< The vector of detected symmetry axes.
    std::vector < std::vector< proshade_double > > allDetectedCAxes; //!< The vector of all detected cyclic symmetry axes.
    std::vector < std::vector< proshade_unsign > > allDetectedDAxes; //!< The vector of all detected dihedral symmetry axes indices in allDetectedCAxes.
    std::vector < proshade_unsign > allDetectedTAxes; //!< The vector of all detected tetrahedral symmetry axes indices in allDetectedCAxes.
    std::vector < proshade_unsign > allDetectedOAxes; //!< The vector of all detected octahedral symmetry axes indices in allDetectedCAxes.
    std::vector < proshade_unsign > allDetectedIAxes; //!< The vector of all detected icosahedral symmetry axes indices in allDetectedCAxes.
    
public: // maybe make this protected?
    //================================================ Variable modifying functions
    void determineBandwidthFromAngle                  ( proshade_double uncertainty );
    void determineBandwidth                           ( proshade_unsign circumference );
    void determineSphereDistances                     ( proshade_single maxMapRange );
    void determineIntegrationOrder                    ( proshade_single maxMapRange );
    void determineAllSHValues                         ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_double xDimAngs,
                                                        proshade_double yDimAngs, proshade_double zDimAngs );
    void setVariablesLeftOnAuto                       ( void );
    
public:
    //================================================ Constructors / Destructors
    ProSHADE_settings                                 ( void );
    ProSHADE_settings                                 ( ProSHADE_Task task );
   ~ProSHADE_settings                                 ( void );
    
    //================================================ Variable setting functions
    void addStructure                                 ( std::string structure );
    void setResolution                                ( proshade_single resolution );
    void setPDBBFactor                                ( proshade_double newBF );
    void setNormalisation                             ( bool normalise );
    void setMapInversion                              ( bool mInv );
    void setVerbosity                                 ( proshade_signed verbosity );
    void setMaskBlurFactor                            ( proshade_single blurFac );
    void setMaskIQR                                   ( proshade_single noIQRs );
    void setMasking                                   ( bool mask );
    void setCorrelationMasking                        ( bool corMask );
    void setTypicalNoiseSize                          ( proshade_single typNoi );
    void setMinimumMaskSize                           ( proshade_single minMS );
    void setMaskSaving                                ( bool savMsk );
    void setMaskFilename                              ( std::string mskFln );
    void setMapReboxing                               ( bool reBx );
    void setBoundsSpace                               ( proshade_single boundsExSp );
    void setBoundsThreshold                           ( proshade_signed boundsThres );
    void setSameBoundaries                            ( bool sameB );
    void setOutputFilename                            ( std::string oFileName );
    void setMapResolutionChange                       ( bool mrChange );
    void setMapResolutionChangeTriLinear              ( bool mrChange );
    void setMapCentering                              ( bool com );
    void setExtraSpace                                ( proshade_single exSpace );
    void setBandwidth                                 ( proshade_unsign band );
    void setSphereDistances                           ( proshade_single sphDist );
    void setIntegrationOrder                          ( proshade_unsign intOrd );
    void setTaylorSeriesCap                           ( proshade_unsign tayCap );
    void setProgressiveSphereMapping                  ( bool progSphMap );
    void setEnergyLevelsComputation                   ( bool enLevDesc );
    void setTraceSigmaComputation                     ( bool trSigVal );
    void setRotationFunctionComputation               ( bool rotfVal );
    void setPeakNeighboursNumber                      ( proshade_unsign pkS );
    void setPeakNaiveNoIQR                            ( proshade_double noIQRs );
    void setPhaseUsage                                ( bool phaseUsage );
    void setEnLevShellWeight                          ( proshade_double mPower );
    void setGroupingSmoothingFactor                   ( proshade_double smFact );
    void setMissingPeakThreshold                      ( proshade_double mpThres );
    void setAxisComparisonThreshold                   ( proshade_double axThres );
    void setAxisComparisonThresholdBehaviour          ( bool behav );
    void setMinimumPeakForAxis                        ( proshade_double minSP );
    void setRecommendedSymmetry                       ( std::string val );
    void setRecommendedFold                           ( proshade_unsign val );
    void setRequestedSymmetry                         ( std::string val );
    void setRequestedFold                             ( proshade_unsign val );
    void setDetectedSymmetry                          ( proshade_double* sym );
    void setOverlaySaveFile                           ( std::string filename );
    void setOverlayJsonFile                           ( std::string filename );
    void setSymmetryRotFunPeaks                       ( bool rotFunPeaks );
    void setBicubicInterpolationSearch                ( bool bicubPeaks );
    void setMaxSymmetryFold                           ( proshade_unsign maxFold );
    
    //================================================ Command line options parsing
    void getCommandLineParams                         ( int argc, char** argv );
    
    //================================================ Debugging
    void printSettings                                ( void );
};

#endif
