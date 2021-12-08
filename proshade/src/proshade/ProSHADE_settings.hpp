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
    \version   0.7.6.2
    \date      DEC 2021
 */

//==================================================== ProSHADE
#include "ProSHADE_maths.hpp"

//==================================================== Overinclusion protection
#ifndef PROSHADE_SETTINGS
#define PROSHADE_SETTINGS

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
    bool removeNegativeDensity;                       //!< Should the negative density be removed from input files?
    
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
    std::string appliedMaskFileName;                  //!< The filename from which mask data will be read from.
    
    //================================================ Settings regarding map weighting
    std::string fourierWeightsFileName;               //!< The filename from which Fourier weights data will be read from.
    
    //================================================ Settings regarding re-boxing
    bool reBoxMap;                                    //!< This switch decides whether re-boxing is needed.
    proshade_single boundsExtraSpace;                 //!< The number of extra angstroms to be added to all re-boxing bounds just for safety.
    proshade_signed boundsSimilarityThreshold;        //!< Number of indices which can be added just to make sure same size in indices is achieved.
    bool useSameBounds;                               //!< Switch to say that the same boundaries as  used for the first should be used for all input maps.
    proshade_signed* forceBounds;                     //!< These will be the boundaries to be forced upon the map.
    
    //================================================ Settings regarding COM
    bool moveToCOM;                                   //!< Logical value stating whether the structure should be moved to have its Centre Of Mass (COM) in the middle.
    std::vector< proshade_double > boxCentre;         //!< If box centre is to be in any other location, this variable will hold the real space location that should be it.
    
    //================================================ Settings regarding extra cell space
    proshade_single addExtraSpace;                    //!< If this value is non-zero, this many angstroms of empty space will be added to the internal map.
    proshade_single coOrdsExtraSpace;                 //!< This number of Angstroms will be added before and any co-ordinates to make sure there is no atom directly at the edge of the map.
    
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
    bool findSymCentre;                               //!< Should phase-less map be used to determine centre of symmetry?
    proshade_double symMissPeakThres;                 //!< Percentage of peaks that could be missing that would warrant starting the missing peaks search procedure.
    proshade_double axisErrTolerance;                 //!< Allowed error on vector axis in in dot product ( acos ( 1 - axErr ) is the allowed difference in radians ).
    bool axisErrToleranceDefault;
    proshade_double minSymPeak;                       //!< Minimum average peak for symmetry axis to be considered as "real".
    std::string recommendedSymmetryType;              //!< The symmetry type that ProSHADE finds the best fitting for the structure. Possible values are "" for none, "C" for cyclic, "D" for Dihedral, "T" for Tetrahedral, "O" for Octahedral and "I" for Icosahedral. C and D types also have fold value associated.
    proshade_unsign recommendedSymmetryFold;          //!< The fold of the recommended symmetry C or D type, 0 otherwise.
    std::string requestedSymmetryType;                //!< The symmetry  type requested by the user. Allowed values are C, D, T, O and I.
    proshade_unsign requestedSymmetryFold;            //!< The fold of the requested symmetry (only applicable to C and D symmetry types).
    bool useBiCubicInterpolationOnPeaks;              //!< This variable switch decides whether best symmetry is detected from peak indices, or whether bicubic interpolation is done to seatch for better axis between indices.
    proshade_unsign maxSymmetryFold;                  //!< The highest symmetry fold to search for.
    proshade_double fscThreshold;                     //!< The threshold for FSC value under which the axis is considered to be likely noise.
    proshade_double peakThresholdMin;                 //!< The threshold for peak height above which axes are considered possible.
    bool fastISearch;                                 //!< Should FSC be computed for all possible I matches, or just for the best one according to FR?
    
    //================================================ Settings regarding centre of map
    std::vector < proshade_double > centrePosition;   //!< The position of the centre of the map in "real space" co-ordinates.
    
    //================================================ Settings regarding the structure overlay
    std::string overlayStructureName;                 //!< The filename to which the rotated and translated moving structure is to be saved.
    std::string rotTrsJSONFile;                       //!< The filename to which the rotation and translation operations are to be saved into.
    
    //================================================ Settings regarding verbosity of the program
    proshade_signed verbose;                          //!< Should the software report on the progress, or just be quiet? Value between -1 (nothing) and 4 (loud)
    proshade_signed messageShift;                     //!< This value allows shifting the messages to create more readable log for sub-processes.
        
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
    void determineAllSHValues                         ( proshade_unsign xDim, proshade_unsign yDim, proshade_single xDimAngs,
                                                        proshade_single yDimAngs, proshade_single zDimAngs );
    void setVariablesLeftOnAuto                       ( void );
    
public:
    //================================================ Constructors / Destructors
#if defined ( _WIN64 ) || defined ( _WIN32 )
    __declspec(dllexport)  ProSHADE_settings          ( );
    __declspec(dllexport)  ProSHADE_settings          ( ProSHADE_settings* settings );
    __declspec(dllexport)  ProSHADE_settings          ( ProSHADE_Task task );
    __declspec(dllexport) ~ProSHADE_settings          ( );
#else
    ProSHADE_settings                                 ( void );
    ProSHADE_settings                                 ( ProSHADE_settings* settings );
    ProSHADE_settings                                 ( ProSHADE_Task task );
   ~ProSHADE_settings                                 ( void );
#endif
    
    //================================================ Variable setting functions
#if defined ( _WIN64 ) || defined ( _WIN32 )
    void __declspec(dllexport) addStructure                                   ( std::string structure );
    void __declspec(dllexport) setResolution                                  ( proshade_single resolution );
    void __declspec(dllexport) setPDBBFactor                                  ( proshade_double newBF );
    void __declspec(dllexport) setNormalisation                               ( bool normalise );
    void __declspec(dllexport) setMapInversion                                ( bool mInv );
    void __declspec(dllexport) setVerbosity                                   ( proshade_signed verbosity );
    void __declspec(dllexport) setMaskBlurFactor                              ( proshade_single blurFac );
    void __declspec(dllexport) setMaskIQR                                     ( proshade_single noIQRs );
    void __declspec(dllexport) setMasking                                     ( bool mask );
    void __declspec(dllexport) setCorrelationMasking                          ( bool corMask );
    void __declspec(dllexport) setTypicalNoiseSize                            ( proshade_single typNoi );
    void __declspec(dllexport) setMinimumMaskSize                             ( proshade_single minMS );
    void __declspec(dllexport) setMaskSaving                                  ( bool savMsk );
    void __declspec(dllexport) setMaskFilename                                ( std::string mskFln );
    void __declspec(dllexport) setAppliedMaskFilename                         ( std::string mskFln );
    void __declspec(dllexport) setFourierWeightsFilename                      ( std::string fWgFln );
    void __declspec(dllexport) setMapReboxing                                 ( bool reBx );
    void __declspec(dllexport) setBoundsSpace                                 ( proshade_single boundsExSp );
    void __declspec(dllexport) setBoundsThreshold                             ( proshade_signed boundsThres );
    void __declspec(dllexport) setSameBoundaries                              ( bool sameB );
    void __declspec(dllexport) setOutputFilename                              ( std::string oFileName );
    void __declspec(dllexport) setMapResolutionChange                         ( bool mrChange );
    void __declspec(dllexport) setMapResolutionChangeTriLinear                ( bool mrChange );
    void __declspec(dllexport) setMapCentering                                ( bool com );
    void __declspec(dllexport) setExtraSpace                                  ( proshade_single exSpace );
    void __declspec(dllexport) setCoordExtraSpace                             ( proshade_single exSpace );
    void __declspec(dllexport) setBoxCentering                                ( proshade_double xPos, proshade_double yPos, proshade_double zPos );
    void __declspec(dllexport) setBandwidth                                   ( proshade_unsign band );
    void __declspec(dllexport) setSphereDistances                             ( proshade_single sphDist );
    void __declspec(dllexport) setIntegrationOrder                            ( proshade_unsign intOrd );
    void __declspec(dllexport) setTaylorSeriesCap                             ( proshade_unsign tayCap );
    void __declspec(dllexport) setProgressiveSphereMapping                    ( bool progSphMap );
    void __declspec(dllexport) setEnergyLevelsComputation                     ( bool enLevDesc );
    void __declspec(dllexport) setTraceSigmaComputation                       ( bool trSigVal );
    void __declspec(dllexport) setRotationFunctionComputation                 ( bool rotfVal );
    void __declspec(dllexport) setPeakNeighboursNumber                        ( proshade_unsign pkS );
    void __declspec(dllexport) setPeakNaiveNoIQR                              ( proshade_double noIQRs );
    void __declspec(dllexport) setPhaseUsage                                  ( bool phaseUsage );
    void __declspec(dllexport) setEnLevShellWeight                            ( proshade_double mPower );
    void __declspec(dllexport) setGroupingSmoothingFactor                     ( proshade_double smFact );
    void __declspec(dllexport) setSymmetryCentreSearch                        ( bool sCen );
    void __declspec(dllexport) setMissingPeakThreshold                        ( proshade_double mpThres );
    void __declspec(dllexport) setAxisComparisonThreshold                     ( proshade_double axThres );
    void __declspec(dllexport) setAxisComparisonThresholdBehaviour            ( bool behav );
    void __declspec(dllexport) setMinimumPeakForAxis                          ( proshade_double minSP );
    void __declspec(dllexport) setRecommendedSymmetry                         ( std::string val );
    void __declspec(dllexport) setRecommendedFold                             ( proshade_unsign val );
    void __declspec(dllexport) setRequestedSymmetry                           ( std::string val );
    void __declspec(dllexport) setRequestedFold                               ( proshade_unsign val );
    void __declspec(dllexport) setDetectedSymmetry                            ( proshade_double* sym );
    void __declspec(dllexport) setOverlaySaveFile                             ( std::string filename );
    void __declspec(dllexport) setOverlayJsonFile                             ( std::string filename );
    void __declspec(dllexport) setBicubicInterpolationSearch                  ( bool bicubPeaks );
    void __declspec(dllexport) setMaxSymmetryFold                             ( proshade_unsign maxFold );
    void __declspec(dllexport) setFSCThreshold                                ( proshade_double fscThr );
    void __declspec(dllexport) setPeakThreshold                               ( proshade_double peakThr );
    void __declspec(dllexport) setNegativeDensity                             ( bool nDens );
#else
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
    void setAppliedMaskFilename                       ( std::string mskFln );
    void setFourierWeightsFilename                    ( std::string fWgFln );
    void setMapReboxing                               ( bool reBx );
    void setBoundsSpace                               ( proshade_single boundsExSp );
    void setBoundsThreshold                           ( proshade_signed boundsThres );
    void setSameBoundaries                            ( bool sameB );
    void setOutputFilename                            ( std::string oFileName );
    void setMapResolutionChange                       ( bool mrChange );
    void setMapResolutionChangeTriLinear              ( bool mrChange );
    void setMapCentering                              ( bool com );
    void setExtraSpace                                ( proshade_single exSpace );
    void setCoordExtraSpace                           ( proshade_single exSpace );
    void setBoxCentering                              ( proshade_double xPos, proshade_double yPos, proshade_double zPos );
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
    void setSymmetryCentreSearch                      ( bool sCen );
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
    void setBicubicInterpolationSearch                ( bool bicubPeaks );
    void setMaxSymmetryFold                           ( proshade_unsign maxFold );
    void setFSCThreshold                              ( proshade_double fscThr );
    void setPeakThreshold                             ( proshade_double peakThr );
    void setNegativeDensity                           ( bool nDens );
#endif
    
    //================================================ Command line options parsing
#if defined ( _WIN64 ) || defined ( _WIN32 )
    void __declspec(dllexport) getCommandLineParams   ( int argc, char** argv );
#else
    void                       getCommandLineParams   ( int argc, char** argv );
#endif
    
    //================================================ Debugging
#if defined ( _WIN64 ) || defined ( _WIN32 )
    void __declspec(dllexport) printSettings          ( void );
#else
    void                       printSettings          ( void );
#endif
};

#endif
