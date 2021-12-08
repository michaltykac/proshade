/*! \file ProSHADE_data.hpp
    \brief This is the header file containing internal data representation and manipulation structures and functions.
 
    This header file contains the ProSHADE_data class declaration as well as the declarations for simple manipulations with the data
    (more complex manipulations are done in dedicated source files) and caller functions for the more complex manipulations.
    The class described here is how ProSHADE stores the structural data internally; however, the user should not need to access
    any of this code manually, as changes to this structure may have large consequences unforseen by the user.
   
    Copyright by Michal Tykac and individual contributors. All rights reserved.
     
    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.
     
    This software is provided by the copyright holders and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are  disclaimed. In     no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or  services, loss of use, data     or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this  software, even if advised of the possibility     of such damage.
     
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.6.2
    \date      DEC 2021
*/

//==================================================== ProSHADE
#include "ProSHADE_peakSearch.hpp"

//==================================================== Overinclusion protection
#ifndef PROSHADE_DATA
#define PROSHADE_DATA

//==================================================== ProSHADE_internal_data Namespace
/*! \namespace ProSHADE_internal_data
    \brief This namespace contains the structure and functions required for data reading and storing their derivates.
 
 
    The ProSHADE_internal_data namespace contains the data structure. It also has the data derivates storing variables,
    but it does not provide the computation code except for the forward declarations. The user should not need to access
    this namespace when using the library.
 */
namespace ProSHADE_internal_data
{    
/*! \class ProSHADE_data
    \brief This class contains all inputed and derived data for a single structure.

    This class codes the object that contains all the information about the input data and the derived information as well. It does not,
    however, provide the computation code as that lives elsewhere, except for the forward declarations.
 */
    class ProSHADE_data
    {
    public:
        //============================================ Variables regarding input file
        std::string fileName;                         //!< This is the original file from which the data were obtained.
        ProSHADE_internal_io::InputType fileType;     //!< This is the type of the input file.
        
        //============================================ Variables regarding map
        proshade_double* internalMap;                 //!< The internal map data representation, which may be amended as the run progresses.
        
        //============================================ Variables regarding map information
        proshade_single xDimSize;                     //!< This is the size of the map cell x dimension in Angstroms.
        proshade_single yDimSize;                     //!< This is the size of the map cell y dimension in Angstroms.
        proshade_single zDimSize;                     //!< This is the size of the map cell z dimension in Angstroms.
        proshade_single aAngle;                       //!< This is the angle a of the map cell in degrees.
        proshade_single bAngle;                       //!< This is the angle b of the map cell in degrees.
        proshade_single cAngle;                       //!< This is the angle c of the map cell in degrees.
        proshade_unsign xDimIndices;                  //!< This is the size of the map cell x dimension in indices.
        proshade_unsign yDimIndices;                  //!< This is the size of the map cell y dimension in indices.
        proshade_unsign zDimIndices;                  //!< This is the size of the map cell z dimension in indices.
        proshade_unsign xGridIndices;                 //!< As far as I know, this is identical to the xDimIndices.
        proshade_unsign yGridIndices;                 //!< As far as I know, this is identical to the yDimIndices.
        proshade_unsign zGridIndices;                 //!< As far as I know, this is identical to the zDimIndices.
        proshade_unsign xAxisOrder;                   //!< This is the order of the x axis.
        proshade_unsign yAxisOrder;                   //!< This is the order of the y axis.
        proshade_unsign zAxisOrder;                   //!< This is the order of the z axis.
        proshade_signed xAxisOrigin;                  //!< This is the origin position along the x axis.
        proshade_signed yAxisOrigin;                  //!< This is the origin position along the y axis.
        proshade_signed zAxisOrigin;                  //!< This is the origin position along the z axis.
        proshade_double xCom;                         //!< The COM of the map after processing along the X-axis.
        proshade_double yCom;                         //!< The COM of the map after processing along the Y-axis.
        proshade_double zCom;                         //!< The COM of the map after processing along the Z-axis.
        
        //============================================ Variables regarding original input values (i.e. these do not change with ProSHADE manipulations)
        proshade_single xDimSizeOriginal;             //!< This is the size of the map cell x dimension in Angstroms.
        proshade_single yDimSizeOriginal;             //!< This is the size of the map cell y dimension in Angstroms.
        proshade_single zDimSizeOriginal;             //!< This is the size of the map cell z dimension in Angstroms.
        proshade_unsign xDimIndicesOriginal;          //!< This is the size of the map cell x dimension in indices.
        proshade_unsign yDimIndicesOriginal;          //!< This is the size of the map cell y dimension in indices.
        proshade_unsign zDimIndicesOriginal;          //!< This is the size of the map cell z dimension in indices.
        proshade_signed xAxisOriginOriginal;          //!< This is the origin position along the x axis.
        proshade_signed yAxisOriginOriginal;          //!< This is the origin position along the y axis.
        proshade_signed zAxisOriginOriginal;          //!< This is the origin position along the z axis.
        proshade_double originalMapXCom;              //!< The COM of the first map to be loaded/computed without any furhter changes being reflacted along the X axis.
        proshade_double originalMapYCom;              //!< The COM of the first map to be loaded/computed without any furhter changes being reflacted along the Y axis.
        proshade_double originalMapZCom;              //!< The COM of the first map to be loaded/computed without any furhter changes being reflacted along the Z axis.
        proshade_double mapMovFromsChangeX;           //!< When the map is translated, the xFrom and xTo values are changed. This variable holds how much they have changed.
        proshade_double mapMovFromsChangeY;           //!< When the map is translated, the yFrom and yTo values are changed. This variable holds how much they have changed.
        proshade_double mapMovFromsChangeZ;           //!< When the map is translated, the zFrom and zTo values are changed. This variable holds how much they have changed.
        proshade_double mapCOMProcessChangeX;         //!< The change in X axis between the creation of the structure (originalMapXCom) and just before rotation.
        proshade_double mapCOMProcessChangeY;         //!< The change in Y axis between the creation of the structure (originalMapYCom) and just before rotation.
        proshade_double mapCOMProcessChangeZ;         //!< The change in Z axis between the creation of the structure (originalMapZCom) and just before rotation.
        
        //============================================ Variables regarding rotation and translation of original input files
        proshade_double originalPdbRotCenX;           //!< The centre of rotation as it relates to the original PDB positions (and not the ProSHADE internal map) along the x-axis.
        proshade_double originalPdbRotCenY;           //!< The centre of rotation as it relates to the original PDB positions (and not the ProSHADE internal map) along the y-axis.
        proshade_double originalPdbRotCenZ;           //!< The centre of rotation as it relates to the original PDB positions (and not the ProSHADE internal map) along the z-axis.
        proshade_double originalPdbTransX;            //!< The optimal translation vector as it relates to the original PDB positions (and not the ProSHADE internal map) along the x-axis.
        proshade_double originalPdbTransY;            //!< The optimal translation vector as it relates to the original PDB positions (and not the ProSHADE internal map) along the y-axis.
        proshade_double originalPdbTransZ;            //!< The optimal translation vector as it relates to the original PDB positions (and not the ProSHADE internal map) along the z-axis.

        //============================================ Variables regarding iterator positions
        proshade_signed xFrom;                        //!< This is the starting index along the x axis.
        proshade_signed yFrom;                        //!< This is the starting index along the y axis.
        proshade_signed zFrom;                        //!< This is the starting index along the z axis.
        proshade_signed xTo;                          //!< This is the final index along the x axis.
        proshade_signed yTo;                          //!< This is the final index along the y axis.
        proshade_signed zTo  ;                        //!< This is the final index along the z axis.
        
        //============================================ Variables regarding SH mapping spheres
        std::vector<proshade_single> spherePos;       //!< Vector of sphere radii from the centre of the map.
        proshade_unsign noSpheres;                    //!< The number of spheres with map projected onto them.
        ProSHADE_internal_spheres::ProSHADE_sphere** spheres; //!< The set of concentric spheres to which the intermal density map has been projected.
        proshade_complex** sphericalHarmonics;        //!< A set of spherical harmonics values arrays for each sphere.
        proshade_complex** rotSphericalHarmonics;     //!< A set of rotated spherical harmonics values arrays for each sphere, used only if map rotation is required.
        proshade_unsign maxShellBand;                 //!< The maximum band for any shell of the object.
        
        //============================================ Variables regarding shape distance computations
        proshade_double*** rrpMatrices;               //!< The energy levels descriptor shell correlation tables.
        proshade_complex*** eMatrices;                //!< The trace sigma and full rotation function c*conj(c) integral tables.
        proshade_double integrationWeight;            //!< The Pearson's c.c. type weighting for the integration.
        proshade_complex* so3Coeffs;                  //!< The coefficients obtained by SO(3) Fourier Transform (SOFT), in this case derived from the E matrices.
        proshade_complex* so3CoeffsInverse;           //!< The inverse coefficients obtained by inverse SO(3) Fourier Transform (SOFT) - i.e. rotation function.
        proshade_complex*** wignerMatrices;           //!< These matrices are computed for a particular rotation to be done in spherical harmonics
        proshade_unsign maxCompBand;                  //!< The largest comparison band - this variable tells how large arrays will be allocated for the comparison.
        proshade_complex* translationMap;             //!< This is where the translation map will be held, if at all used.
        
        //============================================ Variables regarding symmetry detection
        std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*> sphereMappedRotFun;
        
        //============================================ Control variables
        bool isEmpty;                                 //!< This variable stated whether the class contains any information.
        proshade_unsign inputOrder;                   //!< This value is the input order - it is useful to know for writing out files, so that they would not overwrite the same name multiple times.
        
    protected:
        void figureIndexStartStop                     ( void );
        void setPDBMapValues                          ( void );
        void readInMAP                                ( ProSHADE_settings* settings, proshade_double* maskArr = nullptr, proshade_unsign maskXDim = 0, proshade_unsign maskYDim = 0,
                                                        proshade_unsign maskZDim = 0, proshade_double* weightsArr = nullptr, proshade_unsign weigXDim = 0, proshade_unsign weigYDim = 0,
                                                        proshade_unsign weigZDim = 0 );
        void readInPDB                                ( ProSHADE_settings* settings );
        void readInGemmi                              ( gemmi::Structure gemmiStruct, ProSHADE_settings* settings );
        void allocateRRPMemory                        ( );
        
    public:
        //============================================ Constructors / Destructors
        ProSHADE_data                                 ( );
        ProSHADE_data                                 ( std::string strName, double *mapVals, int len, proshade_single xDmSz, proshade_single yDmSz,
                                                        proshade_single zDmSz, proshade_unsign xDmInd, proshade_unsign yDmInd, proshade_unsign zDmInd, proshade_signed xFr,
                                                        proshade_signed yFr, proshade_signed zFr, proshade_signed xT, proshade_signed yT, proshade_signed zT,
                                                        proshade_unsign inputO );
       ~ProSHADE_data                                 ( void );
        
        //============================================ Data I/O functions
        void readInStructure                          ( std::string fName, proshade_unsign inputO, ProSHADE_settings* settings, proshade_double* maskArr = nullptr, proshade_unsign maskXDim = 0,
                                                        proshade_unsign maskYDim = 0, proshade_unsign maskZDim = 0, proshade_double* weightsArr = nullptr, proshade_unsign weigXDim = 0,
                                                        proshade_unsign weigYDim = 0, proshade_unsign weigZDim = 0 );
        void readInStructure                          ( gemmi::Structure gemmiStruct, proshade_unsign inputO, ProSHADE_settings* settings );
        void writeMap                                 ( std::string fName, std::string title = "Created by ProSHADE and written by GEMMI", int mode = 2 );
        void writePdb                                 ( std::string fName, proshade_double euA = 0.0, proshade_double euB = 0.0, proshade_double euG = 0.0,
                                                        proshade_double trsX = 0.0, proshade_double trsY = 0.0, proshade_double trsZ = 0.0, proshade_double rotX = 0.0,
                                                        proshade_double rotY = 0.0, proshade_double rotZ = 0.0, bool firstModel = true );
        void writeGemmi                               ( std::string fName, gemmi::Structure gemmiStruct, proshade_double euA = 0.0, proshade_double euB = 0.0, proshade_double euG = 0.0,
                                                        proshade_double trsX = 0.0, proshade_double trsY = 0.0, proshade_double trsZ = 0.0, proshade_double rotX = 0.0,
                                                        proshade_double rotY = 0.0, proshade_double rotZ = 0.0, bool firstModel = true );
        void writeMask                                ( std::string fName, proshade_double* mask );
        
        //============================================ Data processing functions
        void invertMirrorMap                          ( ProSHADE_settings* settings );
        void normaliseMap                             ( ProSHADE_settings* settings );
        void maskMap                                  ( ProSHADE_settings* settings );
        void getReBoxBoundaries                       ( ProSHADE_settings* settings, proshade_signed*& ret );
        void createNewMapFromBounds                   ( ProSHADE_settings* settings, ProSHADE_data*& newStr, proshade_signed* newBounds );
        void reSampleMap                              ( ProSHADE_settings* settings );
        void centreMapOnCOM                           ( ProSHADE_settings* settings );
        void addExtraSpace                            ( ProSHADE_settings* settings );
        void removePhaseInormation                    ( ProSHADE_settings* settings );
        void shiftToBoxCentre                         ( ProSHADE_settings* settings );
        void shiftToRotationCentre                    ( ProSHADE_settings* settings );
        void processInternalMap                       ( ProSHADE_settings* settings );
        
        //============================================ Data sphere mapping functions
        void getSpherePositions                       ( ProSHADE_settings* settings );
        void mapToSpheres                             ( ProSHADE_settings* settings );
        void computeSphericalHarmonics                ( ProSHADE_settings* settings );
        
        //============================================ Distances pre-computation functions
        bool shellBandExists                          ( proshade_unsign shell, proshade_unsign bandVal );
        void computeRRPMatrices                       ( ProSHADE_settings* settings );
        void allocateEMatrices                        ( proshade_unsign  band );
        void allocateSO3CoeffsSpace                   ( proshade_unsign band );
        void allocateWignerMatricesSpace              ( );
        
        //============================================ Symmetry detection functions
        void computeRotationFunction                  ( ProSHADE_settings* settings );
        void convertRotationFunction                  ( ProSHADE_settings* settings );
        void getRealEMatrixValuesForLM                ( proshade_signed band, proshade_signed order1, double *eMatsLMReal, int len );
        void getImagEMatrixValuesForLM                ( proshade_signed band, proshade_signed order1, double *eMatsLMImag, int len );
        void getRealSO3Coeffs                         ( double *so3CoefsReal, int len );
        void getImagSO3Coeffs                         ( double *so3CoefsImag, int len );
        void getRealRotFunction                       ( double *rotFunReal, int len );
        void getImagRotFunction                       ( double *rotFunImag, int len );
        void getRealTranslationFunction               ( double *trsFunReal, int len );
        void getImagTranslationFunction               ( double *trsFunImag, int len );
        void getRotMatrixFromRotFunInds               ( proshade_signed aI, proshade_signed bI, proshade_signed gI, double *rotMat, int len );
        int so3CoeffsArrayIndex                       ( proshade_signed order1, proshade_signed order2, proshade_signed band );
        std::vector< proshade_double* > getDihedralSymmetriesList    ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList );
        std::vector< proshade_double* > getTetrahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList );
        std::vector< proshade_double* > getOctahedralSymmetriesList  ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList );
        std::vector< proshade_double* > getIcosahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList );
        std::vector < std::vector< proshade_double* > > getPredictedIcosahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList );
        std::vector< proshade_double* > getPredictedOctahedralSymmetriesList  ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList );
        std::vector< proshade_double* > getPredictedTetrahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList );
        void detectSymmetryFromAngleAxisSpace         ( ProSHADE_settings* settings, std::vector< proshade_double* >* axes, std::vector < std::vector< proshade_double > >* allCs );
        std::vector< proshade_double* > getCyclicSymmetriesListFromAngleAxis ( ProSHADE_settings* settings );
        std::vector< proshade_double* > findRequestedCSymmetryFromAngleAxis  ( ProSHADE_settings* settings, proshade_unsign fold, proshade_double* peakThres );
        void saveDetectedSymmetries                   ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSyms, std::vector < std::vector< proshade_double > >* allCs );
        std::string     getRecommendedSymmetryType    ( ProSHADE_settings* settings );
        proshade_unsign getRecommendedSymmetryFold    ( ProSHADE_settings* settings );
        proshade_unsign getNoRecommendedSymmetryAxes  ( ProSHADE_settings* settings );
        std::vector< std::string > getSymmetryAxis    ( ProSHADE_settings* settings, proshade_unsign axisNo );
        void prepareFSCFourierMemory                  ( fftw_complex*& mapData, fftw_complex*& origCoeffs, fftw_complex*& fCoeffs, proshade_signed*& binIndexing,
                                                        proshade_signed* noBins, proshade_double**& bindata, proshade_signed*& binCounts, fftw_plan* planForwardFourier, proshade_double*& fscByBin );
        proshade_double computeFSC                    ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSym, size_t symIndex, fftw_complex* mapData,
                                                        fftw_complex* fCoeffs, fftw_complex* origCoeffs, fftw_plan* planForwardFourier, proshade_signed noBins,
                                                        proshade_signed *binIndexing, proshade_double**& bindata, proshade_signed*& binCounts, proshade_double*& fscByBin );
        proshade_double computeFSC                    ( ProSHADE_settings* settings, proshade_double* sym, fftw_complex* mapData, fftw_complex* fCoeffs, fftw_complex* origCoeffs,
                                                        fftw_plan* planForwardFourier, proshade_signed noBins, proshade_signed *binIndexing, proshade_double**& bindata,
                                                        proshade_signed*& binCounts, proshade_double*& fscByBin );
        void saveRecommendedSymmetry                  ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSym, std::vector< proshade_double* >* DSym,
                                                        std::vector< proshade_double* >* TSym, std::vector< proshade_double* >* OSym,
                                                        std::vector< proshade_double* >* ISym, std::vector< proshade_double* >* axes, fftw_complex* mapData,
                                                        fftw_complex* origCoeffs, fftw_complex* fCoeffs, fftw_plan* planForwardFourier, proshade_signed noBins, proshade_signed* binIndexing,
                                                        proshade_double** bindata, proshade_signed* binCounts, proshade_double*& fscByBin );
        void saveRequestedSymmetryC                   ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSym, std::vector< proshade_double* >* axes );
        void saveRequestedSymmetryD                   ( ProSHADE_settings* settings, std::vector< proshade_double* >* DSym, std::vector< proshade_double* >* axes, fftw_complex* mapData,
                                                        fftw_complex* origCoeffs, fftw_complex* fCoeffs, fftw_plan* planForwardFourier, proshade_signed noBins, proshade_signed* binIndexing,
                                                        proshade_double** bindata, proshade_signed* binCounts, proshade_double*& fscByBin );
        std::vector<std::vector< proshade_double > > getAllGroupElements ( ProSHADE_settings* settings, std::vector< proshade_unsign > axesList, std::string groupType = "", proshade_double matrixTolerance = 0.05 );
        std::vector<std::vector< proshade_double > > getAllGroupElements ( std::vector < std::vector< proshade_double > >* allCs, std::vector< proshade_unsign > axesList, std::string groupType = "", proshade_double matrixTolerance = 0.05 );
        void reportSymmetryResults                    ( ProSHADE_settings* settings );
        
        //============================================ Map overlay functions
        void getOverlayRotationFunction               ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* obj2 );
        std::vector< proshade_double > getBestRotationMapPeaksEulerAngles ( ProSHADE_settings* settings );
        std::vector< proshade_double > getBestTranslationMapPeaksAngstrom ( ProSHADE_internal_data::ProSHADE_data* staticStructure );
        void zeroPaddToDims                           ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim );
        void rotateMapReciprocalSpace                 ( ProSHADE_settings* settings, proshade_double eulerAlpha, proshade_double eulerBeta, proshade_double eulerGamma );
        std::vector< proshade_double > rotateMapRealSpace ( proshade_double axX, proshade_double axY, proshade_double axZ, proshade_double axAng, proshade_double*& map );
        std::vector< proshade_double > rotateMapRealSpaceInPlace ( proshade_double eulA, proshade_double eulB, proshade_double eulG );
        void translateMap                             ( proshade_double trsX, proshade_double trsY, proshade_double trsZ );
        void allocateRotatedSHMemory                  ( void );
        void computeRotatedSH                         ( void );
        void invertSHCoefficients                     ( void );
        void interpolateMapFromSpheres                ( proshade_double*& densityMapRotated );
        void computeTranslationMap                    ( ProSHADE_internal_data::ProSHADE_data* obj1 );
        void findMapCOM                               ( void );
        void writeOutOverlayFiles                     ( ProSHADE_settings* settings, proshade_double eulA, proshade_double eulB, proshade_double eulG, std::vector< proshade_double >* rotCentre,
                                                        std::vector< proshade_double >* ultimateTranslation );
        void reportOverlayResults                     ( ProSHADE_settings* settings, std::vector < proshade_double >* rotationCentre, std::vector < proshade_double >* eulerAngles,
                                                        std::vector < proshade_double >* finalTranslation );
        
        //============================================ Python access functions
        void deepCopyMap                              ( proshade_double*& saveTo, proshade_signed verbose );
        
        //============================================ Accessor functions
        proshade_double getMapValue                   ( proshade_unsign pos );
        proshade_unsign getMaxSpheres                 ( void );
        proshade_unsign getMaxBand                    ( void );
        proshade_double* getRealSphHarmValue          ( proshade_unsign band, proshade_unsign order, proshade_unsign shell );
        proshade_double* getImagSphHarmValue          ( proshade_unsign band, proshade_unsign order, proshade_unsign shell );
        proshade_double getRRPValue                   ( proshade_unsign band, proshade_unsign sh1, proshade_unsign sh2 );
        proshade_double getAnySphereRadius            ( proshade_unsign shell );
        proshade_double getIntegrationWeight          ( void );
        proshade_unsign getShellBandwidth             ( proshade_unsign shell );
        proshade_single getSpherePosValue             ( proshade_unsign shell );
        proshade_complex** getEMatrixByBand           ( proshade_unsign band );
        void getEMatrixValue                          ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_double* valueReal, proshade_double* valueImag );
        proshade_complex* getInvSO3Coeffs             ( void );
        proshade_complex* getSO3Coeffs                ( void );
        proshade_unsign getComparisonBand             ( void );
        void getWignerMatrixValue                     ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_double* valueReal, proshade_double* valueImag );
        proshade_single getXDimSize                   ( void );
        proshade_single getYDimSize                   ( void );
        proshade_single getZDimSize                   ( void );
        proshade_unsign getXDim                       ( void );
        proshade_unsign getYDim                       ( void );
        proshade_unsign getZDim                       ( void );
        proshade_signed* getXFromPtr                  ( void );
        proshade_signed* getYFromPtr                  ( void );
        proshade_signed* getZFromPtr                  ( void );
        proshade_signed* getXToPtr                    ( void );
        proshade_signed* getYToPtr                    ( void );
        proshade_signed* getZToPtr                    ( void );
        proshade_signed* getXAxisOrigin               ( void );
        proshade_signed* getYAxisOrigin               ( void );
        proshade_signed* getZAxisOrigin               ( void );
        proshade_double*& getInternalMap              ( void );
        proshade_complex* getTranslationFnPointer     ( void );
        std::vector< proshade_double > getMapCOMProcessChange ( void );
        
        //============================================ Mutator functions
        void setIntegrationWeight                     ( proshade_double intW );
        void setIntegrationWeightCumul                ( proshade_double intW );
        void setEMatrixValue                          ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_complex val );
        void normaliseEMatrixValue                    ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_double normF );
        void setSO3CoeffValue                         ( proshade_unsign position, proshade_complex val );
        void setWignerMatrixValue                     ( proshade_complex val, proshade_unsign band, proshade_unsign order1, proshade_unsign order2 );
    };

    //================================================ Support functions
    std::vector<std::vector< proshade_double > > computeGroupElementsForGroup ( proshade_double xAx, proshade_double yAx, proshade_double zAx, proshade_signed fold );

    std::vector<std::vector< proshade_double > > joinElementsFromDifferentGroups ( std::vector<std::vector< proshade_double > >* first,
                                                                                   std::vector<std::vector< proshade_double > >* second,
                                                                                   proshade_double matrixTolerance,
                                                                                   bool combine );
}

#endif
