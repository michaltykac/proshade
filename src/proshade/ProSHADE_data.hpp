/*! \file ProSHADE_data.hpp
 \brief ...
 
 ...
 
 This file is part of the ProSHADE library for calculating
 shape descriptors and symmetry operators of protein structures.
 This is a prototype code, which is by no means complete or fully
 tested. Its use is at your own risk only. There is no quarantee
 that the results are correct.
 
 \author    Michal Tykac
 \author    Garib N. Murshudov
 \version   0.7.2
 \date      DEC 2019
 */

//============================================ ProSHADE
#include "ProSHADE_spheres.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE_DATA__
#define __PROSHADE_DATA__

//============================================ ProSHADE_internal_data Namespace
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
        //==================================== Variables regarding input file
        std::string                           fileName;                 //!< This is the original file from which the data were obtained.
        ProSHADE_internal_io::InputType       fileType;                 //!< This is the type of the input file.
        
        //==================================== Variables regarding map
        proshade_double*                      internalMap;              //!< The internal map data representation, which may be amended as the run progresses.
        
        //==================================== Variables regarding map information
        proshade_single                       xDimSize;                 //!< This is the size of the map cell x dimension in Angstroms.
        proshade_single                       yDimSize;                 //!< This is the size of the map cell y dimension in Angstroms.
        proshade_single                       zDimSize;                 //!< This is the size of the map cell z dimension in Angstroms.
        proshade_single                       aAngle;                   //!< This is the angle a of the map cell in degrees.
        proshade_single                       bAngle;                   //!< This is the angle b of the map cell in degrees.
        proshade_single                       cAngle;                   //!< This is the angle c of the map cell in degrees.
        proshade_unsign                       xDimIndices;              //!< This is the size of the map cell x dimension in indices.
        proshade_unsign                       yDimIndices;              //!< This is the size of the map cell y dimension in indices.
        proshade_unsign                       zDimIndices;              //!< This is the size of the map cell z dimension in indices.
        proshade_unsign                       xGridIndices;             //!< As far as I know, this is identical to the xDimIndices.
        proshade_unsign                       yGridIndices;             //!< As far as I know, this is identical to the yDimIndices.
        proshade_unsign                       zGridIndices;             //!< As far as I know, this is identical to the zDimIndices.
        proshade_unsign                       xAxisOrder;               //!< This is the order of the x axis.
        proshade_unsign                       yAxisOrder;               //!< This is the order of the y axis.
        proshade_unsign                       zAxisOrder;               //!< This is the order of the z axis.
        proshade_signed                       xAxisOrigin;              //!< This is the origin position along the x axis.
        proshade_signed                       yAxisOrigin;              //!< This is the origin position along the y axis.
        proshade_signed                       zAxisOrigin;              //!< This is the origin position along the z axis.
        proshade_double                       comMovX;                  //!< How much the map was moved along the X-axis to achieve COM centering.
        proshade_double                       comMovY;                  //!< How much the map was moved along the X-axis to achieve COM centering.
        proshade_double                       comMovZ;                  //!< How much the map was moved along the X-axis to achieve COM centering.

        //==================================== Variables regarding iterator positions
        proshade_signed                       xFrom;                    //!< This is the starting index along the x axis.
        proshade_signed                       yFrom;                    //!< This is the starting index along the y axis.
        proshade_signed                       zFrom;                    //!< This is the starting index along the z axis.
        proshade_signed                       xTo;                      //!< This is the final index along the x axis.
        proshade_signed                       yTo;                      //!< This is the final index along the y axis.
        proshade_signed                       zTo  ;                    //!< This is the final index along the z axis.
        
        //==================================== Variables regarding SH mapping spheres
        std::vector<proshade_single>          spherePos;                //!< Vector of sphere radii from the centre of the map.
        proshade_unsign                       noSpheres;                //!< The number of spheres with map projected onto them.
        ProSHADE_internal_spheres::ProSHADE_sphere** spheres;           //!< The set of concentric spheres to which the intermal density map has been projected.
        proshade_complex**                    sphericalHarmonics;       //!< A set of spherical harmonics values arrays for each sphere.
        proshade_complex**                    rotSphericalHarmonics;    //!< A set of rotated spherical harmonics values arrays for each sphere, used only if map rotation is required.
        proshade_unsign                       maxShellBand;             //!< The maximum band for any shell of the object.
        
        //==================================== Variables regarding shape distance computations
        proshade_double***                    rrpMatrices;              //!< The energy levels descriptor shell correlation tables.
        proshade_complex***                   eMatrices;                //!< The trace sigma and full rotation function c*conj(c) integral tables.
        proshade_double                       integrationWeight;        //!< The Pearson's c.c. type weighting for the integration.
        proshade_complex*                     so3Coeffs;                //!< The coefficients obtained by SO(3) Fourier Transform (SOFT), in this case derived from the E matrices.
        proshade_complex*                     so3CoeffsInverse;         //!< The inverse coefficients obtained by inverse SO(3) Fourier Transform (SOFT) - i.e. rotation function.
        proshade_complex***                   wignerMatrices;           //!< These matrices are computed for a particular rotation to be done in spherical harmonics
        proshade_unsign                       maxCompBand;              //!< The largest comparison band - this variable tells how large arrays will be allocated for the comparison.
        proshade_complex*                     translationMap;           //!< This is where the translation map will be held, if at all used.
        
        //==================================== Control variables
        bool                                  isEmpty;                  //!< This variable stated whether the class contains any information.
        proshade_unsign                       inputOrder;               //!< This value is the input order - it is useful to know for writing out files, so that they would not overwrite the same name multiple times.
        
    protected:
        void figureIndexStartStop             ( void );
        void switchAxes                       ( void );
        void setPDBMapValues                  ( void );
        void readInMAP                        ( ProSHADE_settings* settings );
        void readInPDB                        ( ProSHADE_settings* settings );
        void allocateRRPMemory                ( ProSHADE_settings* settings );
        
    public:
        //==================================== Constructors / Destructors
        ProSHADE_data                         ( ProSHADE_settings* settings );
        ProSHADE_data                         ( ProSHADE_settings* settings, std::string strName, double *mapVals, int len, proshade_single xDmSz, proshade_single yDmSz,
                                                proshade_single zDmSz, proshade_unsign xDmInd, proshade_unsign yDmInd, proshade_unsign zDmInd, proshade_signed xFr,
                                                proshade_signed yFr, proshade_signed zFr, proshade_signed xT, proshade_signed yT, proshade_signed zT,
                                                proshade_unsign inputO );
       ~ProSHADE_data                         ( );
        
        //==================================== Data I/O functions
        void readInStructure                  ( std::string fName, proshade_unsign inputO, ProSHADE_settings* settings );
        void writeMap                         ( std::string fName, std::string title = "" );
        void writeMask                        ( std::string fName, proshade_double* mask );
        int getMapArraySizePython             ( void );
        void getMapPython                     ( double *mapArrayPython, int len );
        void setMapPython                     ( double *mapChangedInPython, int len );
        void setNewMapPython                  ( double *mapChangedInPython, int len );
        
        //==================================== Data processing functions
        void invertMirrorMap                  ( ProSHADE_settings* settings );
        void normaliseMap                     ( ProSHADE_settings* settings );
        void maskMap                          ( ProSHADE_settings* settings );
        void getReBoxBoundaries               ( ProSHADE_settings* settings, proshade_signed*& ret );
        void getReBoxBoundariesPy             ( ProSHADE_settings* settings, int* reBoxBounds, int len );
        void createNewMapFromBounds           ( ProSHADE_settings* settings, ProSHADE_data*& newStr, proshade_signed* newBounds );
        void createNewMapFromBoundsPy         ( ProSHADE_settings* settings, ProSHADE_data* newStr, int* newBounds, int len );
        void reSampleMap                      ( ProSHADE_settings* settings );
        void centreMapOnCOM                   ( ProSHADE_settings* settings );
        void addExtraSpace                    ( ProSHADE_settings* settings );
        void removePhaseInormation            ( ProSHADE_settings* settings );
        void processInternalMap               ( ProSHADE_settings* settings );
        
        //==================================== Data sphere mapping functions
        void getSpherePositions               ( ProSHADE_settings* settings );
        void mapToSpheres                     ( ProSHADE_settings* settings );
        void computeSphericalHarmonics        ( ProSHADE_settings* settings );
        void getRealSphericalHarmonicsForShell ( proshade_unsign shellNo, proshade_signed verbose, double *sphericalHarmsReal, int len );
        void getImagSphericalHarmonicsForShell ( proshade_unsign shellNo, proshade_signed verbose, double *sphericalHarmsImag, int len );
        int sphericalHarmonicsIndex           ( proshade_signed order, proshade_signed band, proshade_signed shell );
        int getSphericalHarmonicsLenForShell  ( proshade_unsign shellNo, proshade_signed verbose );
        
        //==================================== Distances pre-computation functions
        bool shellBandExists                  ( proshade_unsign shell, proshade_unsign bandVal );
        void computeRRPMatrices               ( ProSHADE_settings* settings );
        void allocateEMatrices                ( ProSHADE_settings* settings, proshade_unsign  band );
        void allocateSO3CoeffsSpace           ( proshade_unsign band );
        void allocateWignerMatricesSpace      ( ProSHADE_settings* settings );
        
        //==================================== Symmetry detection functions
        void getRotationFunction              ( ProSHADE_settings* settings );
        void getRealEMatrixValuesForLM        ( proshade_signed band, proshade_signed order1, double *eMatsLMReal, int len );
        void getImagEMatrixValuesForLM        ( proshade_signed band, proshade_signed order1, double *eMatsLMImag, int len );
        void getRealSO3Coeffs                 ( double *so3CoefsReal, int len );
        void getImagSO3Coeffs                 ( double *so3CoefsImag, int len );
        void getRealRotFunction               ( double *rotFunReal, int len );
        void getImagRotFunction               ( double *rotFunImag, int len );
        void getRealTranslationFunction       ( double *trsFunReal, int len );
        void getImagTranslationFunction       ( double *trsFunImag, int len );
        void getRotMatrixFromRotFunInds       ( proshade_signed aI, proshade_signed bI, proshade_signed gI, double *rotMat, int len );
        int so3CoeffsArrayIndex               ( proshade_signed order1, proshade_signed order2, proshade_signed band );
        std::vector< proshade_double* > getCyclicSymmetriesList ( ProSHADE_settings* settings );
        std::vector< proshade_double* > getDihedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList );
        std::vector< proshade_double* > getTetrahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList );
        std::vector< proshade_double* > getOctahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList );
        std::vector< proshade_double* > getIcosahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList );
        void detectSymmetryInStructure        ( ProSHADE_settings* settings, std::vector< proshade_double* >* axes );
        void detectSymmetryInStructurePython  ( ProSHADE_settings* settings );
        std::string     getRecommendedSymmetryType ( ProSHADE_settings* settings );
        proshade_unsign getRecommendedSymmetryFold ( ProSHADE_settings* settings );
        proshade_unsign getNoSymmetryAxes     ( ProSHADE_settings* settings );
        std::vector< std::string > getSymmetryAxis ( ProSHADE_settings* settings, proshade_unsign axisNo );
        proshade_double findBestCScore        ( std::vector< proshade_double* >* CSym, proshade_unsign* symInd );
        proshade_double findBestDScore        ( std::vector< proshade_double* >* DSym, proshade_unsign* symInd );
        proshade_double findTScore            ( std::vector< proshade_double* >* TSym );
        proshade_double findOScore            ( std::vector< proshade_double* >* OSym );
        proshade_double findIScore            ( std::vector< proshade_double* >* ISym );
        void saveRecommendedSymmetry          ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSym, std::vector< proshade_double* >* DSym,
                                                std::vector< proshade_double* >* TSym, std::vector< proshade_double* >* OSym,
                                                std::vector< proshade_double* >* ISym, std::vector< proshade_double* >* axes );
        void saveRequestedSymmetryC           ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSym, std::vector< proshade_double* >* axes );
        void saveRequestedSymmetryD           ( ProSHADE_settings* settings, std::vector< proshade_double* >* DSym, std::vector< proshade_double* >* axes );
        void reportSymmetryResults            ( ProSHADE_settings* settings );
        
        //==================================== Map overlay functions
        void getOverlayRotationFunction       ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* obj2 );
        std::vector< proshade_double > getBestRotationMapPeaksEulerAngles ( ProSHADE_settings* settings );
        std::vector< proshade_double > getBestTranslationMapPeaksAngstrom ( ProSHADE_internal_data::ProSHADE_data* staticStructure );
        void zeroPaddToDims                   ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim );
        void rotateMap                        ( ProSHADE_settings* settings, proshade_double eulerAlpha, proshade_double eulerBeta, proshade_double eulerGamma );
        void translateMap                     ( ProSHADE_settings* settings, proshade_double trsX, proshade_double trsY, proshade_double trsZ );
        void allocateRotatedSHMemory          ( ProSHADE_settings* settings );
        void computeRotatedSH                 ( ProSHADE_settings* settings );
        void inverseSHCoefficients            ( );
        void interpolateMapFromSpheres        ( ProSHADE_settings* settings, proshade_double*& densityMapRotated );
        void computeTranslationMap            ( ProSHADE_internal_data::ProSHADE_data* obj1 );
        
        //==================================== Python access functions
        void deepCopyMap                      ( proshade_double*& saveTo, proshade_unsign verbose );
        
        //==================================== Accessor functions
        proshade_double getMapValue           ( proshade_unsign pos );
        proshade_unsign getMaxSpheres         ( void );
        proshade_unsign getMaxBand            ( void );
        proshade_double* getRealSphHarmValue  ( proshade_unsign band, proshade_unsign order, proshade_unsign shell );
        proshade_double* getImagSphHarmValue  ( proshade_unsign band, proshade_unsign order, proshade_unsign shell );
        proshade_double getRRPValue           ( proshade_unsign band, proshade_unsign sh1, proshade_unsign sh2 );
        proshade_double getAnySphereRadius    ( proshade_unsign shell );
        proshade_double getIntegrationWeight  ( void );
        proshade_unsign getShellBandwidth     ( proshade_unsign shell );
        proshade_double getSpherePosValue     ( proshade_unsign shell );
        proshade_complex** getEMatrixByBand   ( proshade_unsign band );
        void getEMatrixValue                  ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_double* valueReal, proshade_double* valueImag );
        proshade_complex* getInvSO3Coeffs     ( void );
        proshade_complex* getSO3Coeffs        ( void );
        proshade_unsign getComparisonBand     ( void );
        void getWignerMatrixValue             ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_double* valueReal, proshade_double* valueImag );
        proshade_single getXDimSize           ( void );
        proshade_single getYDimSize           ( void );
        proshade_single getZDimSize           ( void );
        proshade_unsign getXDim               ( void );
        proshade_unsign getYDim               ( void );
        proshade_unsign getZDim               ( void );
        proshade_signed* getXFromPtr          ( void );
        proshade_signed* getYFromPtr          ( void );
        proshade_signed* getZFromPtr          ( void );
        proshade_signed* getXToPtr            ( void );
        proshade_signed* getYToPtr            ( void );
        proshade_signed* getZToPtr            ( void );
        proshade_signed* getXAxisOrigin       ( void );
        proshade_signed* getYAxisOrigin       ( void );
        proshade_signed* getZAxisOrigin       ( void );
        proshade_double*& getInternalMap      ( void );
        proshade_complex* getTranslationFnPointer ( void );
        
        //==================================== Mutator functions
        void setIntegrationWeight             ( proshade_double intW );
        void setIntegrationWeightCumul        ( proshade_double intW );
        void setEMatrixValue                  ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_complex val );
        void normaliseEMatrixValue            ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_double normF );
        void setSO3CoeffValue                 ( proshade_unsign position, proshade_complex val );
        void setWignerMatrixValue             ( proshade_complex val, proshade_unsign band, proshade_unsign order1, proshade_unsign order2 );
    };
}

#endif
