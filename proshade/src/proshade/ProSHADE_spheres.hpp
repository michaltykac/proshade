/*! \file ProSHADE_spheres.hpp
    \brief This header file contains the declarations for the ProSHADE_sphere class.
 
    The ProSHADE_sphere class and the ProSHADE_internal_spheres namespace are both declated in this header file. The functionalities that these two objects
    provide to ProSHADE are mapping a 3D map onto a surface of a sphere with a given radius, object for holding such a mapping and onto which the spherical
    harmonics decomposition could be done as well as support function for data access and processing.
 
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
#include "ProSHADE_sphericalHarmonics.hpp"

//==================================================== Overinclusion protection
#ifndef PROSHADE_SPHERES
#define PROSHADE_SPHERES

//==================================================== ProSHADE_internal_spheres Namespace
/*! \namespace ProSHADE_internal_spheres
    \brief This namespace contains the structure and functions required for storing internal map projections onto a set of concentric spheres.
 
 
    The ProSHADE_internal_spheres namespace contains the the structure and functions required to map the internal map data onto a set of concentric
    spheres as required by ProSHADE. It also has the ability to store sphere specific versions of the values such as bandwidth and angular resolution,
    as well as the functionality for computing spherical harmonics. Finally, the automatic spherical harmonics computation variables determination function
    live here.
 */
namespace ProSHADE_internal_spheres
{
    /*! \class ProSHADE_sphere
     \brief This class contains all inputed and derived data for a single sphere.
     
     This class codes the object that contains all the information about a single concentric sphere as well as all the functionality required
     to process such data.
     */
    class ProSHADE_sphere
    {
    private:
        //============================================ General sphere spherical harmonics computation variables
        proshade_unsign localBandwidth;               //!< The bandwidth of spherical harmonics decomposition for thihs sphere.
        proshade_unsign localAngRes;                  //!< The maximum number of indices on the angular grids forming the spheres for this sphere.
        proshade_single sphereWidth;                  //!< The width of this sphere (the same as sphere distance).
        proshade_double sphereRadius;                 //!< The distance from map centre to the shell.
        
        //============================================ General map fraction values
        proshade_single maxSphereRange;               //!< The maximum sphere diagonal distance.
        proshade_unsign shellOrder;                   //!< The order of the shell in the shell vector - 0 is the smallest shell.
        
        //============================================ General map dealing values
        proshade_single xDimSampling;                 //!< The sampling rate (distance between indices) of the x axis.
        proshade_single yDimSampling;                 //!< The sampling rate (distance between indices) of the y axis.
        proshade_single zDimSampling;                 //!< The sampling rate (distance between indices) of the z axis.
        
        //============================================ Mapped data
        proshade_double* mappedData;                  //!< The mapped density data.
        proshade_double* mappedDataRot;               //!< The rotated structure mapped data.
        
    protected:
        proshade_unsign getMaxCircumference           ( proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax, proshade_single maxRange );
        bool getMapPoint                              ( proshade_double* map, proshade_unsign xDimMax, proshade_unsign yDimMax,
                                                        proshade_unsign zDimMax, proshade_signed xPos,
                                                        proshade_signed yPos, proshade_signed zPos, std::vector<proshade_double>* interpVec );
        void getLongitudeCutoffs                      ( std::vector<proshade_double>* lonCO );
        void getLattitudeCutoffs                      ( std::vector<proshade_double>* latCO );
        void getInterpolationXYZ                      ( proshade_double* x, proshade_double* y, proshade_double* z,
                                                        proshade_double thetaIt, std::vector<proshade_double>* lonCO,
                                                        proshade_unsign phiIt, std::vector<proshade_double>* latCO );
        void getXYZTopBottoms                         ( proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax,
                                                        proshade_double x, proshade_double y, proshade_double z,
                                                        proshade_signed* xBottom, proshade_signed* yBottom, proshade_signed* zBottom,
                                                        proshade_signed* xTop,
                                                        proshade_signed* yTop, proshade_signed* zTop );
        void interpolateAlongFirst                    ( std::vector<proshade_double> c000, std::vector<proshade_double> c001, std::vector<proshade_double> c010,
                                                        std::vector<proshade_double> c011, std::vector<proshade_double> c100, std::vector<proshade_double> c101,
                                                        std::vector<proshade_double> c110, std::vector<proshade_double> c111,
                                                        std::vector<proshade_double>* c00, std::vector<proshade_double>* c01,
                                                        std::vector<proshade_double>* c10, std::vector<proshade_double>* c11,
                                                        proshade_double xd );
        void interpolateAlongSecond                   ( std::vector<proshade_double> c00, std::vector<proshade_double> c01,
                                                        std::vector<proshade_double> c10, std::vector<proshade_double> c11,
                                                        std::vector<proshade_double>* c0, std::vector<proshade_double>* c1, proshade_double yd );
        void mapData                                  ( proshade_double* map, proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax );
        
    public:
        //============================================ Constructors / Destructors
        ProSHADE_sphere                               ( proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax, proshade_single xSize,
                                                        proshade_single ySize, proshade_single zSize, proshade_unsign shOrder,
                                                        std::vector<proshade_single>* spherePos,
                                                        bool progressiveMapping, proshade_unsign band, proshade_double* map,
                                                        proshade_unsign* maxShellBand );
       ~ProSHADE_sphere                               ( void );
        
        //============================================ Rotated map data related functions
        void allocateRotatedMap                       ( void );
        
        //============================================ Accessor/Mutator functions
        proshade_unsign getLocalBandwidth             ( void );
        proshade_unsign getLocalAngRes                ( void );
        proshade_double* getMappedData                ( void );
        proshade_double getShellRadius                ( void );
        void setRotatedMappedData                     ( proshade_unsign pos, proshade_double value );
        proshade_double getRotatedMappedData          ( proshade_unsign pos );
    };
    
    //================================================ Sphere resolution and computation decisions
    proshade_unsign autoDetermineBandwidth            ( proshade_unsign circumference );
    proshade_single autoDetermineSphereDistances      ( proshade_single maxMapRange, proshade_single resolution );
    proshade_unsign autoDetermineIntegrationOrder     ( proshade_single maxMapRange, proshade_single sphereDist );

/*! \class ProSHADE_rotFun_spherePeakGroup
    \brief This class contains peak groups detected in the rotation function mapped spheres.
 
    This class codes the object that contains all the information about a single group of peaks found in the set of
    ProSHADE_rotFun_sphere objects with mapped rotation function values.
 */
    class ProSHADE_rotFun_spherePeakGroup
    {
    public:
        proshade_double latSampling;
        proshade_double lonSampling;
        proshade_unsign dimension;
        proshade_double latFrom;
        proshade_double latTo;
        proshade_double lonFrom;
        proshade_double lonTo;
        proshade_double latFromInds;
        proshade_double latToInds;
        proshade_double lonFromInds;
        proshade_double lonToInds;
        
        std::vector<proshade_unsign> spherePositions;
        
        proshade_double* latMinLonMinXYZ;
        proshade_double* latMaxLonMinXYZ;
        proshade_double* latMinLonMaxXYZ;
        proshade_double* latMaxLonMaxXYZ;
        
    protected:
        void computeCornerPositions                   ( void );
        proshade_signed angularDistanceWithBorders    ( proshade_signed origLat, proshade_signed testedLat );
        void getAllAngleDifferences                   ( std::vector< proshade_double >* angDiffs, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*> sphereVals );
        void getAllPossibleFolds                      ( std::vector< proshade_double >* angDiffs, std::vector< proshade_unsign >* foldsToTry );
        void getSpheresFormingFold                    ( proshade_unsign foldToTry, std::vector< proshade_unsign >* spheresFormingFold,
                                                        std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*> sphereVals, proshade_double sphereAngleTolerance );
        void getBestIndexForFold                      ( proshade_double* bestPosVal, proshade_double* bestLatInd, proshade_double* bestLonInd, std::vector< proshade_unsign >* spheresFormingFold,
                                                        std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*> sphereVals );
        
    public:
        ProSHADE_rotFun_spherePeakGroup               ( proshade_double lat, proshade_double lon, proshade_unsign sphPos, proshade_unsign angDim );
       ~ProSHADE_rotFun_spherePeakGroup               ( void );
        
    public:
        bool checkIfPeakBelongs                       ( proshade_double lat, proshade_double lon, proshade_unsign sphPos, proshade_double cosTol, proshade_signed verbose, proshade_signed messageShift );
        void findCyclicPointGroupsGivenFold           ( std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*> sphereVals, std::vector < proshade_double* >* detectedCs,
                                                        bool bicubicInterp, proshade_unsign fold, proshade_signed verbose, proshade_signed messageShift );
        
    public:
        proshade_double getLatFromIndices             ( void );
        proshade_double getLatToIndices               ( void );
        proshade_double getLonFromIndices             ( void );
        proshade_double getLonToIndices               ( void );
        std::vector<proshade_unsign> getSpherePositions ( void );
    };
}


#endif
