/*! \file ProSHADE_spheres.hpp
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
#include "ProSHADE_sphericalHarmonics.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE_SPHERES__
#define __PROSHADE_SPHERES__

//============================================ ProSHADE_internal_spheres Namespace
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
        //==================================== General sphere spherical harmonics computation variables
        proshade_unsign                       localBandwidth;            //!< The bandwidth of spherical harmonics decomposition for thihs sphere.
        proshade_unsign                       localAngRes;               //!< The maximum number of indices on the angular grids forming the spheres for this sphere.
        proshade_single                       sphereWidth;               //!< The width of this sphere (the same as sphere distance).
        proshade_double                       sphereRadius;              //!< The distance from map centre to the shell.
        
        //==================================== General map fraction values
        proshade_single                       maxSphereRange;            //!< The maximum sphere diagonal distance.
        proshade_unsign                       shellOrder;                //!< The order of the shell in the shell vector - 0 is the smallest shell.
        
        //==================================== General map dealing values
        proshade_single                       xDimSampling;              //!< The sampling rate (distance between indices) of the x axis.
        proshade_single                       yDimSampling;              //!< The sampling rate (distance between indices) of the y axis.
        proshade_single                       zDimSampling;              //!< The sampling rate (distance between indices) of the z axis.
        
        //==================================== Mapped data
        proshade_double*                      mappedData;                //!< The mapped density data.
        proshade_double*                      mappedDataRot;             //!< The rotated structure mapped data.
        
    protected:
        proshade_unsign getMaxCircumference   ( proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax, proshade_double maxRange,
                                                proshade_single xSize, proshade_single ySize, proshade_single zSize );
        bool getMapPoint                      ( proshade_double* map, proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax, proshade_signed xPos,
                                                proshade_signed yPos, proshade_signed zPos, std::vector<proshade_double>* interpVec );
        void getLongitudeCutoffs              ( std::vector<proshade_double>* lonCO );
        void getLattitudeCutoffs              ( std::vector<proshade_double>* latCO );
        void getInterpolationXYZ              ( proshade_double* x, proshade_double* y, proshade_double* z, proshade_double thetaIt, std::vector<proshade_double>* lonCO,
                                                proshade_unsign phiIt, std::vector<proshade_double>* latCO );
        void getXYZTopBottoms                 ( proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax, proshade_double x, proshade_double y,
                                                proshade_double z, proshade_signed* xBottom, proshade_signed* yBottom, proshade_signed* zBottom, proshade_signed* xTop,
                                                proshade_signed* yTop, proshade_signed* zTop );
        void interpolateAlongFirst            ( std::vector<proshade_double> c000, std::vector<proshade_double> c001, std::vector<proshade_double> c010, std::vector<proshade_double> c011,
                                                std::vector<proshade_double> c100, std::vector<proshade_double> c101, std::vector<proshade_double> c110, std::vector<proshade_double> c111,
                                                std::vector<proshade_double>* c00, std::vector<proshade_double>* c01, std::vector<proshade_double>* c10, std::vector<proshade_double>* c11,
                                                proshade_double xd );
        void interpolateAlongSecond           ( std::vector<proshade_double> c00, std::vector<proshade_double> c01, std::vector<proshade_double> c10, std::vector<proshade_double> c11,
                                                std::vector<proshade_double>* c0, std::vector<proshade_double>* c1, proshade_double yd );
        void mapData                          ( proshade_double* map, proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax );
        
    public:
        //==================================== Constructors / Destructors
        ProSHADE_sphere                       ( proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax, proshade_single xSize,
                                                proshade_single ySize, proshade_single zSize, proshade_unsign shOrder, std::vector<proshade_single>* spherePos,
                                                bool progressiveMapping, proshade_unsign band, proshade_unsign angR, proshade_double* map,
                                                proshade_unsign* maxShellBand );
       ~ProSHADE_sphere                       ( );
        
        //==================================== Rotated map data related functions
        void allocateRotatedMap               ( void );
        
        //==================================== Accessor/Mutator functions
        proshade_unsign getLocalBandwidth     ( void );
        proshade_unsign getLocalAngRes        ( void );
        proshade_double* getMappedData        ( void );
        proshade_double getShellRadius        ( void );
        void setRotatedMappedData             ( proshade_unsign pos, proshade_double value );
        proshade_double getRotatedMappedData  ( proshade_unsign pos );
    };
    
    //======================================== Sphere resolution and computation decisions
    proshade_unsign autoDetermineBandwidth    ( proshade_unsign circumference );
    proshade_unsign autoDetermineAngularResolution ( proshade_unsign circumference );
    proshade_single autoDetermineSphereDistances ( proshade_single maxMapRange, proshade_single resolution );
    proshade_unsign autoDetermineIntegrationOrder ( proshade_single maxMapRange, proshade_single sphereDist );
}

#endif
