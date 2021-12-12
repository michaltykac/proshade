/*! \file ProSHADE_spheres.cpp
    \brief This source file contains function related to the ProSHADE_sphere class, which generally serve to prepare a shell for SH decomposition.
 
    This source file contains function related to the ProSHADE_sphere class and the ProSHADE_internal_spheres namespace. These functions provide a density map
    mapping onto a sphere (one particular sphere with a particular radius) and also do all the preparations and allow the information access required for spherical harmonics
    decomposition to be done on this shell.
 
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

//============================================ ProSHADE
#include "ProSHADE_spheres.hpp"

/*! \brief Constructor for getting empty ProSHADE_sphere class.
 
    This constructor takes all the information required by the shell and proceeds to set the shell object accordingly. It can either set the
    spherical harmonics variables according to shell size (if the progressive mapping is allowed), or it just uses the whole map defaults (in
    the opposite case). It also maps the density onto the spherical grid of the shell using tri-linear interpolation.
 
    \param[in] xDimMax The internal map maximum index in the x dimension.
    \param[in] yDimMax The internal map maximum index in the y dimension.
    \param[in] zDimMax The internal map maximum index in the z dimension.
    \param[in] xSize The integral map size in angstroms along the x-axis.
    \param[in] ySize The integral map size in angstroms along the y-axis.
    \param[in] zSize The integral map size in angstroms along the z-axis.
    \param[in] shOrder The order of the shell in the shell vector - 0 is the smallest shell.
    \param[in] spherePos A Pointer to vector of proshade_single's with the shell positions.
    \param[in] progressiveMapping Should the bandwidth and angular resolution be set by the actual number of points, or just same for all shells?
    \param[in] band The bandwidth to be set for conservative mapping.
    \param[in] map A pointer to the internal map which should be mapped to the sphere.
    \param[in] maxShellBand This pointer reference will take the shell band, if it is higher than the already known one.
    \param[out] X Data object with all values set and the appropriate map part mapped.
 */
ProSHADE_internal_spheres::ProSHADE_sphere::ProSHADE_sphere ( proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax, proshade_single xSize, proshade_single ySize, proshade_single zSize, proshade_unsign shOrder, std::vector<proshade_single>* spherePos, bool progressiveMapping, proshade_unsign band, proshade_double* map, proshade_unsign* maxShellBand )
{
    //================================================ Save inputs
    this->shellOrder                                  = shOrder;
    this->sphereWidth                                 = static_cast< proshade_single > ( ( spherePos->at(0) + spherePos->at(1) ) / 2.0f );
    this->sphereRadius                                = static_cast< proshade_double > ( spherePos->at(shOrder) );
    
    //================================================ Determine shell ranges in angstroms
    proshade_double maxDist                           = 0.0;
    if ( shOrder == static_cast<proshade_unsign> ( spherePos->size() - 1 ) ) { maxDist = static_cast<proshade_double> ( spherePos->at(spherePos->size()-1) + ( spherePos->at(1) - spherePos->at(0) ) ); }
    else { maxDist = static_cast<proshade_double> ( ( spherePos->at(shOrder) + spherePos->at(shOrder+1) ) / 2.0f ); }
    
    //================================================ Set the max range
    this->maxSphereRange                              = static_cast<proshade_single> ( 2.0 * maxDist );
    
    //================================================ Set map sampling rates
    this->xDimSampling                                = xSize / static_cast<proshade_single> (xDimMax);
    this->yDimSampling                                = ySize / static_cast<proshade_single> (yDimMax);
    this->zDimSampling                                = zSize / static_cast<proshade_single> (zDimMax);
    
    //================================================ Get maximum circumference
    proshade_unsign maxCircumference                  = this->getMaxCircumference ( xDimMax, yDimMax, zDimMax, this->maxSphereRange );
    
    //================================================ Get spherical harmonics calculation values
    if ( progressiveMapping )
    {
        this->localBandwidth                          = std::min ( autoDetermineBandwidth ( maxCircumference ), band );
        this->localAngRes                             = this->localBandwidth * 2;
    }
    else
    {
        this->localBandwidth                          = band;
        this->localAngRes                             = this->localBandwidth * 2;
    }
    
    //================================================ Save the maximum shell band for later
    if ( *maxShellBand < this->localBandwidth ) {  *maxShellBand  = this->localBandwidth; }
    
    //================================================ Allocate memory for sphere mapping
    this->mappedData                                  = new proshade_double[this->localAngRes * this->localAngRes];
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->mappedData, __FILE__, __LINE__, __func__ );
    
    //================================================ Set rotated mapped data
    this->mappedDataRot                               = nullptr;
    
    //================================================ Map the data to the sphere
    this->mapData                                     ( map, xDimMax, yDimMax, zDimMax );
    
}

/*! \brief Destructor for the ProSHADE_sphere class.
 
    This destructor is responsible for releasing all memory used by the data storing object
 
    \param[out] X N/A.
 */
ProSHADE_internal_spheres::ProSHADE_sphere::~ProSHADE_sphere ( )
{
    delete[] this->mappedData;
    
    if ( this->mappedDataRot != nullptr )
    {
        delete[] this->mappedDataRot;
    }
}

/*! \brief This function determines the maximum circumference of a shell.
 
    This function takes the input variables and finds the maximum circumference of the sphere data.
 
    \param[in] xDimMax The internal map maximum index in the x dimension.
    \param[in] yDimMax The internal map maximum index in the y dimension.
    \param[in] zDimMax The internal map maximum index in the z dimension.
    \param[in] maxRange The maximum range of shell in angstroms.
 */
proshade_unsign ProSHADE_internal_spheres::ProSHADE_sphere::getMaxCircumference ( proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax, proshade_single maxRange )
{
    //================================================ Find from and to limits on indices
    proshade_signed xFromInd                          = static_cast<proshade_signed> ( xDimMax / 2 ) -
                                                        static_cast<proshade_signed> ( (maxRange/2) / this->xDimSampling );
    proshade_signed yFromInd                          = static_cast<proshade_signed> ( yDimMax / 2 ) -
                                                        static_cast<proshade_signed> ( (maxRange/2) / this->yDimSampling );
    proshade_signed zFromInd                          = static_cast<proshade_signed> ( zDimMax / 2 ) -
                                                        static_cast<proshade_signed> ( (maxRange/2) / this->zDimSampling );
            
    proshade_signed xToInd                            = static_cast<proshade_signed> ( xDimMax / 2 ) +
                                                        static_cast<proshade_signed> ( (maxRange/2) / this->xDimSampling );
    proshade_signed yToInd                            = static_cast<proshade_signed> ( yDimMax / 2 ) +
                                                        static_cast<proshade_signed> ( (maxRange/2) / this->yDimSampling );
    proshade_signed zToInd                            = static_cast<proshade_signed> ( zDimMax / 2 ) +
                                                        static_cast<proshade_signed> ( (maxRange/2) / this->zDimSampling );
    
    //================================================ Check for bounds
    if ( xFromInd < 0 ) { xFromInd = 0; }
    if ( yFromInd < 0 ) { yFromInd = 0; }
    if ( zFromInd < 0 ) { zFromInd = 0; }
    if ( xToInd > static_cast<proshade_signed> ( xDimMax ) ) { xToInd = static_cast<proshade_signed> ( xDimMax ); }
    if ( yToInd > static_cast<proshade_signed> ( yDimMax ) ) { yToInd = static_cast<proshade_signed> ( yDimMax ); }
    if ( zToInd > static_cast<proshade_signed> ( zDimMax ) ) { zToInd = static_cast<proshade_signed> ( zDimMax ); }
    
    //================================================ Get dim sizes
    proshade_unsign xDimSZ                            = static_cast< proshade_unsign > ( xToInd - xFromInd );
    proshade_unsign yDimSZ                            = static_cast< proshade_unsign > ( yToInd - yFromInd );
    proshade_unsign zDimSZ                            = static_cast< proshade_unsign > ( zToInd - zFromInd );
    
    //================================================ check for too sparse sampling
    if ( ( xDimSZ == 0 ) && ( yDimSZ == 0 ) && ( zDimSZ == 0 ) ) { xDimSZ = 1; yDimSZ = 1; zDimSZ = 1; }
    
    //================================================ Find max and mid dims
    std::vector<proshade_unsign> dimSizes;
    ProSHADE_internal_misc::addToUnsignVector         ( &dimSizes, xDimSZ );
    ProSHADE_internal_misc::addToUnsignVector         ( &dimSizes, yDimSZ );
    ProSHADE_internal_misc::addToUnsignVector         ( &dimSizes, zDimSZ );
    std::sort                                         ( dimSizes.begin(), dimSizes.end() );
    proshade_unsign maxDim                            = dimSizes.at(2);
    proshade_unsign midDim                            = dimSizes.at(1);
    
    //================================================ Return max circumference
    return                                            ( maxDim + midDim );
}

/*! \brief This function maps the internal map to the specific sphere.
 
    This function deals with the map values mapping onto a sphere in abstract terms, or in programatical terms it simply interpolates
    the appropriate map data onto a spherical grid a given by the sphere object calling it. It uses the tri-linear interpolation in the
    order XYZ and saves all results internally.
 
    \param[in] map The density map pointer - the map which should be mapped to the data.
    \param[in] xDimMax The internal map maximum index in the x dimension.
    \param[in] yDimMax The internal map maximum index in the y dimension.
    \param[in] zDimMax The internal map maximum index in the z dimension.
 */
void ProSHADE_internal_spheres::ProSHADE_sphere::mapData ( proshade_double* map, proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax )
{
    //================================================ Initialise local variavles
    proshade_double x, y, z, xRelative, yRelative, zRelative;
    proshade_signed xBottom, yBottom, zBottom, xTop, yTop, zTop;
    std::vector<proshade_double> lonCO                ( this->localAngRes + 1 );
    std::vector<proshade_double> latCO                ( this->localAngRes + 1 );
    std::vector<proshade_double> c000                 = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c001                 = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c010                 = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c011                 = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c100                 = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c101                 = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c110                 = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c111                 = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c00                  = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c01                  = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c10                  = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c11                  = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c0                   = std::vector<proshade_double> ( 4 );
    std::vector<proshade_double> c1                   = std::vector<proshade_double> ( 4 );
    
    //================================================ Find pixelisation cutOffs
    this->getLongitudeCutoffs                         ( &lonCO );
    this->getLattitudeCutoffs                         ( &latCO );
    
    //================================================ Interpolate the map onto this shell
    for ( unsigned int thIt = 0; thIt < this->localAngRes; thIt++ )
    {
        for ( unsigned int phIt = 0; phIt < this->localAngRes; phIt++ )
        {
            //======================================== Get grid point x, y and z
            this->getInterpolationXYZ                 ( &x, &y, &z, thIt, &lonCO, phIt, &latCO );
            
            //======================================== Find 8 closest point around the grid point in indices
            this->getXYZTopBottoms                    ( xDimMax, yDimMax, zDimMax, x, y, z, &xBottom, &yBottom, &zBottom, &xTop, &yTop, &zTop );

            //======================================== Get the 8 closest points interpolation vectors (or set to 0 if out-of-bounds)
            if ( !getMapPoint ( map, xDimMax, yDimMax, zDimMax, xBottom, yBottom, zBottom, &c000 ) ) { this->mappedData[phIt * this->localAngRes + thIt] = 0.0; continue; }
            if ( !getMapPoint ( map, xDimMax, yDimMax, zDimMax, xBottom, yBottom, zTop   , &c001 ) ) { this->mappedData[phIt * this->localAngRes + thIt] = 0.0; continue; }
            if ( !getMapPoint ( map, xDimMax, yDimMax, zDimMax, xBottom, yTop   , zBottom, &c010 ) ) { this->mappedData[phIt * this->localAngRes + thIt] = 0.0; continue; }
            if ( !getMapPoint ( map, xDimMax, yDimMax, zDimMax, xBottom, yTop   , zTop   , &c011 ) ) { this->mappedData[phIt * this->localAngRes + thIt] = 0.0; continue; }
            if ( !getMapPoint ( map, xDimMax, yDimMax, zDimMax, xTop   , yBottom, zBottom, &c100 ) ) { this->mappedData[phIt * this->localAngRes + thIt] = 0.0; continue; }
            if ( !getMapPoint ( map, xDimMax, yDimMax, zDimMax, xTop   , yBottom, zTop   , &c101 ) ) { this->mappedData[phIt * this->localAngRes + thIt] = 0.0; continue; }
            if ( !getMapPoint ( map, xDimMax, yDimMax, zDimMax, xTop   , yTop   , zBottom, &c110 ) ) { this->mappedData[phIt * this->localAngRes + thIt] = 0.0; continue; }
            if ( !getMapPoint ( map, xDimMax, yDimMax, zDimMax, xTop   , yTop   , zTop   , &c111 ) ) { this->mappedData[phIt * this->localAngRes + thIt] = 0.0; continue; }

            //======================================== Interpolate along X axis
            xRelative                                 = ( x - ( static_cast< proshade_double > ( xBottom - static_cast<proshade_signed> ( ( ( xDimMax ) / 2 ) ) ) * static_cast< proshade_double > ( this->xDimSampling ) ) ) / static_cast< proshade_double > ( this->xDimSampling );
            this->interpolateAlongFirst               ( c000, c001, c010, c011, c100, c101, c110, c111, &c00, &c01, &c10, &c11, xRelative );
            
            //======================================== Interpolate along Y axis
            yRelative                                 = ( y - ( static_cast< proshade_double > ( yBottom - static_cast<proshade_signed> ( ( ( yDimMax ) / 2 ) ) ) * static_cast< proshade_double > ( this->yDimSampling ) ) ) / static_cast< proshade_double > ( this->yDimSampling );
            this->interpolateAlongSecond              ( c00, c01, c10, c11, &c0, &c1, yRelative );

            //======================================== Save the resulting value
            zRelative                                 = ( z - ( static_cast< proshade_double > ( zBottom - static_cast<proshade_signed> ( ( ( zDimMax ) / 2 ) ) ) * static_cast< proshade_double > ( this->zDimSampling ) ) ) / static_cast< proshade_double > ( this->zDimSampling );
            this->mappedData[phIt * this->localAngRes + thIt] = ( c0.at(3) * ( 1.0 - zRelative ) ) + ( c1.at(3) * zRelative );
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function fills in the interpolation vector for a single map point.
 
    This function takes the index positions of a map point and a pointer to a vector of 4 proshade_double's and fills in the
    vector so that interpolation can easily be done using its values. It also returns boolean value, true if the point is
    in bounds and everything had succeeded and false if the point was out of bounds and therefore the interpolation is not
    required as the mapped value should just be 0.
 
    \param[in] map The density map pointer - the map which should be mapped to the data.
    \param[in] xDimMax The internal map maximum index in the x dimension.
    \param[in] yDimMax The internal map maximum index in the y dimension.
    \param[in] zDimMax The internal map maximum index in the z dimension.
    \param[in] xPos The x axis index position of the point of which interpolation vector is to be obtained.
    \param[in] yPos The y axis index position of the point of which interpolation vector is to be obtained.
    \param[in] zPos The z axis index position of the point of which interpolation vector is to be obtained.
    \param[in] interpVec The variable to which the interpolation values are to be saved (vector of proshade_double's initialised to length 4).
 */
bool ProSHADE_internal_spheres::ProSHADE_sphere::getMapPoint ( proshade_double* map, proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax, proshade_signed xPos, proshade_signed yPos, proshade_signed zPos, std::vector<proshade_double>* interpVec )
{
    //================================================ Initialise local variables
    proshade_signed posIter                           = zPos + static_cast<proshade_signed> ( zDimMax ) * ( yPos  + static_cast<proshade_signed> ( yDimMax ) * xPos );
    
    //================================================ Check of out-of-bounds shells
    if ( ( xPos < 0 ) || ( xPos >= static_cast<proshade_signed> ( xDimMax ) ) ) { return ( false ); }
    if ( ( yPos < 0 ) || ( yPos >= static_cast<proshade_signed> ( yDimMax ) ) ) { return ( false ); }
    if ( ( zPos < 0 ) || ( zPos >= static_cast<proshade_signed> ( zDimMax ) ) ) { return ( false ); }
    
    //================================================ Get the interpolation values
    interpVec->at(0)                                  = static_cast<proshade_double> ( xPos ) * static_cast<proshade_double> ( this->xDimSampling );
    interpVec->at(1)                                  = static_cast<proshade_double> ( yPos ) * static_cast<proshade_double> ( this->yDimSampling );
    interpVec->at(2)                                  = static_cast<proshade_double> ( zPos ) * static_cast<proshade_double> ( this->zDimSampling );
    interpVec->at(3)                                  = map[posIter];
    
    //================================================ Done
    return                                            ( true );
    
}

/*! \brief This function fills in the vector of longitudal bin boarder values.
 
    This is a simple abstraction function which finds the boarders of the angular grid binds along the longitudal angle.
    
    \param[in] lonCO A pointer for vector of proshade_double's to which the cutt-off values are to be saved. The vector needs to be initialised for this->localAngRes + 1 size.
*/
void ProSHADE_internal_spheres::ProSHADE_sphere::getLongitudeCutoffs ( std::vector<proshade_double>* lonCO )
{
    //================================================ Get the cut-offs
    for ( proshade_unsign iter = 0; iter <= this->localAngRes; iter++ )
    {
        lonCO->at(iter)                               = static_cast<proshade_double> ( iter ) *
                                                        ( ( static_cast<proshade_double> ( M_PI ) * 2.0 ) / static_cast<proshade_double> ( this->localAngRes ) ) -
                                                          ( static_cast<proshade_double> ( M_PI ) );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function fills in the vector of lattitudal bin boarder values.
 
    This is a simple abstraction function which finds the boarders of the angular grid binds along the lattitudal angle.
 
    \param[in] latCO A pointer for vector of proshade_double's to which the cutt-off values are to be saved. The vector needs to be initialised for this->localAngRes + 1 size.
 */
void ProSHADE_internal_spheres::ProSHADE_sphere::getLattitudeCutoffs ( std::vector<proshade_double>* latCO )
{
    //================================================ Get the cut-offs
    for ( proshade_unsign iter = 0; iter <= this->localAngRes; iter++ )
    {
        latCO->at(iter)                               = ( static_cast<proshade_double> ( iter ) *
                                                        ( static_cast<proshade_double> ( M_PI ) / static_cast<proshade_double> ( this->localAngRes ) ) -
                                                        ( static_cast<proshade_double> ( M_PI ) / 2.0 ) );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds the x, y and z positions of supplied shell point.
 
    This function converts the angular grid longitude and lattitude position (defined as in between the cutt-offs of the given index and its
    neighbour) to the X, Y and Z axis positions, which it saves in the supplied variables.
 
    \param[in] x A pointer to variable where the x position is to be saved.
    \param[in] y A pointer to variable where the y position is to be saved.
    \param[in] z A pointer to variable where the z position is to be saved.
    \param[in] thetaIt The longitude cut-off index of the grid point of which the XYZ are to be computed.
    \param[in] lonCO A pointer for vector of proshade_double's where the longitude cutt-off values are saved
    \param[in] phiIt The lattitude cut-off index of the grid point of which the XYZ are to be computed.
    \param[in] lonCO A pointer for vector of proshade_double's where the lattitude cutt-off values are saved
 */
void ProSHADE_internal_spheres::ProSHADE_sphere::getInterpolationXYZ ( proshade_double* x, proshade_double* y, proshade_double* z, proshade_double thetaIt, std::vector<proshade_double>* lonCO, proshade_unsign phiIt, std::vector<proshade_double>* latCO )
{
    //================================================ Compute and save XYZ interpolation positions
   *x                                                 = this->sphereRadius * std::cos ( ( lonCO->at( static_cast<proshade_unsign> ( thetaIt ) ) + lonCO->at( static_cast<proshade_unsign> ( thetaIt+1 ) ) ) / 2.0 ) *
                                                                             std::cos ( ( latCO->at( static_cast<proshade_unsign> ( phiIt ) )   + latCO->at( static_cast<proshade_unsign> ( phiIt+1 ) )   ) / 2.0 );
   *y                                                 = this->sphereRadius * std::sin ( ( lonCO->at( static_cast<proshade_unsign> ( thetaIt ) ) + lonCO->at( static_cast<proshade_unsign> ( thetaIt+1 ) ) ) / 2.0 ) *
                                                                             std::cos ( ( latCO->at( static_cast<proshade_unsign> ( phiIt ) )   + latCO->at( static_cast<proshade_unsign> ( phiIt+1 ) )   ) / 2.0 );
   *z                                                 = this->sphereRadius * std::sin ( ( latCO->at( static_cast<proshade_unsign> ( phiIt ) )   + latCO->at( static_cast<proshade_unsign> ( phiIt+1 ) )   ) / 2.0 );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function fills in the interpolation vector for a single map point.
 
    This function takes the XYZ position of a point and finds the top and bottom map points which surround the XYZ position.
 
    \param[in] xDimMax The internal map maximum index in the x dimension.
    \param[in] yDimMax The internal map maximum index in the y dimension.
    \param[in] zDimMax The internal map maximum index in the z dimension.
    \param[in] x The x position of the point for which the surrounding indices are to be found.
    \param[in] y The y position of the point for which the surrounding indices are to be found.
    \param[in] z The z position of the point for which the surrounding indices are to be found.
    \param[in] xBottom A pointer to the variable which should now hold the x position of the lower surrounding point.
    \param[in] yBottom A pointer to the variable which should now hold the y position of the lower surrounding point.
    \param[in] zBottom A pointer to the variable which should now hold the z position of the lower surrounding point.
    \param[in] xTop A pointer to the variable which should now hold the x position of the upper surrounding point.
    \param[in] yTop A pointer to the variable which should now hold the y position of the upper surrounding point.
    \param[in] zTop A pointer to the variable which should now hold the z position of the upper surrounding point.
 */
void ProSHADE_internal_spheres::ProSHADE_sphere::getXYZTopBottoms ( proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax, proshade_double x, proshade_double y, proshade_double z, proshade_signed* xBottom, proshade_signed* yBottom, proshade_signed* zBottom, proshade_signed* xTop, proshade_signed* yTop, proshade_signed* zTop )
{
    //================================================ Get the values
   *xBottom                                           = static_cast< proshade_signed > ( std::floor ( ( x / static_cast< proshade_double > ( this->xDimSampling ) ) ) + ( static_cast< proshade_double > ( xDimMax ) / 2.0 ) );
   *yBottom                                           = static_cast< proshade_signed > ( std::floor ( ( y / static_cast< proshade_double > ( this->yDimSampling ) ) ) + ( static_cast< proshade_double > ( yDimMax ) / 2.0 ) );
   *zBottom                                           = static_cast< proshade_signed > ( std::floor ( ( z / static_cast< proshade_double > ( this->zDimSampling ) ) ) + ( static_cast< proshade_double > ( zDimMax ) / 2.0 ) );
            
   *xTop                                              = *xBottom + 1;
   *yTop                                              = *yBottom + 1;
   *zTop                                              = *zBottom + 1;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function returns the local bandwidth.
 
    This is a simple accessor function so that the spherical harmonics can allocate the proper memory given the local
    bandwidth of the sphere.
 
    \param[out] localBandwidth The value of the local bandwidth for this particular shell.
 */
proshade_unsign ProSHADE_internal_spheres::ProSHADE_sphere::getLocalBandwidth ( void )
{
    //================================================ Return the value
    return                                            ( this->localBandwidth );
}

/*! \brief This function returns the local angular resolution.
 
    This is a simple accessor function so that the inverse sphere to Cartesian interpolation knows the proper angular
    cut-offs for each sphere.
 
    \param[out] localAngRes The value of the local angular resolution for this particular shell.
 */
proshade_unsign ProSHADE_internal_spheres::ProSHADE_sphere::getLocalAngRes ( void )
{
    //================================================ Return the value
    return                                            ( this->localAngRes );
}

/*! \brief This function returns the mapped data  array.
 
    This is a simple accessor function so that the spherical harmonics can access the mapped data in an
    array.
 
    \param[out] mappedData proshade_double pointer to an array of mapped data.
 */
proshade_double* ProSHADE_internal_spheres::ProSHADE_sphere::getMappedData ( void )
{
    //================================================ Return the value
    return                                            ( this->mappedData );
}

/*! \brief This function interpolates along the first (X) dimension.
 
    This function takes the 8 interpolation vectors for the neighbouring map points and proceeds to interpolate out the X axis
    from all of them, thus resulting in 4 new interpolation vectors which contain the Y and Z dimension values. This is the first
    step of the tri-linear interpolation procedure.
 
    \param[in] c000 The interpolation vector for map point less less less (than the point of interest).
    \param[in] c001 The interpolation vector for map point less less more (than the point of interest).
    \param[in] c010 The interpolation vector for map point less more less (than the point of interest).
    \param[in] c011 The interpolation vector for map point less more more (than the point of interest).
    \param[in] c100 The interpolation vector for map point more less less (than the point of interest).
    \param[in] c101 The interpolation vector for map point more less more (than the point of interest).
    \param[in] c110 The interpolation vector for map point more more less (than the point of interest).
    \param[in] c111 The interpolation vector for map point more more more (than the point of interest).
    \param[in] c00 Pointer to the variable where the less less second interporetation vector will be saved.
    \param[in] c01 Pointer to the variable where the less more second interporetation vector will be saved.
    \param[in] c10 Pointer to the variable where the more less second interporetation vector will be saved.
    \param[in] c11 Pointer to the variable where the more more second interporetation vector will be saved.
    \param[in] xd The relative distance of the point of interest to the lower map point along the x axis.
 */
void ProSHADE_internal_spheres::ProSHADE_sphere::interpolateAlongFirst ( std::vector<proshade_double> c000, std::vector<proshade_double> c001, std::vector<proshade_double> c010, std::vector<proshade_double> c011, std::vector<proshade_double> c100, std::vector<proshade_double> c101, std::vector<proshade_double> c110, std::vector<proshade_double> c111, std::vector<proshade_double>* c00, std::vector<proshade_double>* c01, std::vector<proshade_double>* c10, std::vector<proshade_double>* c11, proshade_double xd )
{
    //================================================ Interpolate for the less less point
    c00->at(0)                                        = ( static_cast< proshade_double > ( this->xDimSampling ) * xd ) + c000.at(0);
    c00->at(1)                                        = c000.at(1);
    c00->at(2)                                        = c000.at(2);
    c00->at(3)                                        = ( c000.at(3) * ( 1.0 - xd ) ) + ( c100.at(3) * xd );
    
    //================================================ Interpolate for the less more point
    c01->at(0)                                        = ( static_cast< proshade_double > ( this->xDimSampling ) * xd ) + c001.at(0);
    c01->at(1)                                        = c001.at(1);
    c01->at(2)                                        = c001.at(2);
    c01->at(3)                                        = ( c001.at(3) * ( 1.0 - xd ) ) + ( c101.at(3) * xd );
    
    //================================================ Interpolate for the more less point
    c10->at(0)                                        = ( static_cast< proshade_double > ( this->xDimSampling ) * xd ) + c010.at(0);
    c10->at(1)                                        = c010.at(1);
    c10->at(2)                                        = c010.at(2);
    c10->at(3)                                        = ( c010.at(3) * ( 1.0 - xd ) ) + ( c110.at(3) * xd );
    
    //================================================ Interpolate for the more more point
    c11->at(0)                                        = ( static_cast< proshade_double > ( this->xDimSampling ) * xd ) + c011.at(0);
    c11->at(1)                                        = c011.at(1);
    c11->at(2)                                        = c011.at(2);
    c11->at(3)                                        = ( c011.at(3) * ( 1.0 - xd ) ) + ( c111.at(3) * xd );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function interpolates along the second (Y) dimension.
 
    This function takes the 4 interpolation vectors for the neighbouring map points with the X axis already interpolate (that is, as
    outputted from the interpolateAlongFirst() function) and proceeds to interpolate out the Y axis from all of them, thus resulting
    in 2 new interpolation vectors which contain the Z dimension values only. This is the second step of the tri-linear interpolation
    procedure.
 
    \param[in] c00 The interpolation vector for map point less less (than the point of interest) as returned from interpolateAlongFirst().
    \param[in] c01 The interpolation vector for map point less more (than the point of interest) as returned from interpolateAlongFirst().
    \param[in] c10 The interpolation vector for map point more less (than the point of interest) as returned from interpolateAlongFirst().
    \param[in] c11 The interpolation vector for map point more more (than the point of interest) as returned from interpolateAlongFirst().
    \param[in] c0 Pointer to the variable where the less third interporetation vector will be saved.
    \param[in] c1 Pointer to the variable where the more third interporetation vector will be saved.
    \param[in] yd The relative distance of the point of interest to the lower map point along the y axis.
 */
void ProSHADE_internal_spheres::ProSHADE_sphere::interpolateAlongSecond ( std::vector<proshade_double> c00, std::vector<proshade_double> c01, std::vector<proshade_double> c10, std::vector<proshade_double> c11, std::vector<proshade_double>* c0, std::vector<proshade_double>* c1, proshade_double yd )
{
    //================================================ Interpolate for the less point
    c0->at(0)                                         = c00.at(0);
    c0->at(1)                                         = ( static_cast< proshade_double > ( this->yDimSampling ) * yd ) + c00.at(1);
    c0->at(2)                                         = c00.at(2);
    c0->at(3)                                         = ( c00.at(3) * ( 1.0 - yd ) ) + ( c10.at(3) * yd );
    
    //================================================ Interpolate for the more point
    c1->at(0)                                         = c01.at(0);
    c1->at(1)                                         = ( static_cast< proshade_double > ( this->yDimSampling ) * yd ) + c01.at(1);
    c1->at(2)                                         = c01.at(2);
    c1->at(3)                                         = ( c01.at(3) * ( 1.0 - yd ) ) + ( c11.at(3) * yd );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function determines the bandwidth for the spherical harmonics computation.
 
    This function is here to automstically determine the bandwidth to which the spherical harmonics computations should be done.
    It accomplishes this by setting it to half of the maximum circumference of the map/sphere, in indices as recommended by
    Kostelec and Rockmore (2007).
 
    \param[in] circumference The maximum circumference of the map/sphere.
 */
proshade_unsign ProSHADE_internal_spheres::autoDetermineBandwidth ( proshade_unsign circumference )
{
    //================================================ Determine and return
    if ( static_cast<proshade_unsign> ( std::ceil ( circumference / 2 ) ) % 2 == 0 )
    {
        return                                        ( static_cast<proshade_unsign> ( std::ceil ( circumference / 2 ) ) );
    }
    else
    {
        return                                        ( static_cast<proshade_unsign> ( std::ceil ( circumference / 2 ) ) + 1 );
    }
}

/*! \brief This function determines the sphere distances for sphere mapping.
 
    This function determines the distance between two consecutive spheres in the sphere mappin galgorithm. It sets it
    as the sampling rate (distance between any two map points). It then checks that there will be at least 10 spheres
    and if not, it changes the sphere distance until at least 10 spheres are to be produced.
 
    \param[in] maxMapRange The maximum diagonal distance of the map in Angstroms.
    \param[in] resolution The resolution to which the computations are to be done.
 */
proshade_single ProSHADE_internal_spheres::autoDetermineSphereDistances ( proshade_single maxMapRange, proshade_single resolution )
{
    //================================================ Get starting point
    proshade_single ret                               = resolution / 2.0f;
    
    //================================================ Make sure at least 10 shells will exist
    while ( std::floor ( maxMapRange / ret ) < 10 )
    {
        ret                                          /= 2.0f;
    }
    
    //================================================ Done
    return                                            ( ret );
}

/*! \brief This function determines the integration order for the between spheres integration.
 
    This function determines the order of the Gauss-Legendre integration which needs to be done between the spheres. To do
    this, it uses the pre-coputed values of maxium distance between integration points for each order and the maxium distance
    between spheres expressed as a fraction of the total.
 
    \param[in] maxMapRange The maximum diagonal distance of the map in Angstroms.
    \param[in] sphereDist The distance between spheres.
 */
proshade_unsign ProSHADE_internal_spheres::autoDetermineIntegrationOrder ( proshade_single maxMapRange, proshade_single sphereDist )
{
    //================================================ Initialise local variables
    proshade_double sphereDistanceAsFractionOfTotal = static_cast<proshade_double> ( sphereDist / ( maxMapRange / 2.0f ) );
    proshade_unsign ret                               = 0;
    
    //================================================ Compare to precomputed values
    for ( proshade_unsign iter = 2; iter < static_cast<proshade_unsign> ( 10000 ); iter++ )
    {
        if ( ProSHADE_internal_precomputedVals::glIntMaxDists[iter] >= sphereDistanceAsFractionOfTotal )
        {
            ret                                       = iter;
        }
    }
    
    //================================================ Return largest passing value
    return                                            ( ret );
    
}

/*! \brief This function returns the radius of the shell in question.
 
    This is a simple accessor function so that the sphere radius can be accessed when need be.
 
    \param[out] sphereRadius The distance of the shell to the centre of the coordinates.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_sphere::getShellRadius ( void )
{
    //================================================ Done
    return                                            ( this->sphereRadius );
    
}

/*! \brief This function allocates the rotated map memory.
 */
void ProSHADE_internal_spheres::ProSHADE_sphere::allocateRotatedMap ( void )
{
    //================================================ Allocate memory for sphere mapping
    this->mappedDataRot                               = nullptr;
    this->mappedDataRot                               = new proshade_double[this->localAngRes * this->localAngRes];
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->mappedDataRot, __FILE__, __LINE__, __func__ );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function sets the rotated mapped data value to the given position.
 
    \param[in] pos The array position to which the value will be set to.
    \param[in] value The value which should be set to the given position.
 */
void ProSHADE_internal_spheres::ProSHADE_sphere::setRotatedMappedData ( proshade_unsign pos, proshade_double value )
{
    //================================================ Set the value to the position
    this->mappedDataRot[pos]                          = value;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function gets the rotated mapped data value for a particular position.
 
    \param[in] pos The array position from which the value should be retrieved.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_sphere::getRotatedMappedData ( proshade_unsign pos )
{
    //================================================ Done
    return                                            ( this->mappedDataRot[pos]  );
    
}

/*! \brief Constructor for getting empty ProSHADE_rotFun_sphere class.
 
    This function simply creates an object of the ProSHADE_rotFun_sphere class and allocates the memory as required.
 
    \param[in] rad The radius of this sphere.
    \param[in] radRange The range in the radial dimension covered by this sphere.
    \param[in] dim The required size of the sampling grid.
    \param[in] repAng The angle represented by this sphere.
    \param[in] sphNo The number of this sphere in the spheres vector.
    \param[out] X Data object with all values set and ready to accept the rotation function values.
 */
ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::ProSHADE_rotFun_sphere ( proshade_double rad, proshade_double radRange, proshade_unsign dim, proshade_double repAng, proshade_unsign sphNo )
{
    //================================================ Set internal values
    this->radius                                      = rad;
    this->angularDim                                  = dim;
    this->radiusMin                                   = this->radius - ( radRange / 2.0 );
    this->radiusMax                                   = this->radius + ( radRange / 2.0 );
    this->representedAngle                            = repAng;
    this->sphereNumber                                = sphNo;
    
    //================================================ Allocate the axis field
    this->axesValues                                  = new proshade_double[dim*dim];
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->axesValues, __FILE__, __LINE__, __func__ );
    
    //================================================ Fill axis field with zeroes
    for ( proshade_unsign iter = 0; iter < ( dim * dim ); iter++ ) { this->axesValues[iter] = 0.0; }
    
}

/*! \brief Destructor for releasing memory from the ProSHADE_rotFun_sphere class.
 
    \param[out] X N/A.
 */
ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::~ProSHADE_rotFun_sphere (  )
{
    //================================================ Release the allocated memory
    if ( this->axesValues != nullptr )
    {
        delete[] this->axesValues;
    }
}

/*! \brief Accessor function for the private variable radius.
 
    \param[out] radius The value of the radius of this specific concentric shell.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::getRadius ( void )
{
    //================================================ Done
    return                                            ( this->radius );
}

/*! \brief Accessor function for the private variable maximum radius.
 
    \param[out] radius The value of the maximum radius of this specific concentric shell.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::getMaxRadius ( void )
{
    //================================================ Done
    return                                            ( this->radiusMax );
    
}

/*! \brief Accessor function for the private variable angular dim.
 
    \param[out] radius The dimension size of the angular sampling grid.
 */
proshade_unsign ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::getAngularDim ( void )
{
    //================================================ Done
    return                                            ( this->angularDim );
    
}

/*! \brief Accessor function for the private variable minimal radius.
 
    \param[out] radius The value of the minimal radius of this specific concentric shell.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::getMinRadius ( void )
{
    //================================================ Done
    return                                            ( this->radiusMin );
    
}

/*! \brief Accessor function for the private variable represented angle.
 
    \param[out] radius The value of the angle represented by this sphere.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::getRepresentedAngle ( void )
{
    //================================================ Done
    return                                            ( this->representedAngle );
    
}

/*! \brief Accessor function for the private variable sphere number.
 
    \param[out] sphereNumber The sphere number.
 */
proshade_unsign ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::getSphereNumber ( void )
{
    //================================================ Done
    return                                            ( this->sphereNumber );
    
}

/*! \brief Accessor function for the private variable containing all detected peaks.
 
    \param[out] peaks The vector containing all detected peaks.
 */
std::vector<std::pair<proshade_unsign,proshade_unsign>> ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::getPeaks ( void )
{
    //================================================ Done
    return                                            ( this->peaks );
    
}

/*! \brief Function for interpolating the sphere grid values from angle-axis converted rotation function.
 
    This function starts by converting each of the sphere sampling points lattitude and longitude to the XYZ position and by adding the represented
    rotation angle, obtains the angle-axis representation for the given point. It then proceeds to locate such points exact position in the indices space
    of the supplied rotation map. From there, it interpolates the exact correlation value for the given point, thus effectivelly re-sampling the rotation
    function space onto the sphere.
 
    \param[in] rotFun proshade_complex pointer to the rotation function values.
 */
void ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::interpolateSphereValues ( proshade_complex* rotFun )
{
    //================================================ Initialise variables
    proshade_double lonSampling                       = ( M_PI       ) / static_cast< proshade_double > ( this->angularDim );
    proshade_double latSampling                       = ( M_PI * 2.0 ) / static_cast< proshade_double > ( this->angularDim );
    
    proshade_double lat, lon, cX, cY, cZ, c000, c001, c010, c011, c100, c101, c110, c111, c00, c01, c10, c11, c0, c1, xRelative, yRelative, zRelative, eulerAlpha, eulerBeta, eulerGamma, mapX, mapY, mapZ;
    proshade_signed xBottom, xTop, yBottom, yTop, zBottom, zTop, mapIndex;
    
    //================================================ For each sphere grid position
    for ( proshade_signed lonIt = 0; lonIt < static_cast<proshade_signed> ( this->angularDim ); lonIt++ )
    {
        for ( proshade_signed latIt = 0; latIt < static_cast<proshade_signed> ( this->angularDim ); latIt++ )
        {   
            //======================================== Convert to XYZ position on unit sphere. The radius here is not important, as it does not change the direction of the vector.
            lon                                       = static_cast<proshade_double> ( lonIt ) * lonSampling;
            lat                                       = static_cast<proshade_double> ( latIt ) * latSampling;
            cX                                        = 1.0 * std::sin ( lon ) * std::cos ( lat );
            cY                                        = 1.0 * std::sin ( lon ) * std::sin ( lat );
            cZ                                        = 1.0 * std::cos ( lon );
            
            //======================================== Convert to ZXZ Euler angles
            ProSHADE_internal_maths::getEulerZXZFromAngleAxis ( cX, cY, cZ, this->representedAngle, &eulerAlpha, &eulerBeta, &eulerGamma );
            
            //======================================== Convert to SOFT map position (decimal, not indices)
            ProSHADE_internal_maths::getSOFTPositionFromEulerZXZ ( this->angularDim / 2, eulerAlpha, eulerBeta, eulerGamma, &mapX, &mapY, &mapZ );
            
            //======================================== Find lower and higher points and deal with boundaries
            xBottom = static_cast< proshade_signed > ( std::floor ( mapX ) ); if ( xBottom < 0 ) { xBottom += this->angularDim; } if ( xBottom >= static_cast<proshade_signed> ( this->angularDim ) ) { xBottom -= static_cast<proshade_signed> ( this->angularDim ); }
            yBottom = static_cast< proshade_signed > ( std::floor ( mapY ) ); if ( yBottom < 0 ) { yBottom += this->angularDim; } if ( yBottom >= static_cast<proshade_signed> ( this->angularDim ) ) { yBottom -= static_cast<proshade_signed> ( this->angularDim ); }
            zBottom = static_cast< proshade_signed > ( std::floor ( mapZ ) ); if ( zBottom < 0 ) { zBottom += this->angularDim; } if ( zBottom >= static_cast<proshade_signed> ( this->angularDim ) ) { zBottom -= static_cast<proshade_signed> ( this->angularDim ); }
            xTop = static_cast< proshade_signed > ( std::ceil ( mapX ) ); if ( xTop < 0 ) { xTop += this->angularDim; } if ( xTop >= static_cast<proshade_signed> ( this->angularDim ) ) { xTop -= static_cast<proshade_signed> ( this->angularDim ); }
            yTop = static_cast< proshade_signed > ( std::ceil ( mapY ) ); if ( yTop < 0 ) { yTop += this->angularDim; } if ( yTop >= static_cast<proshade_signed> ( this->angularDim ) ) { yTop -= static_cast<proshade_signed> ( this->angularDim ); }
            zTop = static_cast< proshade_signed > ( std::ceil ( mapZ ) ); if ( zTop < 0 ) { zTop += this->angularDim; } if ( zTop >= static_cast<proshade_signed> ( this->angularDim ) ) { zTop -= static_cast<proshade_signed> ( this->angularDim ); }
            
            //======================================== Start X interpolation - bottom, bottom, bottom
            mapIndex                                  = zBottom + static_cast< proshade_signed > ( this->angularDim ) * ( yBottom + static_cast< proshade_signed > ( this->angularDim ) * xBottom );
            c000                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - bottom, bottom, top
            mapIndex                                  = zTop    + static_cast< proshade_signed > ( this->angularDim ) * ( yBottom + static_cast< proshade_signed > ( this->angularDim ) * xBottom );
            c001                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - bottom, top, bottom
            mapIndex                                  = zBottom + static_cast< proshade_signed > ( this->angularDim ) * ( yTop    + static_cast< proshade_signed > ( this->angularDim ) * xBottom );
            c010                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - bottom, top, top
            mapIndex                                  = zTop    + static_cast< proshade_signed > ( this->angularDim ) * ( yTop    + static_cast< proshade_signed > ( this->angularDim ) * xBottom );
            c011                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - top, bottom, bottom
            mapIndex                                  = zBottom + static_cast< proshade_signed > ( this->angularDim ) * ( yBottom + static_cast< proshade_signed > ( this->angularDim ) * xTop    );
            c100                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - top, bottom, top
            mapIndex                                  = zTop    + static_cast< proshade_signed > ( this->angularDim ) * ( yBottom + static_cast< proshade_signed > ( this->angularDim ) * xTop    );
            c101                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - top, top, bottom
            mapIndex                                  = zBottom + static_cast< proshade_signed > ( this->angularDim ) * ( yTop    + static_cast< proshade_signed > ( this->angularDim ) * xTop    );
            c110                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - top, top, top
            mapIndex                                  = zTop    + static_cast< proshade_signed > ( this->angularDim ) * ( yTop    + static_cast< proshade_signed > ( this->angularDim ) * xTop    );
            c111                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== Solve for X
            xRelative                                 = mapX - std::floor( mapX );
            c00                                       = ( c000 * ( 1.0 - xRelative ) ) + ( c100 * xRelative );
            c01                                       = ( c001 * ( 1.0 - xRelative ) ) + ( c101 * xRelative );
            c10                                       = ( c010 * ( 1.0 - xRelative ) ) + ( c110 * xRelative );
            c11                                       = ( c011 * ( 1.0 - xRelative ) ) + ( c111 * xRelative );
            
            //======================================== Solve for Y
            yRelative                                 = mapY - std::floor( mapY );
            c0                                        = ( c00 * ( 1.0 - yRelative ) ) + ( c10 * yRelative );
            c1                                        = ( c01 * ( 1.0 - yRelative ) ) + ( c11 * yRelative );
            
            //======================================== Solve for Z
            zRelative                                 = mapZ - std::floor( mapZ );
            
            //======================================== Save result
            mapIndex                                  = lonIt + ( latIt * static_cast< proshade_signed > ( this->angularDim ) );
            this->axesValues[mapIndex]                = ( c0 * ( 1.0 - zRelative ) ) + ( c1 * zRelative );
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Accessor function for specific lattitude and longitude point of the sphere sampling grid.
 
    \param[in] lattitude The lattitude index of the requested sampling grid point.
    \param[in] longitude The longitude index of the requested sampling grid point.
    \param[out] radius The value of the sampling grid point at given lattitude and longitude position.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::getSphereLatLonPosition ( proshade_unsign lattitude, proshade_unsign longitude )
{
    //================================================ Done
    return                                            ( this->axesValues[longitude + ( lattitude * static_cast<proshade_unsign> ( this->angularDim ) )] );
}

/*! \brief Function for obtaining sphere values outside of the grid points.
 
    This function computes linear interpolation value for any grid position of the sampling grid of this sphere.
 
    \param[in] lattitude The lattitude index decimal place of the requested sampling grid position.
    \param[in] longitude The longitude index of the requested sampling grid position.
    \param[out] radius The value of the sampling grid position at given lattitude and longitude.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::getSphereLatLonLinearInterpolationPos ( proshade_double lattitude, proshade_double longitude )
{
    //================================================ Initialise variables
    proshade_double c00, c01, c10, c11, c0, c1, latRelative, lonRelative;
    proshade_signed latTop, latBottom, lonTop, lonBottom, gridIndex;
    
    //================================================ Find lower and higher indices and deal with boundaries
    latBottom = static_cast< proshade_signed > ( std::floor ( lattitude ) ); if ( latBottom < 0 ) { latBottom += this->angularDim; } if ( latBottom >= static_cast<proshade_signed> ( this->angularDim ) ) { latBottom -= this->angularDim; }
    lonBottom = static_cast< proshade_signed > ( std::floor ( longitude ) ); if ( lonBottom < 0 ) { lonBottom += this->angularDim; } if ( lonBottom >= static_cast<proshade_signed> ( this->angularDim ) ) { lonBottom -= this->angularDim; }
    latTop    = static_cast< proshade_signed > ( std::ceil  ( lattitude ) ); if ( latTop    < 0 ) { latTop    += this->angularDim; } if ( latTop    >= static_cast<proshade_signed> ( this->angularDim ) ) { latTop    -= this->angularDim; }
    lonTop    = static_cast< proshade_signed > ( std::ceil  ( longitude ) ); if ( lonTop    < 0 ) { lonTop    += this->angularDim; } if ( lonTop    >= static_cast<proshade_signed> ( this->angularDim ) ) { lonTop    -= this->angularDim; }
    
    //================================================ Interpolate
    gridIndex                                         = lonBottom + ( latBottom * static_cast< proshade_signed > ( this->angularDim ) );
    c00                                               = this->axesValues[gridIndex];
    
    gridIndex                                         = lonBottom + ( latTop    * static_cast< proshade_signed > ( this->angularDim ) );
    c01                                               = this->axesValues[gridIndex];
    
    gridIndex                                         = lonTop    + ( latBottom * static_cast< proshade_signed > ( this->angularDim ) );
    c10                                               = this->axesValues[gridIndex];
    
    gridIndex                                         = lonTop    + ( latTop    * static_cast< proshade_signed > ( this->angularDim ) );
    c11                                               = this->axesValues[gridIndex];
    
    //================================================ Solve for longitude
    lonRelative                                       = longitude - std::floor( longitude );
    c0                                                = ( c00 * ( 1.0 - lonRelative ) ) + ( c10 * lonRelative );
    c1                                                = ( c01 * ( 1.0 - lonRelative ) ) + ( c11 * lonRelative );
    
    //================================================ Solve for lattitude
    latRelative                                       = lattitude - std::floor ( lattitude );
    proshade_double res                               = ( c0 * ( 1.0 - latRelative ) ) + ( c1 * latRelative );
    
    //================================================ Done
    return                                            ( res );
    
}

/*! \brief Function for obtaining a copy of all sphere values.
 
    This function gives access to the mapped sphere values.
 
    \param[out] ret A vector of vectors of doubles indexed longitude, latitude containing all sphere mapped  RF values.
 */
std::vector< std::vector< proshade_double > > ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::getCopyOfValues ( )
{
    //================================================ Initialise variables
    std::vector< std::vector< proshade_double > > ret ( this->angularDim );
    std::vector< proshade_double > retHlp             ( this->angularDim, 0.0 );
    
    //================================================ Fill in the values
    for ( size_t lonIt = 0; lonIt < this->angularDim; lonIt++ )
    {
        for ( size_t latIt = 0; latIt < this->angularDim; latIt++ )
        {
            retHlp.at(latIt)                          = this->getSphereLatLonPosition ( static_cast< proshade_unsign > ( latIt ), static_cast< proshade_unsign > ( lonIt ) );
        }
        
        ret.at(lonIt)                                 = retHlp;
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief Function for finding all peaks in the sampling grid.
 
    This function takes the values on the sampling grid and does a naive peak search, saving the peak position into an internal variable.
 
    \param[in] noSmNeighbours The number of surrounding peaks in any direction that need to be smaller for a value to be a peak.
    \param[in] allHeights A vector to which all detected non-peaks heights will be saved into. This will later be used to determine the threshold for "small" peaks.
 */
void ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::findAllPeaks ( proshade_signed noSmNeighbours, std::vector< proshade_double >* allHeights )
{
    //================================================ Initialise local variables
    proshade_double currentHeight;
    proshade_signed nbLat, nbLon;
    bool isPeak;
    
    //================================================ Find all peaks
    for ( proshade_signed latIt = 0; latIt < static_cast<proshade_signed> ( this->angularDim ); latIt++ )
    {
        for ( proshade_signed lonIt = 0; lonIt < static_cast<proshade_signed> ( this->angularDim ); lonIt++ )
        {
            //======================================== Initialise peak search
            currentHeight                             = this->getSphereLatLonPosition ( static_cast< proshade_unsign > ( latIt ), static_cast< proshade_unsign > ( lonIt ) );
            isPeak                                    = true;
            
            //======================================== Find all neighbours in the same sphere
            for ( proshade_signed latRound = -noSmNeighbours; latRound <= noSmNeighbours; latRound++ )
            {
                for ( proshade_signed lonRound = -noSmNeighbours; lonRound <= noSmNeighbours; lonRound++ )
                {
                    //================================ Ignore same point
                    if ( latRound == 0 && lonRound == 0 ) { continue; }
                    
                    //================================ Get neighbour height
                    nbLat                             = latIt + latRound;
                    nbLon                             = lonIt + lonRound;
                    if ( nbLat < 0 ) { nbLat += this->angularDim; } if ( nbLat >= static_cast<proshade_signed> ( this->angularDim ) ) { nbLat -= this->angularDim; }
                    if ( nbLon < 0 ) { nbLon += this->angularDim; } if ( nbLon >= static_cast<proshade_signed> ( this->angularDim ) ) { nbLon -= this->angularDim; }
                    
                    //================================ If this value is larger than the tested one, no peak
                    if ( this->getSphereLatLonPosition ( static_cast< proshade_unsign > ( nbLat ), static_cast< proshade_unsign > ( nbLon ) ) > currentHeight ) { isPeak = false; break; }
                }
                
                if ( !isPeak ) { break; }
            }
            
            if ( isPeak )
            {
                //==================================== Save!
                this->peaks.emplace_back              ( std::pair<proshade_unsign,proshade_unsign> ( latIt, lonIt ) );
            }
            else
            {
                ProSHADE_internal_misc::addToDoubleVector ( allHeights, currentHeight );
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for removing peaks with too small height.
 
    This function takes the threshold for peaks being too small (computed from all peaks as returned by findAllPeaks() function, but could be any other value) and removes
    all peaks with height below this threshold.
 
    \param[in] peakThres The height above which a peak needs to be in order not to be deleted.
    \param[in] minThres The minimum threshold that needs to be passed for a peak to be believed in.
 */
void ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::removeSmallPeaks( proshade_double peakThres )
{
    //================================================ Initialise variables
    proshade_double curHeight;
    std::vector< proshade_unsign > dels ( 0, static_cast< proshade_unsign > ( this->peaks.size() ) );
    
    //================================================ For each peak in this sphere
    for ( proshade_unsign peakIt = 0; peakIt < static_cast<proshade_unsign> ( this->peaks.size() ); peakIt++ )
    {
        //============================================ Find the peak height
        curHeight                                     = this->getSphereLatLonPosition ( this->peaks.at(peakIt).first, this->peaks.at(peakIt).second );
        
        //============================================ Should this peak be deleted?
        if ( curHeight < peakThres )
        {
            ProSHADE_internal_misc::addToUnsignVector ( &dels, peakIt );
        }
    }
    
    //================================================ Descending sort to avoid changing the order
    std::sort                                         ( dels.begin(), dels.end(), std::greater <proshade_unsign>() );
    
    //================================================ Delete the low peaks
    for ( proshade_unsign delIt = 0; delIt < static_cast< proshade_unsign > ( dels.size() ); delIt++ )
    {
        this->peaks.erase                             ( this->peaks.begin() + static_cast< long int > ( dels.at(delIt) ) );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Constructor for getting empty ProSHADE_rotFun_spherePeakGroup class.
 
    This function simply creates an object of the ProSHADE_rotFun_spherePeakGroup class and fills in the initial data.
 
    \param[in] lat The lattitude value of the first peak of the group.
    \param[in] lon The longitude value of the first peak of the group.
    \param[in] sphPos The sphere number of the peak.
    \param[in] angDim The dimensionality of the sphere grid that we are processing.
    \param[out] X Data object with all values set and ready to add new data or search for point groups in the supplied peaks.
 */
ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::ProSHADE_rotFun_spherePeakGroup ( proshade_double lat, proshade_double lon, proshade_unsign sphPos, proshade_unsign angDim )
{
    //================================================ Compute the run-specific values
    this->dimension                                   = angDim;
    this->lonSampling                                 = ( M_PI       ) / static_cast<proshade_double> ( this->dimension );
    this->latSampling                                 = ( M_PI * 2.0 ) / static_cast<proshade_double> ( this->dimension );
    
    //================================================ The constructor is called when firstt peak of the group is found. Save the values of this initial peak.
    this->latFrom                                     = static_cast<proshade_double> ( lat ) * this->latSampling;
    this->latTo                                       = static_cast<proshade_double> ( lat ) * this->latSampling;
    this->lonFrom                                     = static_cast<proshade_double> ( lon ) * this->lonSampling;
    this->lonTo                                       = static_cast<proshade_double> ( lon ) * this->lonSampling;
    this->latFromInds                                 = lat;
    this->latToInds                                   = lat;
    this->lonFromInds                                 = lon;
    this->lonToInds                                   = lon;
    ProSHADE_internal_misc::addToUnsignVector         ( &this->spherePositions, sphPos );
    
    //================================================ Allocate memory for similarity positions
    this->latMinLonMinXYZ                             = new proshade_double[3];
    this->latMaxLonMinXYZ                             = new proshade_double[3];
    this->latMinLonMaxXYZ                             = new proshade_double[3];
    this->latMaxLonMaxXYZ                             = new proshade_double[3];
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->latMinLonMinXYZ, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->latMaxLonMinXYZ, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->latMinLonMaxXYZ, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->latMaxLonMaxXYZ, __FILE__, __LINE__, __func__ );
    
    //================================================ Compute corner vectors
    this->computeCornerPositions                      ( );
    
    //================================================ Done
    
}

/*! \brief Destructor for the ProSHADE_rotFun_spherePeakGroup class.
 
    This function releases all memory allocated by the ProSHADE_rotFun_spherePeakGroup object.
 */
ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::~ProSHADE_rotFun_spherePeakGroup ( )
{
    //================================================ Release the XYZ arrays
    delete[] this->latMinLonMinXYZ;
    delete[] this->latMaxLonMinXYZ;
    delete[] this->latMinLonMaxXYZ;
    delete[] this->latMaxLonMaxXYZ;
}

/*! \brief This function computes the group corner vectors, saving results into internal variables.
 */
void ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::computeCornerPositions ( )
{
    //================================================ Compute corner vectors
    this->latMinLonMinXYZ[0]                          = 1.0 * std::sin ( this->lonFrom ) * std::cos ( this->latFrom );
    this->latMinLonMinXYZ[1]                          = 1.0 * std::sin ( this->lonFrom ) * std::sin ( this->latFrom );
    this->latMinLonMinXYZ[2]                          = 1.0 * std::cos ( this->lonFrom );
    
    this->latMaxLonMinXYZ[0]                          = 1.0 * std::sin ( this->lonFrom ) * std::cos ( this->latTo   );
    this->latMaxLonMinXYZ[1]                          = 1.0 * std::sin ( this->lonFrom ) * std::sin ( this->latTo   );
    this->latMaxLonMinXYZ[2]                          = 1.0 * std::cos ( this->lonFrom );
    
    this->latMinLonMaxXYZ[0]                          = 1.0 * std::sin ( this->lonTo )   * std::cos ( this->latFrom );
    this->latMinLonMaxXYZ[1]                          = 1.0 * std::sin ( this->lonTo )   * std::sin ( this->latFrom );
    this->latMinLonMaxXYZ[2]                          = 1.0 * std::cos ( this->lonTo );
    
    this->latMaxLonMaxXYZ[0]                          = 1.0 * std::sin ( this->lonTo )   * std::cos ( this->latTo );
    this->latMaxLonMaxXYZ[1]                          = 1.0 * std::sin ( this->lonTo )   * std::sin ( this->latTo );
    this->latMaxLonMaxXYZ[2]                          = 1.0 * std::cos ( this->lonTo );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes two lattitude or longitude positions and finds the smallest distance between them considering the border periodicity.
 
    \param[in] newAngul The first lattitude or longitude value.
    \param[in] currentAngul The second lattitude or longitude value.
    \param[out] ret The smallest distance between the first and the second lattitude or longitude values.
 */
proshade_signed ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::angularDistanceWithBorders ( proshade_signed newAngul, proshade_signed currentAngul )
{
    //================================================ Initialise variables
    proshade_signed smallerAngul                      = newAngul - static_cast< proshade_signed > ( this->dimension );
    proshade_signed largerAngul                       = newAngul + static_cast< proshade_signed > ( this->dimension );
    
    //================================================ Find the smallest distance
    proshade_signed ret                               = std::min ( std::abs ( currentAngul - newAngul ), std::min ( std::abs ( currentAngul - smallerAngul ), std::abs ( currentAngul - largerAngul ) ) );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function takes a new prospective peak and tests if it belongs to this peak group or not.
 
    This function takes a new peak position in terms of lattitude and longitude and proceeds to convert these to XYZ position. It then checks this
    XYZ position against this group's "corners" (i.e. the group's lattitude and longitude minimum and maximum borders). If the tested
    position belongs to the group (i.e. it has small cosine distance to one of the corners), then it is added and the group corners are updated. Otherwise,
    false is returned and nothing changes in the group.
 
    \param[in] lat The lattitude value of the first peak of the group.
    \param[in] lon The longitude value of the first peak of the group.
    \param[in] sphPos The sphere number of the peak.
    \param[in] cosTol The tolerance for cosine distance similarity to consider the two vectors similar.
    \param[in] verbose How verbose should the run be? Use -1 if you do not want any standard output output.
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
    \param[out] res Boolean value signifying if the peak was added.
 */
bool ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::checkIfPeakBelongs ( proshade_double lat, proshade_double lon, proshade_unsign sphPos, proshade_double cosTol, proshade_signed verbose, proshade_signed messageShift )
{
    //================================================ Initialise local variables
    bool peakAdded                                    = false;
    std::stringstream hlpSS;
    std::stringstream hlpSS2;
    
    //================================================ Compute peaks XYZ and its cosine distance to group corners
    proshade_double xPos                              = 1.0 * std::sin ( lon * this->lonSampling ) * std::cos ( lat * this->latSampling );
    proshade_double yPos                              = 1.0 * std::sin ( lon * this->lonSampling ) * std::sin ( lat * this->latSampling );
    proshade_double zPos                              = 1.0 * std::cos ( lon * this->lonSampling );
    hlpSS2 << "Peak " << xPos << " ; " << yPos << " ; " << zPos << " is close enough to group with corner ";
    
    if ( ProSHADE_internal_maths::vectorOrientationSimilaritySameDirection ( xPos, yPos, zPos, this->latMinLonMinXYZ[0], this->latMinLonMinXYZ[1], this->latMinLonMinXYZ[2], cosTol ) ) { peakAdded = true; hlpSS2 << this->latMinLonMinXYZ[0] << " ; " << this->latMinLonMinXYZ[1] << " ; " << this->latMinLonMinXYZ[2]; }
    if ( ProSHADE_internal_maths::vectorOrientationSimilaritySameDirection ( xPos, yPos, zPos, this->latMaxLonMinXYZ[0], this->latMaxLonMinXYZ[1], this->latMaxLonMinXYZ[2], cosTol ) && !peakAdded ) { peakAdded = true; hlpSS2 << this->latMaxLonMinXYZ[0] << " ; " << this->latMaxLonMinXYZ[1] << " ; " << this->latMaxLonMinXYZ[2]; }
    if ( ProSHADE_internal_maths::vectorOrientationSimilaritySameDirection ( xPos, yPos, zPos, this->latMinLonMaxXYZ[0], this->latMinLonMaxXYZ[1], this->latMinLonMaxXYZ[2], cosTol ) && !peakAdded ) { peakAdded = true; hlpSS2 << this->latMinLonMaxXYZ[0] << " ; " << this->latMinLonMaxXYZ[1] << " ; " << this->latMinLonMaxXYZ[2];  }
    if ( ProSHADE_internal_maths::vectorOrientationSimilaritySameDirection ( xPos, yPos, zPos, this->latMaxLonMaxXYZ[0], this->latMaxLonMaxXYZ[1], this->latMaxLonMaxXYZ[2], cosTol ) && !peakAdded ) { peakAdded = true; hlpSS2 << this->latMaxLonMaxXYZ[0] << " ; " << this->latMaxLonMaxXYZ[1] << " ; " << this->latMaxLonMaxXYZ[2];  }
    
    //================================================ If peak within corners, add it
    if ( peakAdded )
    {
        //============================================ Report progress
        hlpSS << "Peak group dimensions changed from LAT " << this->latFromInds << " - " << this->latToInds << " and LON " << this->lonFromInds << " - " << this->lonToInds << " to ";
        ProSHADE_internal_messages::printProgressMessage  ( verbose, 6, hlpSS2.str(), messageShift );
        
        //============================================ Initialise local variables
        proshade_signed largerCorner, smallerCorner;
        bool latCornersDone                           = false;
        bool lonCornersDone                           = false;
        
        //============================================ Check if lattitude boundaries need to be modified
        if ( ( this->latFromInds <= this->latToInds ) && !( ( lat >= this->latFromInds ) && ( lat <= this->latToInds ) ) )
        {
            //======================================== Lattitude is outside of group boundaries
            smallerCorner                             = angularDistanceWithBorders ( static_cast< proshade_signed > ( lat ), static_cast< proshade_signed > ( this->latFromInds ) );
            largerCorner                              = angularDistanceWithBorders ( static_cast< proshade_signed > ( lat ), static_cast< proshade_signed > ( this->latToInds   ) );
            
            if ( smallerCorner < largerCorner )
            {
                this->latFromInds                     = lat;
                latCornersDone                        = true;
            }
            if ( smallerCorner > largerCorner )
            {
                this->latToInds                       = lat;
                latCornersDone                        = true;
            }
            if ( smallerCorner == largerCorner )
            {
                if      ( lat < this->latFromInds )   { this->latFromInds = lat; latCornersDone = true; }
                else if ( lat > this->latToInds   )   { this->latToInds   = lat; latCornersDone = true; }
            }
        }
        
        if ( ( this->latFromInds >  this->latToInds ) && !( ( lat >= this->latFromInds ) || ( lat <= this->latToInds ) ) )
        {
            //======================================== Lattitude is outside of group boundaries
            smallerCorner                             = angularDistanceWithBorders ( static_cast< proshade_signed > ( lat ), static_cast< proshade_signed > ( this->latFromInds ) );
            largerCorner                              = angularDistanceWithBorders ( static_cast< proshade_signed > ( lat ), static_cast< proshade_signed > ( this->latToInds   ) );
            
            if ( smallerCorner < largerCorner )
            {
                this->latFromInds                     = lat;
                latCornersDone                        = true;
            }
            if ( smallerCorner > largerCorner )
            {
                this->latToInds                       = lat;
                latCornersDone                        = true;
            }
            if ( smallerCorner == largerCorner )
            {
                if      ( lat < this->latFromInds )   { this->latFromInds = lat; latCornersDone = true; }
                else if ( lat > this->latToInds   )   { this->latToInds   = lat; latCornersDone = true; }
            }
        }
        
        
        //============================================ Check if longitude boundaries need to be modified
        if ( ( this->lonFromInds <= this->lonToInds ) && !( ( lon >= this->lonFromInds ) && ( lon <= this->lonToInds ) ) )
        {
            //======================================== Lattitude is outside of group boundaries
            smallerCorner                             = angularDistanceWithBorders ( static_cast< proshade_signed > ( lon ), static_cast< proshade_signed > ( this->lonFromInds ) );
            largerCorner                              = angularDistanceWithBorders ( static_cast< proshade_signed > ( lon ), static_cast< proshade_signed > ( this->lonToInds   ) );
            
            if ( smallerCorner < largerCorner )
            {
                this->lonFromInds                     = lon;
                lonCornersDone                        = true;
            }
            if ( smallerCorner > largerCorner )
            {
                this->lonToInds                       = lon;
                lonCornersDone                        = true;
            }
            if ( smallerCorner == largerCorner )
            {
                if      ( lon < this->lonFromInds )   { this->lonFromInds = lon; lonCornersDone = true; }
                else if ( lon > this->lonToInds   )   { this->lonToInds   = lon; lonCornersDone = true; }
            }
        }
        
        if ( ( this->lonFromInds >  this->lonToInds ) && !( ( lon >= this->lonFromInds ) || ( lon <= this->lonToInds ) ) )
        {
            //======================================== Lattitude is outside of group boundaries
            smallerCorner                             = angularDistanceWithBorders ( static_cast< proshade_signed > ( lon ), static_cast< proshade_signed > ( this->lonFromInds ) );
            largerCorner                              = angularDistanceWithBorders ( static_cast< proshade_signed > ( lon ), static_cast< proshade_signed > ( this->lonToInds   ) );
            
            if ( smallerCorner < largerCorner )
            {
                this->lonFromInds                     = lon;
                lonCornersDone                        = true;
            }
            if ( smallerCorner > largerCorner )
            {
                this->lonToInds                       = lon;
                lonCornersDone                        = true;
            }
            if ( smallerCorner == largerCorner )
            {
                if      ( lon < this->lonFromInds )   { this->lonFromInds = lon; lonCornersDone = true; }
                else if ( lon > this->lonToInds   )   { this->lonToInds   = lon; lonCornersDone = true; }
            }
        }
        
        //============================================ Modify corner positions
        if ( latCornersDone )
        {
            this->latFrom                             = static_cast<proshade_double> ( this->latFromInds ) * this->latSampling;
            this->latTo                               = static_cast<proshade_double> ( this->latToInds   ) * this->latSampling;
        }
        
        if ( lonCornersDone )
        {
            this->lonFrom                             = static_cast<proshade_double> ( this->lonFromInds ) * this->lonSampling;
            this->lonTo                               = static_cast<proshade_double> ( this->lonToInds   ) * this->lonSampling;
        }
        
        //============================================ Compute corner vectors
        this->computeCornerPositions                  ( );
        hlpSS << "LAT " << this->latFromInds << " - " << this->latToInds << " and LON " << this->lonFromInds << " - " << this->lonToInds << " ( peak position LAT " << lat << " LON " << lon << " )";
        ProSHADE_internal_messages::printProgressMessage  ( verbose, 7, hlpSS.str(), messageShift );
        
        //============================================ If new sphere, add it to the list
        bool isSphereNew                              = true;
        for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( this->spherePositions.size() ); iter++ ) { if ( this->spherePositions.at(iter) == sphPos ) { isSphereNew = false; } }
        if ( isSphereNew ) { ProSHADE_internal_misc::addToUnsignVector ( &this->spherePositions, sphPos ); }
    }
    
    //================================================ Done
    return                                            ( peakAdded );
    
}

/*! \brief Accessor function for the private variable latFromInds.
 
    \param[out] latFromInds The lattitude index start for the group.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::getLatFromIndices ( void )
{
    //================================================ Done
    return                                            ( this->latFromInds );
    
}

/*! \brief Accessor function for the private variable latToInds.
 
    \param[out] latToInds The lattitude index end for the group.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::getLatToIndices ( void )
{
    //================================================ Done
    return                                            ( this->latToInds );
    
}

/*! \brief Accessor function for the private variable lonFromInds.
 
    \param[out] lonFromInds The longitude index start for the group.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::getLonFromIndices ( void )
{
    //================================================ Done
    return                                            ( this->lonFromInds );
    
}

/*! \brief Accessor function for the private variable lonToInds.
 
    \param[out] lonToInds The longitude index end for the group.
 */
proshade_double ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::getLonToIndices ( void )
{
    //================================================ Done
    return                                            ( this->lonToInds );
    
}

/*! \brief Accessor function for the private variable spherePositions.
 
    \param[out] spherePositions A vector of all angles (spheres) indices present in this group.
 */
std::vector<proshade_unsign> ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::getSpherePositions ( void )
{
    //================================================ Done
    return                                            ( this->spherePositions );
    
}

/*! \brief Function detecting cyclic point groups with a particular fold in a peak group.
 
    This function is a simplification of the findCyclicPointGroups function for the cases where the required fold is known. It simply assumes that all the
    supplied mapped spheres are to be used to find the fold, i.e. that fold-1 is equal to the length of sphereVals. With this assumption, the function can
    go directly for finding the peak index with highest peak height sum.
 
    At this point, this function can also optionally do bi-cubic interpolation around this index with highest peak sum to try to improve the symmetry
    axis by searching between the lattitude and longitude indices. Finally, this function will create the ProSHADE formatted array of symmetry group
    information and save it into the supplied vector, terminating thereafter.
 
    \warning This function  assumes that the supplied sphereVals argument contains only the spheres relating to the
    required fold and no other spheres - this assumption does not hold if the convertRotationFunction() function was called - consider yourself warned.
 
    \param[in] sphereVals A vector of spheres with mapped rotation function values.
    \param[in] detectedCs A vector of double pointers pointer to which any detected axis will be added in the ProSHADE format - [0] = fold, [1] = x-axis, [2] = y-axis, [3] = z-axis, [4] = angle, [5] = average peak height.
    \param[in] bicubicInterp Should the bicubic interpolation between the peak indices be done?
    \param[in] fold The fold for which we are searching for cyclic point groups.
    \param[in] verbose The verbosity of the run.
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
 */
void ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::findCyclicPointGroupsGivenFold ( std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*> sphereVals, std::vector < proshade_double* >* detectedCs, bool bicubicInterp, proshade_unsign fold, proshade_signed verbose, proshade_signed messageShift )
{
    //================================================ Check that this peak group has all the angles
    if ( ( fold - 1 ) != spherePositions.size() ) { return ; }
    
    //================================================ Initialise variables
    proshade_double bestPosVal, bestLatInd, bestLonInd;
    std::vector< proshade_unsign > spheresFormingFold;
    
    //================================================ Set all supplied spheres to be required to form the fold
    for ( proshade_unsign shIt = 0; shIt < static_cast<proshade_unsign> ( sphereVals.size() ); shIt++ )
    {
        ProSHADE_internal_misc::addToUnsignVector     ( &spheresFormingFold, shIt );
    }
    
    //================================================ Find the index with the highest peak height sum
    this->getBestIndexForFold                         ( &bestPosVal, &bestLatInd, &bestLonInd, &spheresFormingFold, sphereVals );
    
    //================================================ Optimise by bicubic interpolation if required
    if ( bicubicInterp )
    {
        ProSHADE_internal_maths::optimiseAxisBiCubicInterpolation ( &bestLatInd, &bestLonInd, &bestPosVal, &spheresFormingFold, &sphereVals );
    }
    
    //================================================ Create ProSHADE symmetry axis array and save it
    proshade_double* detectedSymmetry                 = new proshade_double[7];
    ProSHADE_internal_misc::checkMemoryAllocation ( detectedSymmetry, __FILE__, __LINE__, __func__ );
    
    detectedSymmetry[0]                               = static_cast<proshade_double> ( fold );
    detectedSymmetry[1]                               = 1.0 * std::sin ( bestLonInd * this->lonSampling ) * std::cos ( bestLatInd * this->latSampling );
    detectedSymmetry[2]                               = 1.0 * std::sin ( bestLonInd * this->lonSampling ) * std::sin ( bestLatInd * this->latSampling );
    detectedSymmetry[3]                               = 1.0 * std::cos ( bestLonInd * this->lonSampling );
    detectedSymmetry[4]                               = ( 2.0 * M_PI ) / detectedSymmetry[0];
    detectedSymmetry[5]                               = ( bestPosVal - 1.0 ) / ( detectedSymmetry[0] - 1 );
    detectedSymmetry[6]                               = -std::numeric_limits < proshade_double >::infinity();
    
    //================================================ Make sure max is positive
    const FloatingPoint< proshade_double > lhs1 ( std::max ( std::abs ( detectedSymmetry[1] ), std::max ( std::abs ( detectedSymmetry[2] ), std::abs ( detectedSymmetry[3] ) ) ) );
    const FloatingPoint< proshade_double > rhs1 ( std::abs ( detectedSymmetry[1] ) );
    const FloatingPoint< proshade_double > rhs2 ( std::abs ( detectedSymmetry[2] ) );
    const FloatingPoint< proshade_double > rhs3 ( std::abs ( detectedSymmetry[3] ) );
    if ( ( lhs1.AlmostEquals ( rhs1 ) && ( detectedSymmetry[1] < 0.0 ) ) ||
         ( lhs1.AlmostEquals ( rhs2 ) && ( detectedSymmetry[2] < 0.0 ) ) ||
         ( lhs1.AlmostEquals ( rhs3 ) && ( detectedSymmetry[3] < 0.0 ) ) )
    {
        detectedSymmetry[1]                          *= -1.0;
        detectedSymmetry[2]                          *= -1.0;
        detectedSymmetry[3]                          *= -1.0;
        detectedSymmetry[4]                          *= -1.0;
    }
    
    //================================================ Check for decimal underflows. They are not an issue as rounding to zero is fine, but having negative sign to zero will switch the signs and cause problems.
    const FloatingPoint< proshade_double > llhs1 ( detectedSymmetry[1] ), llhs2 ( detectedSymmetry[2] ), llhs3 ( detectedSymmetry[3] ), rrhs1 ( 0.0 );
    if ( llhs1.AlmostEquals ( rrhs1 ) ) { detectedSymmetry[1] = 0.0; }
    if ( llhs2.AlmostEquals ( rrhs1 ) ) { detectedSymmetry[2] = 0.0; }
    if ( llhs3.AlmostEquals ( rrhs1 ) ) { detectedSymmetry[3] = 0.0; }
    
    //================================================ Save detected point group
    ProSHADE_internal_misc::addToDblPtrVector         ( detectedCs, detectedSymmetry );
    
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "Detected group with fold " << detectedSymmetry[0] << " along axis " << detectedSymmetry[1] << " ; " << detectedSymmetry[2] << " ; " << detectedSymmetry[3] << " and with peak height " << detectedSymmetry[5];
    ProSHADE_internal_messages::printProgressMessage ( verbose, 4, hlpSS.str(), messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes all angles present in this peak group and finds the set of unique angle differeces.
 
    \param[in] angDiffs A pointer to a vector to which all angle differences will be saved into.
    \param[in] sphereVals A vector of spheres with mapped rotation function values.
 */
void ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::getAllAngleDifferences ( std::vector< proshade_double >* angDiffs, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*> sphereVals )
{
    //================================================ Initialise local variables
    std::vector< proshade_double > angs;
    
    //================================================ Find all present angles
    for ( proshade_unsign shPos = 0; shPos < static_cast<proshade_unsign> ( this->spherePositions.size() ); shPos++ )
    {
        ProSHADE_internal_misc::addToDoubleVector     ( &angs, sphereVals.at(this->spherePositions.at(shPos))->getRadius() );
    }
    
    //================================================ Find all angle differences
    for ( proshade_unsign ang1It = 0; ang1It < static_cast<proshade_unsign> ( angs.size() ); ang1It++ )
    {
        for ( proshade_unsign ang2It = 1; ang2It < static_cast<proshade_unsign> ( angs.size() ); ang2It++ )
        {
            //======================================== Use unique combinations only
            if ( ang1It >= ang2It ) { continue; }
            
            //======================================== Add angle difference rounded to 5 decimal places
            ProSHADE_internal_misc::addToDoubleVector ( angDiffs, std::floor ( std::abs ( angs.at(ang1It) - angs.at(ang2It) ) * 100000.0 ) / 100000.0 );
        }
    }
    
    //================================================ Sort and remove duplicates
    std::sort                                         ( (*angDiffs).begin(), (*angDiffs).end() );
  (*angDiffs).erase                                   ( std::unique ( (*angDiffs).begin(), (*angDiffs).end() ), (*angDiffs).end() );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function angle differences and creates a list of folds that may be present in the group.
 
    The function starts by taking each detected angle difference and checking how well it divides a circle (2 Pi). The remainder of
    this division is then check against a tolerance threshold, which takes into account how many peaks in the rotation function space
    is the distance off from the theoretical exact value.
 
    \param[in] angDiffs A pointer to a vector containing all the unique angle differences for this peak group.
    \param[in] foldsToTry A pointer to a vector to which the predicted fold to try to find are to be saved into.
 */
void ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::getAllPossibleFolds ( std::vector< proshade_double >* angDiffs, std::vector< proshade_unsign >* foldsToTry )
{
    //================================================ Initialise local variables
    proshade_double divRem, divBasis, symmErr, angTolerance, angToleranceNext;
    proshade_double peakErr                           = ( M_PI * 2.0 ) / ( static_cast<proshade_double> ( this->dimension ) );
    
    //================================================ For each angle difference in the group
    for ( proshade_unsign diffIt = 0; diffIt < static_cast<proshade_unsign> ( angDiffs->size() ); diffIt++ )
    {
        //============================================ Find the basis and remainder of the 2pi/dist equation
        divRem                                        = std::modf ( static_cast<proshade_double> ( ( 2.0 * M_PI ) / std::abs ( angDiffs->at(diffIt) ) ), &divBasis );
        
        //============================================ If the remainder would be smaller for larger basis, so change the basis
        if ( divRem > 0.5 )
        {
            divRem                                   -= 1.0;
            divBasis                                 += 1.0;
        }
        
        //============================================ Remove fold 1, that is not really what we are after here ...
        const FloatingPoint< proshade_double > lhs1 ( divBasis ), rhs1 ( 1.0 );
        if ( lhs1.AlmostEquals ( rhs1 ) ) { continue; }
        
        //============================================ Is there enough angles in the group for such a fold?
        if ( static_cast< proshade_double > ( this->spherePositions.size() ) < ( divBasis - 1.0 ) ) { continue; }
        
        //============================================ Determine errors on peaks and on folds
        symmErr                                       = divRem * ( ( 2.0 * M_PI ) / static_cast<proshade_double> ( divBasis ) );
        angTolerance                                  = std::abs( symmErr / peakErr );
        angToleranceNext                              = ( ( ( 2.0 * M_PI ) / static_cast<proshade_double> ( divBasis ) ) - ( ( 2.0 * M_PI ) / static_cast<proshade_double> ( divBasis + 1 ) ) ) / peakErr;

        //============================================ Is remainder small enough?
        if ( angTolerance < std::max ( 3.0, ( 0.1 / peakErr ) ) )
        {
            //======================================== Is the next symmetry close enough? If so, test previous and next folds as well.
            if ( angToleranceNext < std::max ( 1.5, ( 0.1 / peakErr ) ) )
            {
                //==================================== The next fold would pass as well. Use one previous and one following fold as well to cover for errors
                ProSHADE_internal_misc::addToUnsignVector ( foldsToTry, static_cast< proshade_unsign > ( divBasis - 1 ) );
                ProSHADE_internal_misc::addToUnsignVector ( foldsToTry, static_cast< proshade_unsign > ( divBasis + 1 ) );
            }
            
            //======================================== This fold seems reasonable, save it
            ProSHADE_internal_misc::addToUnsignVector ( foldsToTry, static_cast< proshade_unsign > ( divBasis ) );
        }
    }
    
    //================================================ Sort and remove duplicates
    std::sort                                         ( (*foldsToTry).begin(), (*foldsToTry).end(), std::greater <proshade_unsign>() );
    foldsToTry->erase                                 ( std::unique ( (*foldsToTry).begin(), (*foldsToTry).end() ), (*foldsToTry).end() );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function simply finds the indices of the spheres which form the requested form.
 
    \param[in] foldToTry The value of the fold for which all the required spheres are to be sought.
    \param[in] spheresFormingFold A pointer to vector to which the sphere indices of spheres forming this fold will be saved into.
    \param[in] sphereVals A vector of spheres with mapped rotation function values.
    \param[in] sphereAngleTolerance The tolerance for how different the sphere angle can be for the sphere to be still considered.
 */
void ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::getSpheresFormingFold ( proshade_unsign foldToTry, std::vector< proshade_unsign >* spheresFormingFold, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*> sphereVals, proshade_double sphereAngleTolerance )
{
    //================================================ Initialise local variables
    proshade_double soughtAngle, minSphereVal;
    proshade_unsign minSpherePos                      = 0;
    
    //================================================ Generate expected angles and check if close-by sphere exists
    for ( proshade_double fIt = 1.0; fIt < static_cast<proshade_double> ( foldToTry ); fIt += 1.0 )
    {
        //============================================ Set variables for the iteration
        minSphereVal                                  = 999.9;
        soughtAngle                                   = fIt * ( 2.0 * M_PI / static_cast<proshade_double> ( foldToTry ) );
        
        //============================================ Find the closest sphere passing conditions
        for ( proshade_unsign angsIt = 0; angsIt < static_cast<proshade_unsign> ( this->spherePositions.size() ); angsIt++ )
        {
            if ( std::abs ( sphereVals.at(this->spherePositions.at(angsIt))->getRepresentedAngle() - soughtAngle ) < sphereAngleTolerance )
            {
                if ( minSphereVal > 1.0 - std::cos ( std::abs ( sphereVals.at(this->spherePositions.at(angsIt))->getRepresentedAngle() - soughtAngle ) ) )
                {
                    minSphereVal                      = 1.0 - std::cos ( std::abs ( sphereVals.at(this->spherePositions.at(angsIt))->getRepresentedAngle() - soughtAngle ) );
                    minSpherePos                      = angsIt;
                }
            }
        }
       
        //============================================ If no passing sphere, test next fold
        const FloatingPoint< proshade_double > lhs1 ( minSphereVal ), rhs1 ( 999.9 );
        if ( lhs1.AlmostEquals ( rhs1 ) ) { break; }
        
        //============================================ Save best position
        ProSHADE_internal_misc::addToUnsignVector     ( spheresFormingFold, this->spherePositions.at(minSpherePos) );
    }
    
    //================================================ Done
    return ;
    
}


/*! \brief Function which does simple search through all peak groups indices and saves the index with the highest peak height sum over all spheres.
 
    \param[in] bestPosVal Pointer to double where the highest sum of heights will be stored.
    \param[in] bestLatInd Pointer to double where the highest values lattitude index will be held.
    \param[in] bestLonInd Pointer to double where the highest values longitude index will be held.
    \param[in] spheresFormingFold A vector pointer to a vector containing the indices of the spheres forming this fold.
    \param[in] sphereVals A vector of spheres with mapped rotation function values.
 */
void ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup::getBestIndexForFold ( proshade_double* bestPosVal, proshade_double* bestLatInd, proshade_double* bestLonInd, std::vector< proshade_unsign >* spheresFormingFold, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*> sphereVals )
{
    //================================================ Initialise variables
    proshade_double curPosVal;
   *bestPosVal                                        = -1.0;
    if ( this->latFromInds > this->latToInds ) { this->latToInds += static_cast< proshade_double > ( this->dimension ); }
    if ( this->lonFromInds > this->lonToInds ) { this->lonToInds += static_cast< proshade_double > ( this->dimension ); }
    
    //================================================ Compute the best average peak height axis for peak indices
    for ( proshade_unsign latIt = static_cast< proshade_unsign > ( this->latFromInds ); latIt <= static_cast< proshade_unsign > ( this->latToInds ); latIt++ )
    {
        //============================================ Deal with boundaries
        if ( latIt >= this->dimension ) { latIt -= this->dimension; this->latToInds -= static_cast< proshade_double > ( this->dimension ); }
        
        for ( proshade_unsign lonIt = static_cast< proshade_unsign > ( this->lonFromInds ); lonIt <= static_cast< proshade_unsign > ( this->lonToInds ); lonIt++ )
        {
            //======================================== Deal with boundaries
            if ( lonIt >= this->dimension ) { lonIt -= this->dimension; this->lonToInds -= static_cast< proshade_double > ( this->dimension ); }
            
            //======================================== Initialise variables
            curPosVal                                 = 1.0;
            
            //======================================== Find this indices value
            for ( proshade_unsign sphIt = 0; sphIt < static_cast<proshade_unsign> ( spheresFormingFold->size() ); sphIt++ )
            {
                curPosVal                            += sphereVals.at(spheresFormingFold->at(sphIt))->getSphereLatLonPosition ( latIt, lonIt );
            }
            
            //======================================== If best, save it
            if ( curPosVal > *bestPosVal )
            {
               *bestPosVal                            = curPosVal;
               *bestLatInd                            = static_cast< proshade_double > ( latIt );
               *bestLonInd                            = static_cast< proshade_double > ( lonIt );
            }
        }
    }
    
    //================================================ Done
    return ;
    
}
