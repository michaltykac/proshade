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
    \version   0.7.4.4
    \date      OCT 2020
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
    this->sphereWidth                                 = ( spherePos->at(0) + spherePos->at(1) ) / 2.0;
    this->sphereRadius                                = spherePos->at(shOrder);
    
    //================================================ Determine shell ranges in angstroms
    proshade_double maxDist                           = 0.0;
    if ( shOrder == static_cast<proshade_unsign> ( spherePos->size() - 1 ) ) { maxDist = static_cast<proshade_double> ( spherePos->at(spherePos->size()-1) + ( spherePos->at(1) - spherePos->at(0) ) ); }
    else { maxDist = static_cast<proshade_double> ( ( spherePos->at(shOrder) + spherePos->at(shOrder+1) ) / 2.0 ); }
    
    //================================================ Set the max range
    this->maxSphereRange                              = static_cast<proshade_single> ( 2.0 * maxDist );
    
    //================================================ Set map sampling rates
    this->xDimSampling                                = xSize / static_cast<proshade_single> (xDimMax);
    this->yDimSampling                                = ySize / static_cast<proshade_single> (yDimMax);
    this->zDimSampling                                = zSize / static_cast<proshade_single> (zDimMax);
    
    //================================================ Get maximum circumference
    proshade_unsign maxCircumference                  = this->getMaxCircumference ( xDimMax, yDimMax, zDimMax, this->maxSphereRange, xSize, ySize, zSize );
    
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
    this->mappedDataRot                               = NULL;
    
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
    
    if ( this->mappedDataRot != NULL )
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
    \param[in] xSize The size of x axis of whole map in angstroms.
    \param[in] ySize The size of y axis of whole map in angstroms.
    \param[in] zSize The size of z axis of whole map in angstroms.
 */
proshade_unsign ProSHADE_internal_spheres::ProSHADE_sphere::getMaxCircumference ( proshade_unsign xDimMax, proshade_unsign yDimMax, proshade_unsign zDimMax, proshade_double maxRange, proshade_single xSize, proshade_single ySize, proshade_single zSize )
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
    if ( xToInd > static_cast<proshade_signed> ( xDimMax ) ) { xToInd = xDimMax; }
    if ( yToInd > static_cast<proshade_signed> ( yDimMax ) ) { yToInd = yDimMax; }
    if ( zToInd > static_cast<proshade_signed> ( zDimMax ) ) { zToInd = zDimMax; }
    
    //================================================ Get dim sizes
    proshade_unsign xDimSZ                            = xToInd - xFromInd;
    proshade_unsign yDimSZ                            = yToInd - yFromInd;
    proshade_unsign zDimSZ                            = zToInd - zFromInd;
    
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
            xRelative                                 = ( x - ( ( xBottom - static_cast<proshade_signed> ( ( ( xDimMax ) / 2 ) ) ) * this->xDimSampling ) ) / this->xDimSampling;
            this->interpolateAlongFirst               ( c000, c001, c010, c011, c100, c101, c110, c111, &c00, &c01, &c10, &c11, xRelative );
            
            //======================================== Interpolate along Y axis
            yRelative                                 = ( y - ( ( yBottom - static_cast<proshade_signed> ( ( ( yDimMax ) / 2 ) ) ) * this->yDimSampling ) ) / this->yDimSampling;
            this->interpolateAlongSecond              ( c00, c01, c10, c11, &c0, &c1, yRelative );

            //======================================== Save the resulting value
            zRelative                                 = ( z - ( ( zBottom - static_cast<proshade_signed> ( ( ( zDimMax ) / 2 ) ) ) * this->zDimSampling ) ) / this->zDimSampling;
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
    proshade_signed posIter                           = zPos + zDimMax * ( yPos  + yDimMax * xPos );
    
    //================================================ Check of out-of-bounds shells
    if ( ( xPos < 0 ) || ( xPos >= static_cast<proshade_signed> ( xDimMax ) ) ) { return ( false ); }
    if ( ( yPos < 0 ) || ( yPos >= static_cast<proshade_signed> ( yDimMax ) ) ) { return ( false ); }
    if ( ( zPos < 0 ) || ( zPos >= static_cast<proshade_signed> ( zDimMax ) ) ) { return ( false ); }
    
    //================================================ Get the interpolation values
    interpVec->at(0)                                  = static_cast<proshade_double> ( xPos * this->xDimSampling );
    interpVec->at(1)                                  = static_cast<proshade_double> ( yPos * this->yDimSampling );
    interpVec->at(2)                                  = static_cast<proshade_double> ( zPos * this->zDimSampling );
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
   *x                                                 = this->sphereRadius * cos( ( lonCO->at(thetaIt) + lonCO->at(thetaIt+1) ) / 2.0 ) *
                                                                             cos( ( latCO->at(phiIt)   + latCO->at(phiIt+1)   ) / 2.0 );
   *y                                                 = this->sphereRadius * sin( ( lonCO->at(thetaIt) + lonCO->at(thetaIt+1) ) / 2.0 ) *
                                                                             cos( ( latCO->at(phiIt)   + latCO->at(phiIt+1)   ) / 2.0 );
   *z                                                 = this->sphereRadius * sin( ( latCO->at(phiIt)   + latCO->at(phiIt+1)   ) / 2.0 );
    
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
   *xBottom                                           = std::floor ( (x / this->xDimSampling) ) + (xDimMax/2);
   *yBottom                                           = std::floor ( (y / this->yDimSampling) ) + (yDimMax/2);
   *zBottom                                           = std::floor ( (z / this->zDimSampling) ) + (zDimMax/2);
            
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
    c00->at(0)                                        = ( this->xDimSampling * xd ) + c000.at(0);
    c00->at(1)                                        = c000.at(1);
    c00->at(2)                                        = c000.at(2);
    c00->at(3)                                        = ( c000.at(3) * ( 1.0 - xd ) ) + ( c100.at(3) * xd );
    
    //================================================ Interpolate for the less more point
    c01->at(0)                                        = ( this->xDimSampling * xd ) + c001.at(0);
    c01->at(1)                                        = c001.at(1);
    c01->at(2)                                        = c001.at(2);
    c01->at(3)                                        = ( c001.at(3) * ( 1.0 - xd ) ) + ( c101.at(3) * xd );
    
    //================================================ Interpolate for the more less point
    c10->at(0)                                        = ( this->xDimSampling * xd ) + c010.at(0);
    c10->at(1)                                        = c010.at(1);
    c10->at(2)                                        = c010.at(2);
    c10->at(3)                                        = ( c010.at(3) * ( 1.0 - xd ) ) + ( c110.at(3) * xd );
    
    //================================================ Interpolate for the more more point
    c11->at(0)                                        = ( this->xDimSampling * xd ) + c011.at(0);
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
    c0->at(1)                                         = ( this->yDimSampling * yd ) + c00.at(1);
    c0->at(2)                                         = c00.at(2);
    c0->at(3)                                         = ( c00.at(3) * ( 1.0 - yd ) ) + ( c10.at(3) * yd );
    
    //================================================ Interpolate for the more point
    c1->at(0)                                         = c01.at(0);
    c1->at(1)                                         = ( this->yDimSampling * yd ) + c01.at(1);
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
    return                                            ( std::ceil ( circumference / 2 ) );
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
    proshade_single ret                               = static_cast<proshade_single> ( resolution / 2.0 );
    
    //================================================ Make sure at least 10 shells will exist
    while ( std::floor ( maxMapRange / ret ) < 10 )
    {
        ret                                          /= 2.0;
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
    proshade_double sphereDistanceAsFractionOfTotal = static_cast<proshade_double> ( sphereDist ) / ( maxMapRange / 2.0 );
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
    this->mappedDataRot                               = NULL;
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
 
    ...
 
    \param[in] 
    \param[out] X Data object with all values set and ready to accept the rotation function values.
 */
ProSHADE_internal_spheres::ProSHADE_rotFun_sphere::ProSHADE_rotFun_sphere ( proshade_double rad, proshade_double radRange, proshade_unsign dim, proshade_double repAng )
{
    //================================================ Set internal values
    this->radius                                      = rad;
    this->angularDim                                  = dim;
    this->radiusMin                                   = this->radius - ( radRange / 2.0 );
    this->radiusMax                                   = this->radius + ( radRange / 2.0 );
    this->representedAngle                            = repAng;
    
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
    proshade_double lonSampling                       = ( M_PI       ) / static_cast<proshade_double> ( this->angularDim );
    proshade_double latSampling                       = ( M_PI * 2.0 ) / static_cast<proshade_double> ( this->angularDim );
    
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
            ProSHADE_internal_maths::getEulerZXZFromAngleAxis ( cX, cY, cZ, this->representedAngle, &eulerAlpha, &eulerBeta, &eulerGamma, this->angularDim );
            
            //======================================== Convert to SOFT map position (decimal, not indices)
            ProSHADE_internal_maths::getSOFTPositionFromEulerZXZ ( this->angularDim / 2, eulerAlpha, eulerBeta, eulerGamma, &mapX, &mapY, &mapZ );
            
            //======================================== Find lower and higher points and deal with boundaries
            xBottom = std::floor ( mapX ); if ( xBottom < 0.0 ) { xBottom += this->angularDim; } if ( xBottom >= this->angularDim ) { xBottom -= this->angularDim; }
            yBottom = std::floor ( mapY ); if ( yBottom < 0.0 ) { yBottom += this->angularDim; } if ( yBottom >= this->angularDim ) { yBottom -= this->angularDim; }
            zBottom = std::floor ( mapZ ); if ( zBottom < 0.0 ) { zBottom += this->angularDim; } if ( zBottom >= this->angularDim ) { zBottom -= this->angularDim; }
            xTop = std::ceil ( mapX ); if ( xTop < 0.0 ) { xTop += this->angularDim; } if ( xTop >= this->angularDim ) { xTop -= this->angularDim; }
            yTop = std::ceil ( mapY ); if ( yTop < 0.0 ) { yTop += this->angularDim; } if ( yTop >= this->angularDim ) { yTop -= this->angularDim; }
            zTop = std::ceil ( mapZ ); if ( zTop < 0.0 ) { zTop += this->angularDim; } if ( zTop >= this->angularDim ) { zTop -= this->angularDim; }
            
            //======================================== Start X interpolation - bottom, bottom, bottom
            mapIndex                                  = zBottom + this->angularDim * ( yBottom + this->angularDim * xBottom );
            c000                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - bottom, bottom, top
            mapIndex                                  = zTop    + this->angularDim * ( yBottom + this->angularDim * xBottom );
            c001                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - bottom, top, bottom
            mapIndex                                  = zBottom + this->angularDim * ( yTop    + this->angularDim * xBottom );
            c010                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - bottom, top, top
            mapIndex                                  = zTop    + this->angularDim * ( yTop    + this->angularDim * xBottom );
            c011                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - top, bottom, bottom
            mapIndex                                  = zBottom + this->angularDim * ( yBottom + this->angularDim * xTop    );
            c100                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - top, bottom, top
            mapIndex                                  = zTop    + this->angularDim * ( yBottom + this->angularDim * xTop    );
            c101                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - top, top, bottom
            mapIndex                                  = zBottom + this->angularDim * ( yTop    + this->angularDim * xTop    );
            c110                                      = pow( rotFun[mapIndex][0], 2.0 ) + pow( rotFun[mapIndex][1], 2.0 );
            
            //======================================== X interpolation - top, top, top
            mapIndex                                  = zTop    + this->angularDim * ( yTop    + this->angularDim * xTop    );
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
            mapIndex                                  = lonIt + ( latIt * static_cast<proshade_unsign> ( this->angularDim ) );
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
 
    ...
 
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
    latBottom = std::floor ( lattitude ); if ( latBottom < 0.0 ) { latBottom += this->angularDim; } if ( latBottom >= this->angularDim ) { latBottom -= this->angularDim; }
    lonBottom = std::floor ( longitude ); if ( lonBottom < 0.0 ) { lonBottom += this->angularDim; } if ( lonBottom >= this->angularDim ) { lonBottom -= this->angularDim; }
    latTop    = std::ceil  ( lattitude ); if ( latTop    < 0.0 ) { latTop    += this->angularDim; } if ( latTop    >= this->angularDim ) { latTop    -= this->angularDim; }
    lonTop    = std::ceil  ( longitude ); if ( lonTop    < 0.0 ) { lonTop    += this->angularDim; } if ( lonTop    >= this->angularDim ) { lonTop    -= this->angularDim; }
    
    //================================================ Interpolate
    gridIndex                                         = lonBottom + ( latBottom * static_cast<proshade_unsign> ( this->angularDim ) );
    c00                                               = this->axesValues[gridIndex];
    
    gridIndex                                         = lonBottom + ( latTop    * static_cast<proshade_unsign> ( this->angularDim ) );
    c01                                               = this->axesValues[gridIndex];
    
    gridIndex                                         = lonTop    + ( latBottom * static_cast<proshade_unsign> ( this->angularDim ) );
    c10                                               = this->axesValues[gridIndex];
    
    gridIndex                                         = lonTop    + ( latTop    * static_cast<proshade_unsign> ( this->angularDim ) );
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
