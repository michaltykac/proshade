/*! \file ProSHADE_io.cpp
    \brief This source file contains the functions required for specifc data format manipulations.
 
    The functions in this source file are required for detection of supported input and output formats as well as reading and writing into the
    supported formats. They are the lowest level functions ProSHADE has for disk access.
 
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
#include "ProSHADE_io.hpp"

/*! \brief Function determining if the input data type is PDB.
 
    This function checks if the input file is a PDB file and can be read by the gemmi library.
 
    \param[in] fName The file name of the file for which the type should be determined.
    \param[out] X Bool value true if the file is a PDB file readable by gemmi and false otherwise.
 */
bool ProSHADE_internal_io::isFilePDB ( std::string fName )
{
    //================================================ Try reading the file using Gemmi
    try
    {
        gemmi::Structure structure                    = gemmi::read_structure ( gemmi::MaybeGzipped ( fName ) );
    }
    catch ( std::runtime_error& e )
    {
        //============================================ Supress MSVC C4101 Unferenced variable warning / clang and gcc -Wunused-exception-parameter
        (void)e;
        
        //============================================ Read failed. Done
        return                                        ( false );
    }
    
    //================================================ Read successfull. Done
    return                                            ( true );
    
}

/*! \brief Function determining if the input data type is MAP.
 
    This function checks if the input file is a MAP file and can be read by the CMAP library.
 
    \param[in] fName The file name of the file for which the type should be determined.
    \param[out] X Bool value true if the file is a MAP file readable by CMAP and false otherwise.
 */
bool ProSHADE_internal_io::isFileMAP ( std::string fName )
{
    gemmi::Ccp4<float> map;
    try
    {
        map.read_ccp4                                 ( gemmi::MaybeGzipped (fName.c_str() ) );
    }
    catch ( std::runtime_error& e )
    {
        //============================================ Supress MSVC C4101 Unferenced variable warning / clang and gcc -Wunused-exception-parameter
        (void)e;
        
        //============================================ Failed to read the map
        return                                        ( false );
    }

    //================================================ Done
    return                                            ( true );
    
}

/*! \brief This function parses the CCP4 MAP file header as read in by gemmi.
 
    This function uses the gemmi Ccp4 object, which contains all the information read in from a MAP file (including the header), to parse out the ProSHADE required
    information from the header and saving it to the supplied variables.
 
    \param[in] map A gemmi Ccp4 objecct containing all the data read in from a MAP file.
    \param[in] xDimInds Address to a variable to save the x-axis size in indices.
    \param[in] yDimInds Address to a variable to save the y-axis size in indices.
    \param[in] zDimInds Address to a variable to save the z-axis size in indices.
    \param[in] xDim Address to a variable to save the x dimension size in angstroms.
    \param[in] yDim Address to a variable to save the y dimension size in angstroms.
    \param[in] zDim Address to a variable to save the z dimension size in angstroms.
    \param[in] aAng Address to a variable to save the a angle in degrees.
    \param[in] bAng Address to a variable to save the b angle in degrees.
    \param[in] cAng Address to a variable to save the c angle in degrees.
    \param[in] xFrom Address to a variable to save the starting index along the x-axis.
    \param[in] yFrom Address to a variable to save the starting index along the y-axis.
    \param[in] zFrom Address to a variable to save the starting index along the z-axis.
    \param[in] xAxOrigin Address to a variable to save the map origin positon along the x-axis.
    \param[in] yAxOrigin Address to a variable to save the map origin positon along the y-axis.
    \param[in] zAxOrigin Address to a variable to save the map origin positon along the z-axis.
    \param[in] xAxOrder Address to a variable to save the order of x axis.
    \param[in] yAxOrder Address to a variable to save the order of y axis.
    \param[in] zAxOrder Address to a variable to save the order of z axis.
    \param[in] xGridInds Address to a variable to save the grid indices count along the x-axis.
    \param[in] yGridInds Address to a variable to save the grid indices count along the y-axis.
    \param[in] zGridInds Address to a variable to save the grid indices count along the z-axis.
 */
void ProSHADE_internal_io::readInMapHeader ( gemmi::Ccp4<float> *map, proshade_unsign *xDimInds, proshade_unsign *yDimInds, proshade_unsign *zDimInds, proshade_single *xDim, proshade_single *yDim, proshade_single *zDim, proshade_single *aAng, proshade_single *bAng, proshade_single *cAng, proshade_signed *xFrom, proshade_signed *yFrom, proshade_signed *zFrom, proshade_signed *xAxOrigin, proshade_signed *yAxOrigin, proshade_signed *zAxOrigin, proshade_unsign *xAxOrder, proshade_unsign *yAxOrder, proshade_unsign *zAxOrder, proshade_unsign *xGridInds, proshade_unsign *yGridInds, proshade_unsign *zGridInds )
{
    //================================================ Read in the map file header
   *xDimInds                                          = static_cast<proshade_unsign> ( map->header_i32   ( 1  ) );
   *yDimInds                                          = static_cast<proshade_unsign> ( map->header_i32   ( 2  ) );
   *zDimInds                                          = static_cast<proshade_unsign> ( map->header_i32   ( 3  ) );
    
   *xFrom                                             = static_cast<proshade_signed> ( map->header_i32   ( 5  ) );
   *yFrom                                             = static_cast<proshade_signed> ( map->header_i32   ( 6  ) );
   *zFrom                                             = static_cast<proshade_signed> ( map->header_i32   ( 7  ) );
   
   *xDim                                              = static_cast<proshade_single> ( map->header_float ( 11 ) );
   *yDim                                              = static_cast<proshade_single> ( map->header_float ( 12 ) );
   *zDim                                              = static_cast<proshade_single> ( map->header_float ( 13 ) );
   
   *aAng                                              = static_cast<proshade_single> ( map->header_float ( 14 ) );
   *bAng                                              = static_cast<proshade_single> ( map->header_float ( 15 ) );
   *cAng                                              = static_cast<proshade_single> ( map->header_float ( 16 ) );
   
   *xAxOrigin                                         = static_cast<proshade_signed> ( map->header_i32   ( 50 ) ) + (*xFrom);
   *yAxOrigin                                         = static_cast<proshade_signed> ( map->header_i32   ( 51 ) ) + (*yFrom);
   *zAxOrigin                                         = static_cast<proshade_signed> ( map->header_i32   ( 52 ) ) + (*zFrom);
   
   *xAxOrder                                          = static_cast<proshade_unsign> ( map->header_i32   ( 17 ) );
   *yAxOrder                                          = static_cast<proshade_unsign> ( map->header_i32   ( 18 ) );
   *zAxOrder                                          = static_cast<proshade_unsign> ( map->header_i32   ( 19 ) );
   
   *xGridInds                                         = static_cast<proshade_unsign> ( map->header_i32   ( 8  ) );
   *yGridInds                                         = static_cast<proshade_unsign> ( map->header_i32   ( 9  ) );
   *zGridInds                                         = static_cast<proshade_unsign> ( map->header_i32   ( 10 ) );
    
    //================================================ Deal with sampling being different from cell size
    if ( *xGridInds != *xDimInds )
    {
        *xDim                                         = *xDim * ( static_cast<proshade_single> ( *xDimInds ) / static_cast<proshade_single> ( *xGridInds ) );
        *xGridInds                                    = *xDimInds;
    }
    
    if ( *yGridInds != *yDimInds )
    {
        *yDim                                         = *yDim * ( static_cast<proshade_single> ( *yDimInds ) / static_cast<proshade_single> ( *yGridInds ) );
        *yGridInds                                    = *yDimInds;
    }
    
    if ( *zGridInds != *zDimInds )
    {
        *zDim                                         = *zDim * ( static_cast<proshade_single> ( *zDimInds ) / static_cast<proshade_single> ( *zGridInds ) );
        *zGridInds                                    = *zDimInds;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function converts the gemmi Ccp4 object data to ProSHADE internal map representation.
 
    This function firstly allocates the required memory for the ProSHADE internal map representation variable according to the grid size. Then, it iterates over the axes in such a way, so that the resulting ProSHADE
    variable would have XYZ axis order independently on the axis order of the Ccp4 gemmi object. This should not be necessary as the gemmi setup function should have been called by now, but one never knows.
 
    \param[in] gemmiMap Pointer to a gemmi Ccp4 object containing the read in MAP file information.
    \param[in] map Pointer reference to a variable to save the map data.
    \param[in] xDimInds The size of x dimension in indices.
    \param[in] yDimInds The size of y dimension in indices.
    \param[in] zDimInds The size of z dimension in indices.
    \param[in] xAxOrder The order of the x-axis.
    \param[in] yAxOrder The order of the y-axis.
    \param[in] zAxOrder The order of the z-axis.
 */
void ProSHADE_internal_io::readInMapData ( gemmi::Ccp4<float> *gemmiMap, proshade_double*& map, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds, proshade_unsign xAxOrder, proshade_unsign yAxOrder, proshade_unsign zAxOrder )
{
    //================================================ Allocate internal variables
    proshade_unsign *axOrdArr                         = new proshade_unsign[3];
    proshade_unsign *axDimArr                         = new proshade_unsign[3];
    proshade_unsign arrPos                            = 0;
    
    //================================================ Check memory allocation and fill in values
    ProSHADE_internal_misc::checkMemoryAllocation     ( axOrdArr, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( axDimArr, __FILE__, __LINE__, __func__ );
    axDimArr[0]                                       = xDimInds;
    axDimArr[1]                                       = yDimInds;
    axDimArr[2]                                       = zDimInds;
    
    //================================================ Allocate the ProSHADE internal map variable memory
    map                                               = new proshade_double [xDimInds * yDimInds * zDimInds];
    ProSHADE_internal_misc::checkMemoryAllocation     ( map, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy read in data to internal map variable
    for ( axOrdArr[0] = 0; axOrdArr[0] < axDimArr[xAxOrder-1]; axOrdArr[0]++ )
    {
        for ( axOrdArr[1] = 0; axOrdArr[1] < axDimArr[yAxOrder-1]; axOrdArr[1]++ )
        {
            for ( axOrdArr[2] = 0; axOrdArr[2] < axDimArr[zAxOrder-1]; axOrdArr[2]++ )
            {
                arrPos                                = axOrdArr[2]  + axDimArr[zAxOrder-1] * ( axOrdArr[1]  + axDimArr[yAxOrder-1] * axOrdArr[0] );
                map[arrPos]                           = static_cast< proshade_double > ( gemmiMap->grid.get_value_q( static_cast< int > ( axOrdArr[xAxOrder-1] ),
                                                                                                                     static_cast< int > ( axOrdArr[yAxOrder-1] ),
                                                                                                                     static_cast< int > ( axOrdArr[zAxOrder-1] ) ) );
            }
        }
    }
    
    //================================================ Release internal variables memory
    delete[] axDimArr;
    delete[] axOrdArr;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function converts the gemmi Ccp4 object data to ProSHADE mask representation.
 
    This function firstly allocates the required memory for the ProSHADE mask representation variable according to the grid size. Then, it iterates over the axes in such a way, so that the resulting ProSHADE
    variable would have XYZ axis order independently on the axis order of the Ccp4 gemmi object. This should not be necessary as the gemmi setup function should have been called by now, but one never knows.
 
    \param[in] gemmiMap Pointer to a gemmi Ccp4 object containing the read in mask file information.
    \param[in] map Pointer reference to a variable to save the map data.
    \param[in] xDimInds The size of x dimension in indices.
    \param[in] yDimInds The size of y dimension in indices.
    \param[in] zDimInds The size of z dimension in indices.
    \param[in] xAxOrder The order of the x-axis.
    \param[in] yAxOrder The order of the y-axis.
    \param[in] zAxOrder The order of the z-axis.
 */
void ProSHADE_internal_io::readInMapData ( gemmi::Ccp4<int8_t> *gemmiMap, proshade_double*& map, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds, proshade_unsign xAxOrder, proshade_unsign yAxOrder, proshade_unsign zAxOrder )
{
    //================================================ Allocate internal variables
    proshade_unsign *axOrdArr                         = new proshade_unsign[3];
    proshade_unsign *axDimArr                         = new proshade_unsign[3];
    proshade_unsign arrPos                            = 0;
    
    //================================================ Check memory allocation and fill in values
    ProSHADE_internal_misc::checkMemoryAllocation     ( axOrdArr, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( axDimArr, __FILE__, __LINE__, __func__ );
    axDimArr[0]                                       = xDimInds;
    axDimArr[1]                                       = yDimInds;
    axDimArr[2]                                       = zDimInds;
    
    //================================================ Allocate the ProSHADE internal map variable memory
    map                                               = new proshade_double [xDimInds * yDimInds * zDimInds];
    ProSHADE_internal_misc::checkMemoryAllocation     ( map, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy read in data to internal map variable
    for ( axOrdArr[0] = 0; axOrdArr[0] < axDimArr[xAxOrder-1]; axOrdArr[0]++ )
    {
        for ( axOrdArr[1] = 0; axOrdArr[1] < axDimArr[yAxOrder-1]; axOrdArr[1]++ )
        {
            for ( axOrdArr[2] = 0; axOrdArr[2] < axDimArr[zAxOrder-1]; axOrdArr[2]++ )
            {
                arrPos                                = axOrdArr[2]  + axDimArr[zAxOrder-1] * ( axOrdArr[1]  + axDimArr[yAxOrder-1] * axOrdArr[0] );
                map[arrPos]                           = static_cast< proshade_double > ( gemmiMap->grid.get_value_q( static_cast< int > ( axOrdArr[xAxOrder-1] ),
                                                                                                                     static_cast< int > ( axOrdArr[yAxOrder-1] ),
                                                                                                                     static_cast< int > ( axOrdArr[zAxOrder-1] ) ) );
            }
        }
    }
    
    //================================================ Release internal variables memory
    delete[] axDimArr;
    delete[] axOrdArr;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function reads and applies the mask to the map.
 
    This function reads in the mask file and determines, if its dimensions are the same as the input map or not. If not,
    then the mask is re-sampled is such a way as to have the same number of indices as the input map. Once the map
    and mask dimensions are the same, the mask is applied by multiplying all map values by corresponding mask value.
    The map change is done in-place!
 
    \param[in] map Pointer to a the internal variable containing the map (it will be changed in place).
    \param[in] maskFile Filename of where the mask is to be read from.
    \param[in] xDimInds The size of x dimension in indices.
    \param[in] yDimInds The size of y dimension in indices.
    \param[in] zDimInds The size of z dimension in indices.
    \param[in] verbose How much std::out output would you like?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
    \param[in] maskArray An array of mask values (default nullptr) to be used instead of an input file.
    \param[in] maXInds The size of maskArray x dimension in indices (defaults to 0).
    \param[in] maYInds The size of maskArray y dimension in indices (defaults to 0).
    \param[in] maZInds The size of maskArray z dimension in indices (defaults to 0).
 */
void ProSHADE_internal_io::applyMask ( proshade_double*& map, std::string maskFile, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds, proshade_signed verbose, proshade_signed messageShift, proshade_double* maskArray, proshade_unsign maXInds, proshade_unsign maYInds, proshade_unsign maZInds )
{
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "Reading mask " << maskFile;
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, hlpSS.str() );
    
    //================================================ Are we reading from array or from file?
    if ( ( maskArray != nullptr ) && ( maXInds != 0 ) && ( maYInds != 0 ) && ( maZInds != 0 ) )
    {
        //============================================ Array it is!
        ProSHADE_internal_io::applyMaskFromArray      ( map, xDimInds, yDimInds, zDimInds, maskArray, maXInds, maYInds, maZInds, verbose, messageShift );
    }
    else
    {
        //============================================ Check if filename was given
        if ( maskFile == "" )                         { ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "No mask supplied." ); return; }
        
        //============================================ File it is! Open the mask
        gemmi::Ccp4<float> mask;
        mask.read_ccp4                                ( gemmi::MaybeGzipped ( maskFile.c_str() ) );
        
        //============================================ Convert to XYZ and create complete mask, if need be
        mask.setup                                    ( gemmi::GridSetup::ReorderOnly, 0 );
        
        //============================================ Read in the rest of the mask file header
        proshade_unsign xDI, yDI, zDI, xAOR, yAOR, zAOR, xGI, yGI, zGI;
        proshade_single xDS, yDS, zDS, aA, bA, cA;
        proshade_signed xF, yF, zF, xAO, yAO, zAO;
        ProSHADE_internal_io::readInMapHeader         ( &mask,
                                                        &xDI,  &yDI,  &zDI,
                                                        &xDS,  &yDS,  &zDS,
                                                        &aA,   &bA,   &cA,
                                                        &xF,   &yF,   &zF,
                                                        &xAO,  &yAO,  &zAO,
                                                        &xAOR, &yAOR, &zAOR,
                                                        &xGI,  &yGI,  &zGI );

        //============================================ Save the mask values to ProSHADE variable
        proshade_double* internalMask                 = nullptr;
        ProSHADE_internal_io::readInMapData           ( &mask, internalMask, xDI, yDI, zDI, xAOR, yAOR, zAOR );
        
        //============================================ Apply mask from array
        ProSHADE_internal_io::applyMaskFromArray      ( map, xDimInds, yDimInds, zDimInds, internalMask, xDI, yDI, zDI, verbose, messageShift );

        //============================================ Release the memory
        delete[] internalMask;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function applies the mask to the map.
 
    This function firstlydetermines if its dimensions of mask are the same as the input map or not. If not,
    then the mask is re-sampled is such a way as to have the same number of indices as the input map. Once the map
    and mask dimensions are the same, the mask is applied by multiplying all map values by corresponding mask value.
    The map change is done in-place!
 
    \param[in] map Pointer to a the internal variable containing the map (it will be changed in place).
    \param[in] xDimInds The size of the map x dimension in indices.
    \param[in] yDimInds The size of the map y dimension in indices.
    \param[in] zDimInds The size of the map z dimension in indices.
    \param[in] mask A pointer to 1D array of the mask values.
    \param[in] xDimIndsMsk The size of the mask x dimension in indices.
    \param[in] yDimIndsMsk The size of the mask y dimension in indices.
    \param[in] zDimIndsMsk The size of the mask z dimension in indices.
    \param[in] verbose How much std::out output would you like?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
 */
void ProSHADE_internal_io::applyMaskFromArray ( proshade_double*& map, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds, proshade_double*& mask, proshade_unsign xDimIndsMsk, proshade_unsign yDimIndsMsk, proshade_unsign zDimIndsMsk, proshade_signed verbose, proshade_signed messageShift )
{
    //================================================ Initialise local variables
    size_t origVolume                                 = xDimInds * yDimInds * zDimInds;
    size_t newVolume                                  = xDimIndsMsk * yDimIndsMsk * zDimIndsMsk;
    proshade_double* maskFinal;
    
    //================================================ If mask has different number of indices than map, then re-sample mask
    if ( ( xDimIndsMsk != xDimInds ) || ( yDimIndsMsk != yDimInds ) || ( zDimIndsMsk != zDimInds ) )
    {
        //============================================ Initialise variables
        fftw_complex* origCoeffs                      = new fftw_complex [newVolume ];
        fftw_complex* origCoeffsHKL                   = new fftw_complex [newVolume ];
        fftw_complex* modifCoeffs                     = new fftw_complex [origVolume];
        fftw_complex* modifCoeffsHKL                  = new fftw_complex [origVolume];
        fftw_complex* inMap                           = new fftw_complex [newVolume ];
        fftw_complex* outMap                          = new fftw_complex [origVolume];

        //============================================ Check memory allocation
        ProSHADE_internal_misc::checkMemoryAllocation ( inMap,          __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( outMap,         __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( origCoeffs,     __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( modifCoeffs,    __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( origCoeffsHKL,  __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( modifCoeffsHKL, __FILE__, __LINE__, __func__ );

        //============================================ Set array to zeroes
        for ( size_t iter = 0; iter < origVolume; iter++ ) { modifCoeffsHKL[iter][0] = 0.0; modifCoeffsHKL[iter][1] = 0.0; }

        //============================================ Cope mask to Fourier input array
        for ( size_t iter = 0; iter < newVolume; iter++ ) { inMap[iter][0] = mask[iter]; inMap[iter][1] = 0.0; }

        //============================================ Prepare Fourier transform plans
        fftw_plan planForwardFourier                  = fftw_plan_dft_3d ( static_cast< int > ( xDimIndsMsk ), static_cast< int > ( yDimIndsMsk ), static_cast< int > ( zDimIndsMsk ), inMap, origCoeffs, FFTW_FORWARD,  FFTW_ESTIMATE );
        fftw_plan inverseFoourier                     = fftw_plan_dft_3d ( static_cast< int > ( xDimInds ), static_cast< int > ( yDimInds ), static_cast< int > ( zDimInds ), modifCoeffs, outMap, FFTW_BACKWARD, FFTW_ESTIMATE );

        //============================================ Compute pre and post changes
        proshade_signed xPre, yPre, zPre;
        xPre                                          = std::abs ( ( static_cast< proshade_signed > ( xDimInds ) - static_cast< proshade_signed > ( xDimIndsMsk ) ) / 2 );
        yPre                                          = std::abs ( ( static_cast< proshade_signed > ( yDimInds ) - static_cast< proshade_signed > ( yDimIndsMsk ) ) / 2 );
        zPre                                          = std::abs ( ( static_cast< proshade_signed > ( zDimInds ) - static_cast< proshade_signed > ( zDimIndsMsk ) ) / 2 );
        
        if ( ( ( static_cast< proshade_signed > ( xDimInds ) - static_cast< proshade_signed > ( xDimIndsMsk ) ) % 2 ) == 1 ) { xPre -= 1; }
        if ( ( ( static_cast< proshade_signed > ( yDimInds ) - static_cast< proshade_signed > ( yDimIndsMsk ) ) % 2 ) == 1 ) { yPre -= 1; }
        if ( ( ( static_cast< proshade_signed > ( zDimInds ) - static_cast< proshade_signed > ( zDimIndsMsk ) ) % 2 ) == 1 ) { zPre -= 1; }

        //============================================ Run forward Fourier
        fftw_execute                                  ( planForwardFourier );

        //============================================ Initialise local variables
        proshade_signed maskMapIndex                  = 0;
        proshade_signed densMapIndex                  = 0;
        proshade_signed xMaskPos, yMaskPos, zMaskPos, xDensPos, yDensPos, zDensPos;
        proshade_signed maskH, maskK, maskL;

        //============================================ Convert mask to HKL for re-boxing
        for ( proshade_signed xIt = 0; xIt < static_cast< proshade_signed > ( xDimIndsMsk ); xIt++ )
        {
            for ( proshade_signed yIt = 0; yIt < static_cast< proshade_signed > ( yDimIndsMsk ); yIt++ )
            {
                for ( proshade_signed zIt = 0; zIt < static_cast< proshade_signed > ( zDimIndsMsk ); zIt++ )
                {
                    //================================ Convert to HKL
                    maskH                             = xIt + static_cast< proshade_signed > ( (xDimIndsMsk+1) / 2 ); if ( maskH >= static_cast< proshade_signed > ( xDimIndsMsk ) ) { maskH -= xDimIndsMsk; }
                    maskK                             = yIt + static_cast< proshade_signed > ( (yDimIndsMsk+1) / 2 ); if ( maskK >= static_cast< proshade_signed > ( yDimIndsMsk ) ) { maskK -= yDimIndsMsk; }
                    maskL                             = zIt + static_cast< proshade_signed > ( (zDimIndsMsk+1) / 2 ); if ( maskL >= static_cast< proshade_signed > ( zDimIndsMsk ) ) { maskL -= zDimIndsMsk; }

                    //================================ Find the positions
                    maskMapIndex                      = zIt   + static_cast< proshade_signed > ( zDimIndsMsk ) * ( yIt   + static_cast< proshade_signed > ( yDimIndsMsk ) * xIt   );
                    densMapIndex                      = maskL + static_cast< proshade_signed > ( zDimIndsMsk ) * ( maskK + static_cast< proshade_signed > ( yDimIndsMsk ) * maskH );

                    //================================ Save the values
                    origCoeffsHKL[densMapIndex][0]    = origCoeffs[maskMapIndex][0];
                    origCoeffsHKL[densMapIndex][1]    = origCoeffs[maskMapIndex][1];
                }
            }
        }

        //============================================ Rebox
        for ( proshade_signed xIt = 0; xIt < static_cast< proshade_signed > ( xDimInds ); xIt++ )
        {
            for ( proshade_signed yIt = 0; yIt < static_cast< proshade_signed > ( yDimInds ); yIt++ )
            {
                for ( proshade_signed zIt = 0; zIt < static_cast< proshade_signed > ( zDimInds ); zIt++ )
                {
                    //================================ Deal with X
                    if ( xDimIndsMsk >= xDimInds )    { xMaskPos = xIt + xPre; }
                    else                              { xMaskPos = xIt - xPre; }
                    xDensPos                          = xIt;

                    //================================ Deal with Y
                    if ( yDimIndsMsk >= yDimInds )    { yMaskPos = yIt + yPre; }
                    else                              { yMaskPos = yIt - yPre; }
                    yDensPos                          = yIt;

                    //================================ Deal with Z
                    if ( zDimIndsMsk >= zDimInds )    { zMaskPos = zIt + zPre; }
                    else                              { zMaskPos = zIt - zPre; }
                    zDensPos                          = zIt;

                    //================================ Skip if mask value not available (because the modifCoeffsHKL array is zeroed, we do not need to do anything here)
                    if ( ( xMaskPos < 0 ) || ( xMaskPos >= static_cast< proshade_signed > ( xDimIndsMsk ) ) ) { continue; }
                    if ( ( yMaskPos < 0 ) || ( yMaskPos >= static_cast< proshade_signed > ( yDimIndsMsk ) ) ) { continue; }
                    if ( ( zMaskPos < 0 ) || ( zMaskPos >= static_cast< proshade_signed > ( zDimIndsMsk ) ) ) { continue; }

                    //================================ Find the positions
                    maskMapIndex                      = zMaskPos + static_cast< proshade_signed > ( zDimIndsMsk ) * ( yMaskPos + static_cast< proshade_signed > ( yDimIndsMsk ) * xMaskPos );
                    densMapIndex                      = zDensPos + static_cast< proshade_signed > ( zDimInds    ) * ( yDensPos + static_cast< proshade_signed > ( yDimInds    ) * xDensPos );

                    //================================ Copy values
                    modifCoeffsHKL[densMapIndex][0]   = origCoeffsHKL[maskMapIndex][0];
                    modifCoeffsHKL[densMapIndex][1]   = origCoeffsHKL[maskMapIndex][1];
                }
            }
        }

        //============================================ Convert mask back to FFTW order
        for ( proshade_signed xIt = 0; xIt < static_cast< proshade_signed > ( xDimInds ); xIt++ )
        {
            for ( proshade_signed yIt = 0; yIt < static_cast< proshade_signed > ( yDimInds ); yIt++ )
            {
                for ( proshade_signed zIt = 0; zIt < static_cast< proshade_signed > ( zDimInds ); zIt++ )
                {
                    //================================ Convert to HKL
                    maskH                             = xIt + static_cast< proshade_signed > ( xDimInds / 2 ); if ( maskH >= static_cast< proshade_signed > ( xDimInds ) ) { maskH -= xDimInds; }
                    maskK                             = yIt + static_cast< proshade_signed > ( yDimInds / 2 ); if ( maskK >= static_cast< proshade_signed > ( yDimInds ) ) { maskK -= yDimInds; }
                    maskL                             = zIt + static_cast< proshade_signed > ( zDimInds / 2 ); if ( maskL >= static_cast< proshade_signed > ( zDimInds ) ) { maskL -= zDimInds; }
    
                    //================================ Find the positions
                    maskMapIndex                      = zIt   + static_cast< proshade_signed > ( zDimInds ) * ( yIt   + static_cast< proshade_signed > ( yDimInds ) * xIt );
                    densMapIndex                      = maskL + static_cast< proshade_signed > ( zDimInds ) * ( maskK + static_cast< proshade_signed > ( yDimInds ) * maskH );

                    //================================ Save the values
                    modifCoeffs[densMapIndex][0]      = modifCoeffsHKL[maskMapIndex][0];
                    modifCoeffs[densMapIndex][1]      = modifCoeffsHKL[maskMapIndex][1];
                }
            }
        }

        //============================================ Run inverse Fourier on the modified coefficients
        fftw_execute                                  ( inverseFoourier );

        //============================================ Delete old mask and allocate memory for the new, re-sampled mask
        maskFinal                                     = new proshade_double [origVolume];
        ProSHADE_internal_misc::checkMemoryAllocation ( maskFinal, __FILE__, __LINE__, __func__ );

        //======================================== Copy results into a new, properly sampled mask
        for ( size_t iter = 0; iter < origVolume; iter++ ) { maskFinal[iter] = outMap[iter][0]; }

        //======================================== Release remaining memory
        fftw_destroy_plan                         ( planForwardFourier );
        fftw_destroy_plan                         ( inverseFoourier );
        delete[] origCoeffs;
        delete[] modifCoeffs;
        delete[] origCoeffsHKL;
        delete[] modifCoeffsHKL;
        delete[] inMap;
        delete[] outMap;
    }
    else
    {
        maskFinal                                     = new proshade_double [origVolume];
        ProSHADE_internal_misc::checkMemoryAllocation ( maskFinal, __FILE__, __LINE__, __func__ );
        for ( size_t iter = 0; iter < origVolume; iter++ ) { maskFinal[iter] = mask[iter]; }
    }

    //================================================ Apply the mask to the map
    for ( size_t iter = 0; iter < static_cast< size_t > ( xDimInds * yDimInds * zDimInds ); iter++ ) { map[iter] *= maskFinal[iter]; }

    //================================================ Release memory
    delete[] maskFinal;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "Mask read in and applied successfully.", messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function reads and applies the Fourier weights to the map.
 
    This function reads in the weights file and determines, if its dimensions are the same as the input map or not. If not,
    then the weights are re-sampled is such a way as to have the same number of indices as the input map. Once the map
    and weights map dimensions are the same, the weights are applied by multiplying all map values by corresponding weight
    value. The map change is done in-place!
 
    \param[in] map Pointer to a the internal variable containing the map (it will be changed in place).
    \param[in] weightsFile Filename of where the mask is to be read from.
    \param[in] xDimInds The size of x dimension in indices.
    \param[in] yDimInds The size of y dimension in indices.
    \param[in] zDimInds The size of z dimension in indices.
    \param[in] verbose How much std::out output would you like?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
    \param[in] weightsArray An array of weights (default nullptr) to be used instead of input file.
    \param[in] waXInds The size of weightsArray x dimension in indices (defaults to 0).
    \param[in] waYInds The size of weightsArray y dimension in indices (defaults to 0).
    \param[in] waZInds The size of weightsArray z dimension in indices (defaults to 0).
 */
void ProSHADE_internal_io::applyWeights ( proshade_double*& map, std::string weightsFile, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds, proshade_signed verbose, proshade_signed messageShift, proshade_double* weightsArray, proshade_unsign waXInds, proshade_unsign waYInds, proshade_unsign waZInds )
{
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "Reading weights " << weightsFile;
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, hlpSS.str() );
    
    //================================================ Are we reading from file, or from array?
    if ( ( weightsArray != nullptr ) && ( waXInds != 0 ) && ( waYInds != 0 ) && ( waZInds != 0 ) )
    {
        //============================================ From array it is!
        ProSHADE_internal_io::applyWeightsFromArray   ( map, xDimInds, yDimInds, zDimInds, weightsArray, waXInds, waYInds, waZInds, verbose, messageShift );
    }
    else
    {
        //============================================ Check if weights file was given
        if ( weightsFile == "" )                      { ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "No weights supplied. Assuming all weights to be 1.0.", messageShift ); return; }
        
        //============================================ From file it is! Open the weights file
        gemmi::Ccp4<float> weights;
        weights.read_ccp4                             ( gemmi::MaybeGzipped ( weightsFile.c_str() ) );
        
        //============================================ Convert to XYZ and create complete weights, if need be
        weights.setup                                 ( gemmi::GridSetup::ReorderOnly, 0 );
        
        //============================================ Read in the rest of the weights file header
        proshade_unsign xDI, yDI, zDI, xAOR, yAOR, zAOR, xGI, yGI, zGI;
        proshade_single xDS, yDS, zDS, aA, bA, cA;
        proshade_signed xF, yF, zF, xAO, yAO, zAO;
        ProSHADE_internal_io::readInMapHeader         ( &weights,
                                                        &xDI,  &yDI,  &zDI,
                                                        &xDS,  &yDS,  &zDS,
                                                        &aA,   &bA,   &cA,
                                                        &xF,   &yF,   &zF,
                                                        &xAO,  &yAO,  &zAO,
                                                        &xAOR, &yAOR, &zAOR,
                                                        &xGI,  &yGI,  &zGI );

        //============================================ Save the weights values to local variable
        proshade_double* internalWeights              = nullptr;
        ProSHADE_internal_io::readInMapData           ( &weights, internalWeights, xDI, yDI, zDI, xAOR, yAOR, zAOR );
        
        //============================================ Apply weights from array
        ProSHADE_internal_io::applyWeightsFromArray   ( map, xDimInds, yDimInds, zDimInds, internalWeights, xDI, yDI, zDI, verbose, messageShift );

        //============================================ Release the memory
        delete[] internalWeights;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function applies the weights to the map.
 
    This function firstly determines if its dimensions of weights are the same as the input map or not. If not,
    then the weights are re-sampled is such a way as to have the same number of indices as the input map. Once the map
    and weights dimensions are the same, the weights are applied by multiplying all map Fourier coefficients by corresponding
    weights value. The map change is done in-place!
 
    \param[in] map Pointer to a the internal variable containing the map (it will be changed in place).
    \param[in] xDimInds The size of the map x dimension in indices.
    \param[in] yDimInds The size of the map y dimension in indices.
    \param[in] zDimInds The size of the map z dimension in indices.
    \param[in] weights A pointer to 1D array of the weights values.
    \param[in] xDimIndsWgh The size of the weights x dimension in indices.
    \param[in] yDimIndsWgh The size of the weights y dimension in indices.
    \param[in] zDimIndsWgh The size of the weights z dimension in indices.
    \param[in] verbose How much std::out output would you like?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
 */
void ProSHADE_internal_io::applyWeightsFromArray ( proshade_double*& map, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds, proshade_double*& weights, proshade_unsign xDimIndsWgh, proshade_unsign yDimIndsWgh, proshade_unsign zDimIndsWgh, proshade_signed verbose, proshade_signed messageShift )
{
    //================================================ Initialise local variables
    proshade_double* weightsFinal;
    size_t origVolume                                 = xDimInds * yDimInds * zDimInds;
    size_t newVolume                                  = xDimIndsWgh * yDimIndsWgh * zDimIndsWgh;
    
    //================================================ If weights have different number of indices than map, then re-sample weights in supplied space
    if ( ( xDimIndsWgh != xDimInds ) || ( yDimIndsWgh != yDimInds ) || ( zDimIndsWgh != zDimInds ) )
    {
        //============================================ Initialise variables
        fftw_complex* origCoeffs                      = new fftw_complex [newVolume ];
        fftw_complex* origCoeffsHKL                   = new fftw_complex [newVolume ];
        fftw_complex* modifCoeffs                     = new fftw_complex [origVolume];
        fftw_complex* modifCoeffsHKL                  = new fftw_complex [origVolume];
        fftw_complex* inMap                           = new fftw_complex [newVolume ];
        fftw_complex* outMap                          = new fftw_complex [origVolume];

        //============================================ Check memory allocation
        ProSHADE_internal_misc::checkMemoryAllocation ( inMap,          __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( outMap,         __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( origCoeffs,     __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( modifCoeffs,    __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( origCoeffsHKL,  __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( modifCoeffsHKL, __FILE__, __LINE__, __func__ );

        //============================================ Set array to zeroes
        for ( size_t iter = 0; iter < origVolume; iter++ ) { modifCoeffsHKL[iter][0] = 0.0; modifCoeffsHKL[iter][1] = 0.0; }

        //============================================ Copy weights to Fourier input array
        for ( size_t iter = 0; iter < newVolume; iter++ ) { inMap[iter][0] = weights[iter]; inMap[iter][1] = 0.0; }
        
        //============================================ Prepare Fourier transform plans
        fftw_plan planForwardFourier                  = fftw_plan_dft_3d ( static_cast< int > ( xDimIndsWgh ), static_cast< int > ( yDimIndsWgh ), static_cast< int > ( zDimIndsWgh ), inMap, origCoeffs, FFTW_FORWARD,  FFTW_ESTIMATE );
        fftw_plan inverseFoourier                     = fftw_plan_dft_3d ( static_cast< int > ( xDimInds ), static_cast< int > ( yDimInds ), static_cast< int > ( zDimInds ), modifCoeffs, outMap, FFTW_BACKWARD, FFTW_ESTIMATE );

        //============================================ Compute pre and post changes
        proshade_signed xPre, yPre, zPre;
        xPre                                          = std::abs ( ( static_cast< proshade_signed > ( xDimInds ) - static_cast< proshade_signed > ( xDimIndsWgh ) ) / 2 );
        yPre                                          = std::abs ( ( static_cast< proshade_signed > ( yDimInds ) - static_cast< proshade_signed > ( yDimIndsWgh ) ) / 2 );
        zPre                                          = std::abs ( ( static_cast< proshade_signed > ( zDimInds ) - static_cast< proshade_signed > ( zDimIndsWgh ) ) / 2 );
        
        if ( ( ( static_cast< proshade_signed > ( xDimInds ) - static_cast< proshade_signed > ( xDimIndsWgh ) ) % 2 ) == 1 ) { xPre -= 1; }
        if ( ( ( static_cast< proshade_signed > ( yDimInds ) - static_cast< proshade_signed > ( yDimIndsWgh ) ) % 2 ) == 1 ) { yPre -= 1; }
        if ( ( ( static_cast< proshade_signed > ( zDimInds ) - static_cast< proshade_signed > ( zDimIndsWgh ) ) % 2 ) == 1 ) { zPre -= 1; }

        //============================================ Run forward Fourier
        fftw_execute                                  ( planForwardFourier );

        //============================================ Initialise local variables
        proshade_signed maskMapIndex                  = 0;
        proshade_signed densMapIndex                  = 0;
        proshade_signed xMaskPos, yMaskPos, zMaskPos, xDensPos, yDensPos, zDensPos;
        proshade_signed maskH, maskK, maskL;

        //============================================ Convert weights to HKL for re-boxing
        for ( proshade_signed xIt = 0; xIt < static_cast< proshade_signed > ( xDimIndsWgh ); xIt++ )
        {
            for ( proshade_signed yIt = 0; yIt < static_cast< proshade_signed > ( yDimIndsWgh ); yIt++ )
            {
                for ( proshade_signed zIt = 0; zIt < static_cast< proshade_signed > ( zDimIndsWgh ); zIt++ )
                {
                    //================================ Convert to HKL
                    maskH                             = xIt + static_cast< proshade_signed > ( (xDimIndsWgh+1) / 2 ); if ( maskH >= static_cast< proshade_signed > ( xDimIndsWgh ) ) { maskH -= xDimIndsWgh; }
                    maskK                             = yIt + static_cast< proshade_signed > ( (yDimIndsWgh+1) / 2 ); if ( maskK >= static_cast< proshade_signed > ( yDimIndsWgh ) ) { maskK -= yDimIndsWgh; }
                    maskL                             = zIt + static_cast< proshade_signed > ( (zDimIndsWgh+1) / 2 ); if ( maskL >= static_cast< proshade_signed > ( zDimIndsWgh ) ) { maskL -= zDimIndsWgh; }

                    //================================ Find the positions
                    maskMapIndex                      = zIt   + static_cast< proshade_signed > ( zDimIndsWgh ) * ( yIt   + static_cast< proshade_signed > ( yDimIndsWgh ) * xIt   );
                    densMapIndex                      = maskL + static_cast< proshade_signed > ( zDimIndsWgh ) * ( maskK + static_cast< proshade_signed > ( yDimIndsWgh ) * maskH );

                    //================================ Save the values
                    origCoeffsHKL[densMapIndex][0]    = origCoeffs[maskMapIndex][0];
                    origCoeffsHKL[densMapIndex][1]    = origCoeffs[maskMapIndex][1];
                }
            }
        }

        //============================================ Rebox
        for ( proshade_signed xIt = 0; xIt < static_cast< proshade_signed > ( xDimInds ); xIt++ )
        {
            for ( proshade_signed yIt = 0; yIt < static_cast< proshade_signed > ( yDimInds ); yIt++ )
            {
                for ( proshade_signed zIt = 0; zIt < static_cast< proshade_signed > ( zDimInds ); zIt++ )
                {
                    //================================ Deal with X
                    if ( xDimIndsWgh >= xDimInds )    { xMaskPos = xIt + xPre; }
                    else                              { xMaskPos = xIt - xPre; }
                    xDensPos                          = xIt;

                    //================================ Deal with Y
                    if ( yDimIndsWgh >= yDimInds )    { yMaskPos = yIt + yPre; }
                    else                              { yMaskPos = yIt - yPre; }
                    yDensPos                          = yIt;

                    //================================ Deal with Z
                    if ( zDimIndsWgh >= zDimInds )    { zMaskPos = zIt + zPre; }
                    else                              { zMaskPos = zIt - zPre; }
                    zDensPos                          = zIt;

                    //================================ Skip if weights value not available (because the modifCoeffsHKL array is zeroed, we do not need to do anything here)
                    if ( ( xMaskPos < 0 ) || ( xMaskPos >= static_cast< proshade_signed > ( xDimIndsWgh ) ) ) { continue; }
                    if ( ( yMaskPos < 0 ) || ( yMaskPos >= static_cast< proshade_signed > ( yDimIndsWgh ) ) ) { continue; }
                    if ( ( zMaskPos < 0 ) || ( zMaskPos >= static_cast< proshade_signed > ( zDimIndsWgh ) ) ) { continue; }

                    //================================ Find the positions
                    maskMapIndex                      = zMaskPos + static_cast< proshade_signed > ( zDimIndsWgh ) * ( yMaskPos + static_cast< proshade_signed > ( yDimIndsWgh ) * xMaskPos );
                    densMapIndex                      = zDensPos + static_cast< proshade_signed > ( zDimInds    ) * ( yDensPos + static_cast< proshade_signed > ( yDimInds    ) * xDensPos );

                    //================================ Copy values
                    modifCoeffsHKL[densMapIndex][0]   = origCoeffsHKL[maskMapIndex][0];
                    modifCoeffsHKL[densMapIndex][1]   = origCoeffsHKL[maskMapIndex][1];
                }
            }
        }

        //============================================ Convert weights back to FFTW order
        for ( proshade_signed xIt = 0; xIt < static_cast< proshade_signed > ( xDimInds ); xIt++ )
        {
            for ( proshade_signed yIt = 0; yIt < static_cast< proshade_signed > ( yDimInds ); yIt++ )
            {
                for ( proshade_signed zIt = 0; zIt < static_cast< proshade_signed > ( zDimInds ); zIt++ )
                {
                    //================================ Convert to HKL
                    maskH                             = xIt + static_cast< proshade_signed > ( xDimInds / 2 ); if ( maskH >= static_cast< proshade_signed > ( xDimInds ) ) { maskH -= xDimInds; }
                    maskK                             = yIt + static_cast< proshade_signed > ( yDimInds / 2 ); if ( maskK >= static_cast< proshade_signed > ( yDimInds ) ) { maskK -= yDimInds; }
                    maskL                             = zIt + static_cast< proshade_signed > ( zDimInds / 2 ); if ( maskL >= static_cast< proshade_signed > ( zDimInds ) ) { maskL -= zDimInds; }
    
                    //================================ Find the positions
                    maskMapIndex                      = zIt   + static_cast< proshade_signed > ( zDimInds ) * ( yIt   + static_cast< proshade_signed > ( yDimInds ) * xIt );
                    densMapIndex                      = maskL + static_cast< proshade_signed > ( zDimInds ) * ( maskK + static_cast< proshade_signed > ( yDimInds ) * maskH );

                    //================================ Save the values
                    modifCoeffs[densMapIndex][0]      = modifCoeffsHKL[maskMapIndex][0];
                    modifCoeffs[densMapIndex][1]      = modifCoeffsHKL[maskMapIndex][1];
                }
            }
        }

        //============================================ Run inverse Fourier on the modified coefficients
        fftw_execute                                  ( inverseFoourier );

        //============================================ Delete old weights and allocate memory for the new, re-sampled weights
        weightsFinal                                  = new proshade_double [origVolume];
        ProSHADE_internal_misc::checkMemoryAllocation ( weightsFinal, __FILE__, __LINE__, __func__ );

        //============================================ Copy results into a new, properly sampled weights
        for ( size_t iter = 0; iter < origVolume; iter++ ) { weightsFinal[iter] = outMap[iter][0]; }

        //============================================ Release remaining memory
        fftw_destroy_plan                             ( planForwardFourier );
        fftw_destroy_plan                             ( inverseFoourier );
        delete[] origCoeffs;
        delete[] modifCoeffs;
        delete[] origCoeffsHKL;
        delete[] modifCoeffsHKL;
        delete[] inMap;
        delete[] outMap;
    }
    else
    {
        weightsFinal                                  = new proshade_double [origVolume];
        ProSHADE_internal_misc::checkMemoryAllocation ( weightsFinal, __FILE__, __LINE__, __func__ );
        for ( size_t iter = 0; iter < origVolume; iter++ ) { weightsFinal[iter] = weights[iter]; }
    }
    
    //================================================ Allocate memory for map Fourier transform
    fftw_complex* inMap                               = new fftw_complex [origVolume];
    fftw_complex* outMap                              = new fftw_complex [origVolume];
    ProSHADE_internal_misc::checkMemoryAllocation     ( inMap,  __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( outMap, __FILE__, __LINE__, __func__ );
    fftw_plan planForwardFourier                      = fftw_plan_dft_3d ( static_cast< int > ( xDimInds ), static_cast< int > ( yDimInds ), static_cast< int > ( zDimInds ), inMap, outMap, FFTW_FORWARD,  FFTW_ESTIMATE );
    fftw_plan inverseFoourier                         = fftw_plan_dft_3d ( static_cast< int > ( xDimInds ), static_cast< int > ( yDimInds ), static_cast< int > ( zDimInds ), outMap, inMap, FFTW_BACKWARD, FFTW_ESTIMATE );

    //================================================ Set data
    for ( size_t iter = 0; iter < static_cast< size_t > ( origVolume ); iter++ ) { inMap[iter][0] = map[iter]; inMap[iter][1] = 0.0; }
    
    //================================================ Convert map to Fourier space
    fftw_execute                                      ( planForwardFourier );
    
    //================================================ Apply the weights to the map in Fourier space
    proshade_double normFactor                        = static_cast<proshade_double> ( origVolume );
    for ( size_t iter = 0; iter < static_cast< size_t > ( origVolume ); iter++ ) { outMap[iter][0] *= weightsFinal[iter] / normFactor; outMap[iter][1] *= weightsFinal[iter] / normFactor; }
    
    //================================================ Convert weighted map from Fourier space
    fftw_execute                                      ( inverseFoourier );

    //================================================ Copy results to map
    for ( size_t iter = 0; iter < static_cast< size_t > ( xDimInds * yDimInds * zDimInds ); iter++ ) { map[iter] = inMap[iter][0]; }
    
    //================================================ Release memory
    delete[] weightsFinal;
    delete[] inMap;
    delete[] outMap;
    fftw_destroy_plan                                 ( planForwardFourier );
    fftw_destroy_plan                                 ( inverseFoourier );
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "Mask read in and applied successfully.", messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function parses the CCP4 MAP file header as read in by gemmi.
 
    This function uses the gemmi Ccp4 object, which contains all the information read in from a MAP file (including the header), to parse out the ProSHADE required
    information from the header and saving it to the supplied variables.
 
    \param[in] map A gemmi Ccp4 objecct containing all the data read in from a MAP file.
    \param[in] xDimInds Variable holding the x-axis size in indices.
    \param[in] yDimInds Variable holding the y-axis size in indices.
    \param[in] zDimInds Variable holding the z-axis size in indices.
    \param[in] xDim Variable holding the x dimension size in angstroms.
    \param[in] yDim Variable holding the y dimension size in angstroms.
    \param[in] zDim Variable holding the z dimension size in angstroms.
    \param[in] aAng Variable holding the a angle in degrees.
    \param[in] bAng Variable holding the b angle in degrees.
    \param[in] cAng Variable holding the c angle in degrees.
    \param[in] xFrom Variable holding the starting index along the x-axis.
    \param[in] yFrom Variable holding the starting index along the y-axis.
    \param[in] zFrom Variable holding the starting index along the z-axis.
    \param[in] xAxOrigin Variable holding the map origin positon along the x-axis.
    \param[in] yAxOrigin Variable holding the map origin positon along the y-axis.
    \param[in] zAxOrigin Variable holding the map origin positon along the z-axis.
    \param[in] xAxOrder Variable holding the order of x axis.
    \param[in] yAxOrder Variable holding the order of y axis.
    \param[in] zAxOrder Variable holding the order of z axis.
    \param[in] xGridInds Variable holding the grid indices count along the x-axis.
    \param[in] yGridInds Variable holding the grid indices count along the y-axis.
    \param[in] zGridInds Variable holding the grid indices count along the z-axis.
    \param[in] title The title to be written into the MAP file.
    \param[in] mode The variable type of the data, please leave two (float) unless you require any specific other mode.
 */
void ProSHADE_internal_io::writeOutMapHeader ( gemmi::Ccp4<float> *map, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds, proshade_single xDim, proshade_single yDim, proshade_single zDim, proshade_single aAng, proshade_single bAng, proshade_single cAng, proshade_signed xFrom, proshade_signed yFrom, proshade_signed zFrom, proshade_signed xAxOrigin, proshade_signed yAxOrigin, proshade_signed zAxOrigin, proshade_unsign xAxOrder, proshade_unsign yAxOrder, proshade_unsign zAxOrder, proshade_unsign xGridInds, proshade_unsign yGridInds, proshade_unsign zGridInds, std::string title, int mode )
{
    //================================================ Fill in the map file header
    map->set_header_i32                               ( 1 , static_cast<int32_t> ( xDimInds ) );       // Number of columns in 3D data array (fast axis)
    map->set_header_i32                               ( 2 , static_cast<int32_t> ( yDimInds ) );       // Number of columns in 3D data array (medium axis)
    map->set_header_i32                               ( 3 , static_cast<int32_t> ( zDimInds ) );       // Number of columns in 3D data array (slow axis)
    map->set_header_i32                               ( 4 , static_cast<int32_t> ( mode ) );           // Map mode
    map->set_header_i32                               ( 5 , static_cast<int32_t> ( xFrom ) );          // Starting index (fast axis)
    map->set_header_i32                               ( 6 , static_cast<int32_t> ( yFrom ) );          // Starting index (medium axis)
    map->set_header_i32                               ( 7 , static_cast<int32_t> ( zFrom ) );          // Starting index (slow axis)
    map->set_header_i32                               ( 8 , static_cast<int32_t> ( xGridInds ) );      // Grid sampling (fast axis)
    map->set_header_i32                               ( 9 , static_cast<int32_t> ( yGridInds ) );      // Grid sampling (medium axis)
    map->set_header_i32                               ( 10, static_cast<int32_t> ( zGridInds ) );      // Grid sampling (slow axis)
    map->set_header_float                             ( 11, static_cast<float>   ( xDim ) );           // Grid dimension in Angstrom (fast axis)
    map->set_header_float                             ( 12, static_cast<float>   ( yDim ) );           // Grid dimension in Angstrom (medium axis)
    map->set_header_float                             ( 13, static_cast<float>   ( zDim ) );           // Grid dimension in Angstrom (slow axis)
    map->set_header_float                             ( 14, static_cast<float>   ( aAng ) );           // Alpha angle in degrees
    map->set_header_float                             ( 15, static_cast<float>   ( bAng ) );           // Beta angle in degrees
    map->set_header_float                             ( 16, static_cast<float>   ( cAng ) );           // Gamma angle in degrees
    map->set_header_i32                               ( 17, static_cast<int32_t> ( xAxOrder ) );       // MAPC
    map->set_header_i32                               ( 18, static_cast<int32_t> ( yAxOrder ) );       // MAPR
    map->set_header_i32                               ( 19, static_cast<int32_t> ( zAxOrder ) );       // MAPS
    if ( map->grid.spacegroup ) { map->set_header_i32 ( 23, static_cast<int32_t> ( map->grid.spacegroup->ccp4 ) ); } // Space group
    else                        { map->set_header_i32 ( 23, static_cast<int32_t> ( 1 ) ); }
    map->set_header_i32                               ( 24, static_cast<int32_t> ( map->grid.spacegroup->operations().order() * 80 ) ); // NSYMBT - size of extended header (which follows main header) in bytes
    map->set_header_str                               ( 27, "CCP4" );                                  // Code for the type of extended header
    map->set_header_i32                               ( 28, static_cast<int32_t> ( 20140 ) );          // Version
    map->set_header_i32                               ( 50, static_cast<int32_t> ( xAxOrigin ) );      // Origin of the map (fast axis)
    map->set_header_i32                               ( 51, static_cast<int32_t> ( yAxOrigin ) );      // Origin of the map (medium axis)
    map->set_header_i32                               ( 52, static_cast<int32_t> ( zAxOrigin ) );      // Origin of the map (slow axis)
    map->set_header_str                               ( 53, "MAP" );                                   // File format
    if ( gemmi::is_little_endian() ) { map->set_header_i32 ( 54, static_cast<int32_t> ( 0x00004144 ) ); }       // Machine stamp encoding byte ordering of data
    else                             { map->set_header_i32 ( 54, static_cast<int32_t> ( 0x11110000 ) ); }
    map->set_header_i32                               ( 56, static_cast<int32_t> ( 1 ) );                       // Number of labels used
    std::memset                                       ( reinterpret_cast<void*> ( &(map->ccp4_header.at( 56 )) ), ' ', static_cast< size_t > ( 800 + map->grid.spacegroup->operations().order() * 80 ) ); // 56 is used because the vector is indexed from 0
    map->set_header_str                               ( 57, title );                                            // Title
    
    //================================================ Done
    return ;
    
}

/*! \brief Function determining input data type.
 
    This function determines the type of the input structure. The possible outputs are MAP for MRC map files, PDB for mmCIF or PDB formatted data,
    or UNKNOWN if gemmi fail to read the file as co-ordinates as well as map.
 
    \param[in] fName The file name of the file for which the type should be determined.
    \param[out] X ProSHADE InputType variable with values UNKNOWN, MAP or PDB depending on the type of the input file.
 */
ProSHADE_internal_io::InputType ProSHADE_internal_io::figureDataType ( std::string fName )
{
    //================================================ Try readin as PDB
    if ( isFilePDB ( fName ) )
    {
        return                                        ( PDB );
    }
    
    //================================================ If not, try readin as MAP
    if ( isFileMAP ( fName ) )
    {
        return                                        ( MAP );
    }
    
    //================================================ No luck? UNKNOWN it is ...
    return                                            ( UNKNOWN );
    
    //================================================ Done
    
}

/*! \brief Function for writing out the optimal rotation and translation into a JSON file.
 
    This function takes the initial translation (assuming that the centre of rotation is not at the centre of indices), the
    rotation (in terms of the Euler angles) and the translation detected by the overlay task and proceeds to
    write a JSON file containing all of these information in the order in which they should be applied with the exception
    that depending on around which point the rotation should be done, either the translation to origin or to the map centre
    should be applied.
 
    This function assumes that the second translation includes the reverse of the first translation already.
 
    \param[in] trsX1 The translation required to get the rotation centre to origin along the x-axis.
    \param[in] trsY1 The translation required to get the rotation centre to origin along the y-axis.
    \param[in] trsZ1 The translation required to get the rotation centre to origin along the z-axis.
    \param[in] eulA The optimal rotation Euler angle alpha.
    \param[in] eulB The optimal rotation Euler angle beta.
    \param[in] eulG The optimal rotation Euler angle gamma.
    \param[in] trsX2 The optimal translation along the x-axis + reverse of the trsX1.
    \param[in] trsY2 The optimal translation along the y-axis + reverse of the trsY1.
    \param[in] trsZ2 The optimal translation along the z-axis + reverse of the trsZ1.
    \param[in] fileName The file name of the file for which the information should be written into.
 */
void ProSHADE_internal_io::writeRotationTranslationJSON ( proshade_double trsX1, proshade_double trsY1, proshade_double trsZ1, proshade_double eulA, proshade_double eulB, proshade_double eulG, proshade_double trsX2, proshade_double trsY2, proshade_double trsZ2, std::string fileName )
{
    //================================================ Open file for writing
    std::ofstream jsonFile;
    jsonFile.open                                     ( fileName );
    
    //================================================ Check file opening success
    if ( !jsonFile.is_open( ) )
    {
        throw ProSHADE_exception ( "Failed to open JSON output file.", "E000056", __FILE__, __LINE__, __func__, "Failed to open json file to which the rotation and\n                    : translation would be written into. Most likely cause is\n                    : lack of rights to write in the current folder." );
    }
    
    //================================================ Get rotation matrix from Euler angles
    proshade_double* rotMat                           = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( eulA, eulB, eulG, rotMat );
    
    //================================================ Write the info
    jsonFile << "{\n";
    jsonFile << "   \"translationToOrigin\" :           [ " << trsX1 << ", " << trsY1 << ", " << trsZ1 << " ], \n";
    
    jsonFile << "   \"rotationMatrix:\" :               [ " << rotMat[0] << ", " << rotMat[1] << ", " << rotMat[2] << ", \n";
    jsonFile << "                                         " << rotMat[3] << ", " << rotMat[4] << ", " << rotMat[5] << ", \n";
    jsonFile << "                                         " << rotMat[6] << ", " << rotMat[7] << ", " << rotMat[8] << "], \n";

    jsonFile << "   \"translationFromRotCenToOverlay\" : [ " << trsX2 << ", " << trsY2 << ", " << trsZ2 << " ] \n";
    jsonFile << "}\n";
    
    //================================================ Close file
    jsonFile.close                                    ( );
    
    //================================================ Release memory
    delete[] rotMat;
    
    //================================================ Done
    return ;
    
}
