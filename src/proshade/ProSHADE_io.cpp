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
    \version   0.7.3
    \date      JUL 2020
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
   
   *xDim                                              = static_cast<proshade_single> ( map->header_float ( 11 ) );
   *yDim                                              = static_cast<proshade_single> ( map->header_float ( 12 ) );
   *zDim                                              = static_cast<proshade_single> ( map->header_float ( 13 ) );
   
   *aAng                                              = static_cast<proshade_single> ( map->header_float ( 14 ) );
   *bAng                                              = static_cast<proshade_single> ( map->header_float ( 15 ) );
   *cAng                                              = static_cast<proshade_single> ( map->header_float ( 16 ) );
   
   *xFrom                                             = static_cast<proshade_signed> ( map->header_i32   ( 5  ) );
   *yFrom                                             = static_cast<proshade_signed> ( map->header_i32   ( 6  ) );
   *zFrom                                             = static_cast<proshade_signed> ( map->header_i32   ( 7  ) );
   
   *xAxOrigin                                         = static_cast<proshade_signed> ( map->header_i32   ( 50 ) ) + (*xFrom);
   *yAxOrigin                                         = static_cast<proshade_signed> ( map->header_i32   ( 51 ) ) + (*yFrom);
   *zAxOrigin                                         = static_cast<proshade_signed> ( map->header_i32   ( 52 ) ) + (*zFrom);
   
   *xAxOrder                                          = static_cast<proshade_unsign> ( map->header_i32   ( 17 ) );
   *yAxOrder                                          = static_cast<proshade_unsign> ( map->header_i32   ( 18 ) );
   *zAxOrder                                          = static_cast<proshade_unsign> ( map->header_i32   ( 19 ) );
   
   *xGridInds                                         = static_cast<proshade_unsign> ( *xDimInds );
   *yGridInds                                         = static_cast<proshade_unsign> ( *yDimInds );
   *zGridInds                                         = static_cast<proshade_unsign> ( *zDimInds );
    
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
                map[arrPos]                           = gemmiMap->grid.get_value_q( axOrdArr[xAxOrder-1], axOrdArr[yAxOrder-1], axOrdArr[zAxOrder-1] );
            }
        }
    }
    
    //================================================ Release internal variables memory
    delete[] axDimArr;
    delete[] axOrdArr;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for writing out the map cell into the header.
 
 This function takes all the required information and writes the map header with the cell information.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to write to.
 \param[in] xAngs How many angstroms are there along the x dimension.
 \param[in] yAngs How many angstroms are there along the y dimension
 \param[in] zAngs How many angstroms are there along the z dimension.
 \param[in] aAngle The cell a angle in degrees.
 \param[in] bAngle The cell b angle in degrees.
 \param[in] cAngle The cell c angle in degrees.
 */
void ProSHADE_internal_io::writeMapCell ( CMap_io::CMMFile* mapFile, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_single aAngle, proshade_single bAngle, proshade_single cAngle )
{
    //================================================ Declare local variables
    float cell[6];

    //================================================ Initialise writen veriable to correct values
    cell[0]                                           = static_cast<float> ( xAngs );
    cell[1]                                           = static_cast<float> ( yAngs );
    cell[2]                                           = static_cast<float> ( zAngs );
    cell[3]                                           = static_cast<float> ( aAngle );
    cell[4]                                           = static_cast<float> ( bAngle );
    cell[5]                                           = static_cast<float> ( cAngle );
    
    //================================================ Write data to header
    CMap_io::ccp4_cmap_set_cell                       ( mapFile, cell );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for writing out the map grid into the header.
 
 This function takes all the required information and writes the map header with the grid information.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to write to.
 \param[in] xGrid How many indices are there along the x dimension of the grid.
 \param[in] yGrid How many indices are there along the y dimension of the grid.
 \param[in] zGrid How many indices are there along the z dimension of the grid.
 */
void ProSHADE_internal_io::writeMapGrid ( CMap_io::CMMFile* mapFile, proshade_unsign xGrid, proshade_unsign yGrid, proshade_unsign zGrid )
{
    //================================================ Declare local variables
    int grid[3];
    
    //================================================ Initialise writen veriable to correct values
    grid[0]                                           = static_cast<int> ( xGrid );
    grid[1]                                           = static_cast<int> ( yGrid );
    grid[2]                                           = static_cast<int> ( zGrid );
    
    //================================================ Write data to header
    CMap_io::ccp4_cmap_set_grid                       ( mapFile, grid );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for writing out the map axis order into the header.
 
 This function takes all the required information and writes the map header with the axis order information.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to write to.
 \param[in] xAxisOrder The order of the x axis.
 \param[in] yAxisOrder The order of the y axis.
 \param[in] zAxisOrder The order of the z axis.
 */
void ProSHADE_internal_io::writeMapOrder ( CMap_io::CMMFile* mapFile, proshade_unsign xAxisOrder, proshade_unsign yAxisOrder, proshade_unsign zAxisOrder )
{
    //================================================ Declare local variables
    int order[3];
    
    //================================================ Initialise writen veriable to correct values
    order[0]                                          = static_cast<int> ( xAxisOrder );
    order[1]                                          = static_cast<int> ( yAxisOrder );
    order[2]                                          = static_cast<int> ( zAxisOrder );
    
    //================================================ Write data to header
    CMap_io::ccp4_cmap_set_order                      ( mapFile, order );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for writing out the map dimensions into the header.
 
 This function takes all the required information and writes the map header with the dimension sizes information.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to write to.
 \param[in] xDims The number of indiced alon the x axis in the map.
 \param[in] yDims The number of indiced alon the y axis in the map.
 \param[in] zDims The number of indiced alon the z axis in the map.
 */
void ProSHADE_internal_io::writeMapDims ( CMap_io::CMMFile* mapFile, proshade_unsign xDims, proshade_unsign yDims, proshade_unsign zDims )
{
    //================================================ Declare local variables
    int dims[3];
    
    //================================================ Initialise writen veriable to correct values
    dims[0]                                           = static_cast<int> ( xDims );
    dims[1]                                           = static_cast<int> ( yDims );
    dims[2]                                           = static_cast<int> ( zDims );
    
    //================================================ Write data to header
    CMap_io::ccp4_cmap_set_dim                ( mapFile, dims );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for writing out the map origin into the header.
 
 This function takes all the required information and writes the map header with the origin location information.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to write to.
 \param[in] xOrigin The position of the first index along the x axis.
 \param[in] yOrigin The position of the first index along the y axis.
 \param[in] zOrigin The position of the first index along the z axis.
 */
void ProSHADE_internal_io::writeMapOrigin ( CMap_io::CMMFile* mapFile, proshade_unsign xOrigin, proshade_unsign yOrigin, proshade_unsign zOrigin )
{
    //================================================ Declare local variables
    int orig[3];
    
    //================================================ Initialise writen veriable to correct values
    orig[0]                                           = static_cast<int> ( xOrigin );
    orig[1]                                           = static_cast<int> ( yOrigin );
    orig[2]                                           = static_cast<int> ( zOrigin );
    
    //================================================ Write data to header
    CMap_io::ccp4_cmap_set_origin                     ( mapFile, orig );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for writing out the spacegroup, map mode and title into the header.
 
 This function takes all the required information and writes the map header with the spacegroup, map mode and title information.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to write to.
 \param[in] title String which should become the title of the saved map.
 \param[in] mode The map mode to be used. Defaul value is 2.
 \param[in] spaceGroup The spacegroup value to be saved. If this value is not 1, I have no clue what will happen and whether all will fail. Just do not change it.
 */
void ProSHADE_internal_io::writeMapTitleEtc ( CMap_io::CMMFile* mapFile, std::string title, proshade_unsign mode, proshade_unsign spaceGroup )
{
    //================================================ Write obvious data to header
    CMap_io::ccp4_cmap_set_spacegroup                 ( mapFile, static_cast<int> ( spaceGroup ) );
    CMap_io::ccp4_cmap_set_datamode                   ( mapFile, static_cast<int> ( mode ) );
    
    //================================================ Deal with the title
    if ( title == "" )
    {
        const char* titl                              = "ProSHADE genrated map                                                           ";
        CMap_io::ccp4_cmap_set_title                  ( mapFile, titl );
    }
    else
    {
        char buff[80];
        snprintf                                      ( buff, sizeof ( buff ), "%s", title.c_str() );
        CMap_io::ccp4_cmap_set_title                  ( mapFile, buff );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for writing out the map into the data portion of the output file.
 
 This function writes the internal map representation into a MRC MAP file.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to write to.
 \param[in] map The map to be written into the data portion of the file.
 \param[in] xDim The number of indices along the x dimension.
 \param[in] yDim The number of indices along the y dimension.
 \param[in] zDim The number of indices along the z dimension.
 \param[in] xAxisOrder The order of the x axis.
 \param[in] yAxisOrder The order of the y axis.
 \param[in] zAxisOrder The order of the z axis.
 */
void ProSHADE_internal_io::writeMapData ( CMap_io::CMMFile* mapFile, proshade_double* map, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_unsign xAxisOrder, proshade_unsign yAxisOrder, proshade_unsign zAxisOrder )
{
    
    //================================================ Initialise local variables
    proshade_unsign dims[3];      dims[0]      = xDim;       dims[1]      = yDim;       dims[2]      = zDim;
    proshade_unsign axisOrder[3]; axisOrder[0] = xAxisOrder; axisOrder[1] = yAxisOrder; axisOrder[2] = zAxisOrder;
    std::vector<float> section                ( dims[axisOrder[0]-1] * dims[axisOrder[1]-1] );
    proshade_unsign index;
    proshade_unsign iters[3];
    proshade_unsign arrPos;
    
    //================================================ Write out the map data
    for ( iters[2] = 0; iters[2] < dims[axisOrder[2]-1]; iters[2]++ )
    {
        index = 0;
        
        for ( iters[1] = 0; iters[1] < dims[axisOrder[1]-1]; iters[1]++ )
        {
            for ( iters[0] = 0; iters[0] < dims[axisOrder[0]-1]; iters[0]++ )
            {
                arrPos                                = iters[2]  + (dims[axisOrder[2]-1]) * ( iters[1]  + (dims[axisOrder[1]-1]) * iters[0] );
                section[ index++ ]                    = static_cast<float> ( map[arrPos] );
            }
        }
        
        CMap_io::ccp4_cmap_write_section              ( mapFile, &section[0] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Function determining input data type.
 
 This function determines the type of the input structure. The possible outputs are MAP for MRC map files, PDB for mmCIF or PDB formatted data,
 or UNKNOWN if both cmaplib and gemmi fail to read the file.
 
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
