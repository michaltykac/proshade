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
    //================================================ Create CCP4MAPfile object
    CMap_io::CMMFile *mapFile                         = NULL;
    
    //================================================ Read in the file, if possible
    CCP4::ccp4_liberr_verbosity                       ( 0 );
    mapFile                                           = reinterpret_cast<CMap_io::CMMFile*> ( CMap_io::ccp4_cmap_open ( fName.c_str() , O_RDONLY ) );
    
    //================================================ Catch for fail
    if ( mapFile == NULL )
    {
        //============================================ The file is not MAP format
        CMap_io::ccp4_cmap_close                      ( mapFile );
        
        //============================================ Done testing
        return                                        ( false );
    }
    //================================================ The file is MAP format
    CMap_io::ccp4_cmap_close                          ( mapFile );
    return                                            ( true );
    
    //================================================ Done
    
}

/*! \brief Function or reading the map cell from the header.
 
 This function takes an open map file and reads the cell information, saving the read values to the supp;ied
 variables.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to read from.
 \param[in] xDim Address to a variable to save the x dimension size in angstroms.
 \param[in] yDim Address to a variable to save the y dimension size in angstroms.
 \param[in] zDim Address to a variable to save the z dimension size in angstroms.
 \param[in] aAng Address to a variable to save the a angle in degrees.
 \param[in] bAng Address to a variable to save the b angle in degrees.
 \param[in] cAng Address to a variable to save the c angle in degrees.
 */
void ProSHADE_internal_io::readInMapCell ( CMap_io::CMMFile* mapFile, proshade_single* xDim, proshade_single* yDim, proshade_single* zDim, proshade_single* aAng, proshade_single* bAng, proshade_single* cAng )
{
    //================================================ Initialise variables
    float *cell                                       = NULL;
    
    //================================================ Allocate memory
    cell                                              = (float*) malloc ( 6 * sizeof ( float ) );
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( cell, __FILE__, __LINE__, __func__ );
    
    //================================================ Read the data
    CMap_io::ccp4_cmap_get_cell                       ( mapFile, cell );
    
    //================================================ Save the data
   *xDim                                              = static_cast<proshade_single> ( cell[0] );
   *yDim                                              = static_cast<proshade_single> ( cell[1] );
   *zDim                                              = static_cast<proshade_single> ( cell[2] );
   *aAng                                              = static_cast<proshade_single> ( cell[3] );
   *bAng                                              = static_cast<proshade_single> ( cell[4] );
   *cAng                                              = static_cast<proshade_single> ( cell[5] );
    
    //================================================ Check for perpendicular axes - only P1 is supported for now.
    if ( ( *aAng != 90.0 ) || ( *bAng != 90.0 ) || ( *cAng != 90.0 ) )
    {
        throw ProSHADE_exception ( "The map axes are not perpendicular. Only P1 cells are\n                    : supported for now.", "EM00046", __FILE__, __LINE__, __func__, "ProSHADE currently only supports map cells with\n                    : perpendicular (90 degrees angled) axes. Your map seems to\n                    : have differently angled axes and so cannot be processed by\n                    : this ProSHADE version. This feature is coming in future\n                    : update!" );
    }
    
    //================================================ Release memory
    free                                              ( cell );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function or reading the map dimensions from the header.
 
 This function takes an open map file and reads the cell dimensions information, saving the read values to the supp;ied
 variables.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to read from.
 \param[in] xDim Address to a variable to save the x dimension size in indices.
 \param[in] yDim Address to a variable to save the y dimension size in indices.
 \param[in] zDim Address to a variable to save the z dimension size in indices.
 */
void ProSHADE_internal_io::readInMapDim ( CMap_io::CMMFile* mapFile, proshade_unsign* xDim, proshade_unsign* yDim, proshade_unsign* zDim )
{
    //================================================ Initialise variables
    int *dim                                          = NULL;
    
    //================================================ Allocate memory
    dim                                               = (int*) malloc (3 * sizeof ( int ) );
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( dim, __FILE__, __LINE__, __func__ );
    
    //================================================ Read the data
    CMap_io::ccp4_cmap_get_dim                        ( mapFile, dim );
    
    //================================================ Save the data
   *xDim                                              = static_cast<proshade_signed> ( dim[0] );
   *yDim                                              = static_cast<proshade_signed> ( dim[1] );
   *zDim                                              = static_cast<proshade_signed> ( dim[2] );
    
    //================================================ Release memory
    free                                              ( dim );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function or reading the map grid from the header.
 
 This function takes an open map file and reads the cell grid information, saving the read values to the supp;ied
 variables.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to read from.
 \param[in] xGrid Address to a variable to save the x dimension size in angstroms (as far as I believe).
 \param[in] yGrid Address to a variable to save the y dimension size in angstroms (as far as I believe).
 \param[in] zGrid Address to a variable to save the z dimension size in angstroms (as far as I believe).
 */
void ProSHADE_internal_io::readInMapGrid ( CMap_io::CMMFile* mapFile, proshade_unsign* xGrid, proshade_unsign* yGrid, proshade_unsign* zGrid )
{
    //================================================ Initialise variables
    int *grid                                         = NULL;
    
    //================================================ Allocate memory
    grid                                              = (int*) malloc (3 * sizeof ( int ) );
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( grid, __FILE__, __LINE__, __func__ );
    
    //================================================ Read the data
    CMap_io::ccp4_cmap_get_grid                       ( mapFile, grid );
    
    //================================================ Save the data
   *xGrid                                             = static_cast<proshade_signed> ( grid[0] );
   *yGrid                                             = static_cast<proshade_signed> ( grid[1] );
   *zGrid                                             = static_cast<proshade_signed> ( grid[2] );
    
    //================================================ Release memory
    free                                              ( grid );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function or reading the map dimension order from the header.
 
 This function takes an open map file and reads the order of axes, saving the read values to the supp;ied
 variables.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to read from.
 \param[in] xAxisOrder Address to a variable to save the x axis order.
 \param[in] yAxisOrder Address to a variable to save the y axis order.
 \param[in] zAxisOrder Address to a variable to save the z axis order.
 */
void ProSHADE_internal_io::readInMapOrder ( CMap_io::CMMFile* mapFile, proshade_unsign* xAxisOrder, proshade_unsign* yAxisOrder, proshade_unsign* zAxisOrder )
{
    //================================================ Initialise variables
    int *order                                        = NULL;
    
    //================================================ Allocate memory
    order                                             = (int*) malloc (3 * sizeof ( int ) );
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( order, __FILE__, __LINE__, __func__ );
    
    //================================================ Read the data
    CMap_io::ccp4_cmap_get_order                      ( mapFile, order );
    
    //================================================ Save the data
   *xAxisOrder                                        = static_cast<proshade_signed> ( order[0] );
   *yAxisOrder                                        = static_cast<proshade_signed> ( order[1] );
   *zAxisOrder                                        = static_cast<proshade_signed> ( order[2] );
    
    //================================================ Release memory
    free                                              ( order );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for reading the map origin from the header.
 
 This function takes an open map file and reads the origin of the map, saving the read values to the supp;ied
 variables.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to read from.
 \param[in] xOrigin Address to a variable to save the x axis origin position.
 \param[in] yOrigin Address to a variable to save the y axis origin position.
 \param[in] zOrigin Address to a variable to save the z axis origin position.
 */
void ProSHADE_internal_io::readInMapOrigin ( CMap_io::CMMFile* mapFile, proshade_signed* xOrigin, proshade_signed* yOrigin, proshade_signed* zOrigin )
{
    //================================================ Initialise variables
    int *origin                                       = NULL;
    
    //================================================ Allocate memory
    origin                                            = (int*) malloc (3 * sizeof ( int ) );
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( origin, __FILE__, __LINE__, __func__ );
    
    //================================================ Read the data
    CMap_io::ccp4_cmap_get_origin                     ( mapFile, origin );
    
    //================================================ Save the data
   *xOrigin                                           = static_cast<proshade_signed> ( origin[0] );
   *yOrigin                                           = static_cast<proshade_signed> ( origin[1] );
   *zOrigin                                           = static_cast<proshade_signed> ( origin[2] );
    
    //================================================ Release memory
    free                                              ( origin );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for reading the map data from a file.
 
 This function takes an open map file and reads the map data into the supplied map pointer.
 
 \param[in] mapFile Pointer to an open CMap_io::CMMFile object to read from.
 \param[in] map Pointer reference to a variable to save the map data.
 \param[in] xIndexStart Pointer to the x index start value.
 \param[in] yIndexStart Pointer to the y index start value.
 \param[in] zIndexStart Pointer to the z index start value.
 \param[in] xDimInIndices Pointer to the size of x dimension in indices.
 \param[in] yDimInIndices Pointer to the size of y dimension in indices.
 \param[in] zDimInIndices Pointer to the size of z dimension in indices.
 \param[in] xAxisOrder Pointer to the x axis order.
 \param[in] yAxisOrder Pointer to the y axis order.
 \param[in] zAxisOrder Pointer to the z axis order.
 */
void ProSHADE_internal_io::readInMapData ( CMap_io::CMMFile* mapFile, proshade_double*& map, proshade_signed xIndexStart, proshade_signed yIndexStart, proshade_signed zIndexStart, proshade_signed xDimInIndices, proshade_signed yDimInIndices, proshade_signed zDimInIndices, proshade_signed xAxisOrder, proshade_signed yAxisOrder, proshade_signed zAxisOrder )
{
    //================================================ Allocate map memory
    map                                               = new proshade_double [xDimInIndices * yDimInIndices * zDimInIndices];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( map, __FILE__, __LINE__, __func__ );
    
    //================================================ Check the map mode
    int mapMode                                       = CMap_io::ccp4_cmap_get_datamode ( mapFile );
    if ( ( mapMode != 0 ) && ( mapMode != 2 ) )
    {
        throw ProSHADE_exception ( "The map file mode is not supported.", "EM00009", __FILE__, __LINE__, __func__, "ProSHADE currently only supports map modes 0 and 2. The\n                    : supplied map seems to have different map mode and\n                    : therefore cannot be read in." );
    }
    
    //================================================ Read in map data
    // ... Prepare the variables
    proshade_signed index;
    proshade_signed iters[3];
    proshade_signed maxLim[3];
    proshade_signed XYZOrder[3];
    proshade_signed newU, newV, newW;
    proshade_signed arrPos;
    proshade_signed maxMapV = 0, maxMapW = 0;
    proshade_signed order[3]; order[0] = xAxisOrder;    order[1] = yAxisOrder;    order[2] = zAxisOrder;
    proshade_signed orig[3];  orig[0]  = xIndexStart;   orig[1]  = yIndexStart;   orig[2]  = zIndexStart;
    proshade_signed dim[3];   dim[0]   = xDimInIndices; dim[1]   = yDimInIndices; dim[2]   = zDimInIndices;
    
    // ... Set the dimensions for XYZ indexing
    if ( order[1] == 1 ) { maxMapV = dim[0] - 1; }
    if ( order[1] == 2 ) { maxMapV = dim[1] - 1; }
    if ( order[1] == 3 ) { maxMapV = dim[2] - 1; }
    
    if ( order[2] == 1 ) { maxMapW = dim[0] - 1; }
    if ( order[2] == 2 ) { maxMapW = dim[1] - 1; }
    if ( order[2] == 3 ) { maxMapW = dim[2] - 1; }
    
    // ... Set the XYZ indexing order and indices
    for ( proshade_signed iter = 0; iter < 3; iter++ )
    {
        maxLim[iter]                                  = orig[iter] + dim[iter] - 1;
        XYZOrder[order[iter]-1]                       = iter;
    }
    
    // ... Solve the dimensions and sizes for reading
    proshade_signed fastDimSize                       = ( maxLim[0] - orig[0] + 1 );
    proshade_signed midDimSize                        = ( maxLim[1] - orig[1] + 1 ) * fastDimSize;
    std::vector < float > section                     ( midDimSize );
    
    // ... Read in the map data
    for ( iters[2] = orig[2]; iters[2] <= maxLim[2]; iters[2]++ )
    {
        index                                         = 0;
        CMap_io::ccp4_cmap_read_section( mapFile, &section[0] );
        
        if ( mapMode == 0 )
        {
            for ( int iter = ( midDimSize - 1 ); iter >= 0; iter-- )
            {
                section[iter]                         = static_cast < float > ( ( reinterpret_cast<unsigned char*> (&section[0]) )[iter] );
            }
        }
        
        for ( iters[1] = orig[1]; iters[1] <= maxLim[1]; iters[1]++ )
        {
            for ( iters[0] = orig[0]; iters[0] <= maxLim[0]; iters[0]++ )
            {
                newU                                  = iters[XYZOrder[0]] - orig[XYZOrder[0]];
                newV                                  = iters[XYZOrder[1]] - orig[XYZOrder[1]];
                newW                                  = iters[XYZOrder[2]] - orig[XYZOrder[2]];
                arrPos                                = newW + ( maxMapW + 1) * ( newV + ( maxMapV + 1) * newU );
                map[arrPos]                           = static_cast < proshade_double > ( section[ index++ ] );
            }
        }
    }
    
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
