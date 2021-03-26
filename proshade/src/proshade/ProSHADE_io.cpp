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
    \version   0.7.5.4
    \date      MAR 2021
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
    std::memset                                       ( reinterpret_cast<void*> ( &(map->ccp4_header.at( 56 )) ), ' ', 800 + map->grid.spacegroup->operations().order() * 80); // 56 is used because the vector is indexed from 0
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
