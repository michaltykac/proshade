/*! \file ProSHADE_io.hpp
    \brief This header file declares all the functions required for low level file format access.
 
    This header file declares the ProSHADE_internal_io namespace, which groups together all the low level functions required for input and output file format detection, reading and writing. Also,
    most of the interactions with the CMAPLIB CCP4 library are done by these functions as to avoid these being dispersed everywhere.
 
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
#include "ProSHADE_settings.hpp"

//==================================================== Overinclusion protection
#ifndef __PROSHADE_IO__
#define __PROSHADE_IO__

//==================================================== ProSHADE_internal_io Namespace
/*! \namespace ProSHADE_internal_io
 \brief This namespace contains the internal input/output functions. None of these should be directly accessed by the user.
 
 The ProSHADE_internal_io namespace contains the helper functions for the data input and output. These should never be
 directly used by the user and these only serve to allow for self-documenting nature of the code. They are called internally
 by more advanced functions from the higher complexity classes.
 */
namespace ProSHADE_internal_io
{
    //================================================ The InputType data type
    enum InputType                                    { UNKNOWN, PDB, MAP };
    
    //================================================ Low level file access functions
    InputType figureDataType                          ( std::string fName );
    bool isFilePDB                                    ( std::string fName );
    bool isFileMAP                                    ( std::string fName );
    void readInMapCell                                ( CMap_io::CMMFile* mapFile, proshade_single* xDim, proshade_single* yDim, proshade_single* zDim,
                                                        proshade_single* aAng, proshade_single* bAng, proshade_single* cAng );
    void readInMapDim                                 ( CMap_io::CMMFile* mapFile, proshade_unsign* xDim, proshade_unsign* yDim, proshade_unsign* zDim );
    void readInMapGrid                                ( CMap_io::CMMFile* mapFile, proshade_unsign* xGrid, proshade_unsign* yGrid, proshade_unsign* zGrid );
    void readInMapOrder                               ( CMap_io::CMMFile* mapFile, proshade_unsign* xAxisOrder, proshade_unsign* yAxisOrder, proshade_unsign* zAxisOrder );
    void readInMapOrigin                              ( CMap_io::CMMFile* mapFile, proshade_signed* xOrigin, proshade_signed* yOrigin, proshade_signed* zOrigin );
    void readInMapData                                ( CMap_io::CMMFile* mapFile, proshade_double*& map, proshade_signed xIndexStart, proshade_signed yIndexStart,
                                                        proshade_signed zIndexStart, proshade_signed xDimInIndices, proshade_signed yDimInIndices,
                                                        proshade_signed zDimInIndices, proshade_signed xAxisOrder, proshade_signed yAxisOrder, proshade_signed zAxisOrder );
    void writeMapCell                                 ( CMap_io::CMMFile* mapFile, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs,
                                                        proshade_single aAngle, proshade_single bAngle, proshade_single cAngle );
    void writeMapGrid                                 ( CMap_io::CMMFile* mapFile, proshade_unsign xGrid, proshade_unsign yGrid, proshade_unsign zGrid );
    void writeMapOrder                                ( CMap_io::CMMFile* mapFile, proshade_unsign xAxisOrder, proshade_unsign yAxisOrder, proshade_unsign zAxisOrder );
    void writeMapDims                                 ( CMap_io::CMMFile* mapFile, proshade_unsign xDims, proshade_unsign yDims, proshade_unsign zDims );
    void writeMapOrigin                               ( CMap_io::CMMFile* mapFile, proshade_unsign xOrigin, proshade_unsign yOrigin, proshade_unsign zOrigin );
    void writeMapTitleEtc                             ( CMap_io::CMMFile* mapFile, std::string title, proshade_unsign mode = 2, proshade_unsign spaceGroup = 1 );
    void writeMapData                                 ( CMap_io::CMMFile* mapFile, proshade_double* map, proshade_unsign xDim, proshade_unsign yDim,
                                                        proshade_unsign zDim, proshade_unsign xAxisOrder, proshade_unsign yAxisOrder, proshade_unsign zAxisOrder );
}

#endif
