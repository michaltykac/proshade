/*! \file ProSHADE_io.hpp
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
#include "ProSHADE_settings.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE_IO__
#define __PROSHADE_IO__

//============================================ ProSHADE_internal_io Namespace
/*! \namespace ProSHADE_internal_io
 \brief This namespace contains the internal input/output functions. None of these should be directly accessed by the user.
 
 The ProSHADE_internal_io namespace contains the helper functions for the data input and output. These should never be
 directly used by the user and these only serve to allow for self-documenting nature of the code. They are called internally
 by more advanced functions from the higher complexity classes.
 */
namespace ProSHADE_internal_io
{
    //======================================== The InputType data type
    enum InputType                            { UNKNOWN, PDB, MAP };
    
    InputType figureDataType                  ( std::string fName );
    bool isFilePDB                            ( std::string fName );
    bool isFileMAP                            ( std::string fName );
    void readInMapCell                        ( CMap_io::CMMFile* mapFile, proshade_single* xDim, proshade_single* yDim, proshade_single* zDim,
                                                proshade_single* aAng, proshade_single* bAng, proshade_single* cAng );
    void readInMapDim                         ( CMap_io::CMMFile* mapFile, proshade_unsign* xDim, proshade_unsign* yDim, proshade_unsign* zDim );
    void readInMapGrid                        ( CMap_io::CMMFile* mapFile, proshade_unsign* xGrid, proshade_unsign* yGrid, proshade_unsign* zGrid );
    void readInMapOrder                       ( CMap_io::CMMFile* mapFile, proshade_unsign* xAxisOrder, proshade_unsign* yAxisOrder, proshade_unsign* zAxisOrder );
    void readInMapOrigin                      ( CMap_io::CMMFile* mapFile, proshade_signed* xOrigin, proshade_signed* yOrigin, proshade_signed* zOrigin );
    void readInMapData                        ( CMap_io::CMMFile* mapFile, proshade_double*& map, proshade_signed xIndexStart, proshade_signed yIndexStart,
                                                proshade_signed zIndexStart, proshade_signed xDimInIndices, proshade_signed yDimInIndices,
                                                proshade_signed zDimInIndices, proshade_signed xAxisOrder, proshade_signed yAxisOrder, proshade_signed zAxisOrder );
    void writeMapCell                         ( CMap_io::CMMFile* mapFile, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs,
                                                proshade_single aAngle, proshade_single bAngle, proshade_single cAngle );
    void writeMapGrid                         ( CMap_io::CMMFile* mapFile, proshade_unsign xGrid, proshade_unsign yGrid, proshade_unsign zGrid );
    void writeMapOrder                        ( CMap_io::CMMFile* mapFile, proshade_unsign xAxisOrder, proshade_unsign yAxisOrder, proshade_unsign zAxisOrder );
    void writeMapDims                         ( CMap_io::CMMFile* mapFile, proshade_unsign xDims, proshade_unsign yDims, proshade_unsign zDims );
    void writeMapOrigin                       ( CMap_io::CMMFile* mapFile, proshade_unsign xOrigin, proshade_unsign yOrigin, proshade_unsign zOrigin );
    void writeMapTitleEtc                     ( CMap_io::CMMFile* mapFile, std::string title, proshade_unsign mode = 2, proshade_unsign spaceGroup = 1 );
    void writeMapData                         ( CMap_io::CMMFile* mapFile, proshade_double* map, proshade_unsign xDim, proshade_unsign yDim,
                                                proshade_unsign zDim, proshade_unsign xAxisOrder, proshade_unsign yAxisOrder, proshade_unsign zAxisOrder );
}

#endif
