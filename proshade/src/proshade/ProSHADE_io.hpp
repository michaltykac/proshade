/*! \file ProSHADE_io.hpp
    \brief This header file declares all the functions required for low level file format access.
 
    This header file declares the ProSHADE_internal_io namespace, which groups together all the low level functions required for input and output file format detection, reading and writing. Also,
    most of the interactions with the Gemmi library are done by these functions as to avoid these being dispersed everywhere.
 
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
#include "ProSHADE_settings.hpp"

//==================================================== Overinclusion protection
#ifndef PROSHADE_IO
#define PROSHADE_IO

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
    enum InputType                                    { UNKNOWN, PDB, MAP, GEMMI };
    
    //================================================ Low level file access functions
    InputType figureDataType                          ( std::string fName );
    bool isFilePDB                                    ( std::string fName );
    bool isFileMAP                                    ( std::string fName );
    void readInMapHeader                              ( gemmi::Ccp4<float> *map, proshade_unsign *xDimInds, proshade_unsign *yDimInds, proshade_unsign *zDimInds, proshade_single *xDim,
                                                        proshade_single *yDim, proshade_single *zDim, proshade_single *aAng, proshade_single *bAng, proshade_single *cAng, proshade_signed *xFrom,
                                                        proshade_signed *yFrom, proshade_signed *zFrom, proshade_signed *xAxOrigin, proshade_signed *yAxOrigin, proshade_signed *zAxOrigin,
                                                        proshade_unsign *xAxOrder, proshade_unsign *yAxOrder, proshade_unsign *zAxOrder, proshade_unsign *xGridInds, proshade_unsign *yGridInds,
                                                        proshade_unsign *zGridInds );
    void readInMapData                                ( gemmi::Ccp4<float> *gemmiMap, proshade_double*& map, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds,
                                                        proshade_unsign xAxOrder, proshade_unsign yAxOrder, proshade_unsign zAxOrder );
    void readInMapData                                ( gemmi::Ccp4<int8_t> *gemmiMap, proshade_double*& map, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds,
                                                        proshade_unsign xAxOrder, proshade_unsign yAxOrder, proshade_unsign zAxOrder );
    void applyMask                                    ( proshade_double*& map, std::string maskFile, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds,
                                                        proshade_signed verbose, proshade_signed messageShift, proshade_double* maskArray = nullptr, proshade_unsign maXInds = 0,
                                                        proshade_unsign maYInds = 0, proshade_unsign maZInds = 0 );
    void applyMaskFromArray                           ( proshade_double*& map, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds, proshade_double*& mask,
                                                        proshade_unsign xDimIndsMsk, proshade_unsign yDimIndsMsk, proshade_unsign zDimIndsMsk, proshade_signed verbose, proshade_signed messageShift );
    void applyWeights                                 ( proshade_double*& map, std::string weightsFile, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds,
                                                        proshade_signed verbose, proshade_signed messageShift, proshade_double* weightsArray = nullptr, proshade_unsign waXInds = 0,
                                                        proshade_unsign waYInds = 0, proshade_unsign waZInds = 0 );
    void applyWeightsFromArray                        ( proshade_double*& map, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds, proshade_double*& weights,
                                                        proshade_unsign xDimIndsWgh, proshade_unsign yDimIndsWgh, proshade_unsign zDimIndsWgh, proshade_signed verbose, proshade_signed messageShift );
    void writeOutMapHeader                            ( gemmi::Ccp4<float> *map, proshade_unsign xDimInds, proshade_unsign yDimInds, proshade_unsign zDimInds, proshade_single xDim,
                                                        proshade_single yDim, proshade_single zDim, proshade_single aAng, proshade_single bAng, proshade_single cAng, proshade_signed xFrom,
                                                        proshade_signed yFrom, proshade_signed zFrom, proshade_signed xAxOrigin, proshade_signed yAxOrigin, proshade_signed zAxOrigin,
                                                        proshade_unsign xAxOrder, proshade_unsign yAxOrder, proshade_unsign zAxOrder, proshade_unsign xGridInds, proshade_unsign yGridInds,
                                                       proshade_unsign zGridInds, std::string title, int mode );
    void writeRotationTranslationJSON                 ( proshade_double trsX1, proshade_double trsY1, proshade_double trsZ1, proshade_double eulA, proshade_double eulB, proshade_double eulG,
                                                        proshade_double trsX2, proshade_double trsY2, proshade_double trsZ2, std::string fileName );
}

#endif
