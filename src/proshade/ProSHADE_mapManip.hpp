/*! \file ProSHADE_mapManip.hpp
    \brief This header file declares the ProSHADE_internal_mapManip namespace, which groups functions for internal map manipulation.
 
    The functions grouped here are used to modify the internal density map, or the MMDB2 library objects from which an internal map will be computed. These functions
    deal with the boundaries and their changes, re-sampling using trilinear interpolation or phase remova, to name a few examples.
 
    Copyright by Michal Tykac and individual contributors. All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
 
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.3
    \date      JUN 2020
 */

//==================================================== ProSHADE
#include "ProSHADE_io.hpp"

//==================================================== Overinclusion protection
#ifndef __PROSHADE_MAPMANIP__
#define __PROSHADE_MAPMANIP__

//==================================================== ProSHADE_internal_mapManip namespace
/*! \namespace ProSHADE_internal_mapManip
 \brief This namespace contains the internal functions for manipulating maps already present in the internal structures.
 
 The ProSHADE_internal_mapManip namespace contains helper functions for map manipulation and processing. However, these functions do make
 minimum assumptions about the map internal organisation and variables and simply perform tasks on generic variables. None of these
 functions should be used directly be the user.
 */
namespace ProSHADE_internal_mapManip
{
    proshade_signed myRound                           ( proshade_double x );
    void determinePDBRanges                           ( clipper::mmdb::CMMDBManager* pdbFile, proshade_single* xFrom, proshade_single* xTo, proshade_single* yFrom,
                                                        proshade_single* yTo, proshade_single* zFrom, proshade_single* zTo, int* noAt );
    void changePDBBFactors                            ( clipper::mmdb::CMMDBManager* pdbFile, proshade_double newBFactorValue );
    void movePDBForClipper                            ( clipper::mmdb::CMMDBManager* pdbFile, proshade_single xMov, proshade_single yMov, proshade_single zMov );
    void generateMapFromPDB                           ( clipper::mmdb::CMMDBManager* pdbFile, proshade_double*& map, proshade_single requestedResolution,
                                                        proshade_single xCell, proshade_single yCell, proshade_single zCell, int noAtoms, proshade_signed* xTo,
                                                        proshade_signed* yTo, proshade_signed* zTo );
    void moveMapByIndices                             ( proshade_single* xMov, proshade_single* yMov, proshade_single* zMov, proshade_single xAngs, proshade_single yAngs,
                                                        proshade_single zAngs, proshade_signed* xFrom, proshade_signed* xTo, proshade_signed* yFrom, proshade_signed* yTo,
                                                        proshade_signed* zFrom, proshade_signed* zTo, proshade_signed* xOrigin, proshade_signed* yOrigin,
                                                        proshade_signed* zOrigin );
    void moveMapByFourier                             ( proshade_double*& map, proshade_single xMov, proshade_single yMov, proshade_single zMov, proshade_single xAngs,
                                                        proshade_single yAngs, proshade_single zAngs, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim );
    void blurSharpenMap                               ( proshade_double*& map, proshade_double*& maskedMap, proshade_unsign xDimS, proshade_unsign yDimS,
                                                        proshade_unsign zDimS, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_single blurringFactor );
    void getMaskFromBlurr                             ( proshade_double*& blurMap, proshade_double*& outMap, proshade_unsign xDimS, proshade_unsign yDimS,
                                                        proshade_unsign zDimS, proshade_single noIQRs );
    void getNonZeroBounds                             ( proshade_double* map, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim, proshade_single xAngs,
                                                        proshade_single yAngs, proshade_single zAngs, proshade_signed*& ret );
    void addExtraBoundSpace                           ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_single xAngs, proshade_single yAngs,
                                                        proshade_single zAngs, proshade_signed*& bounds, proshade_single extraSpace );
    void reSampleMapToResolutionTrilinear             ( proshade_double*& map, proshade_single resolution, proshade_unsign xDimS, proshade_unsign yDimS,
                                                        proshade_unsign zDimS, proshade_single xAngs, proshade_single yAngs,  proshade_single zAngs,
                                                        proshade_single*& corrs );
    void removeMapPhase                               ( fftw_complex*& mapCoeffs, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim );
    void getFakeHalfMap                               ( proshade_double*& map, proshade_double*& fakeHalfMap, proshade_unsign xDimS, proshade_unsign yDimS,
                                                       proshade_unsign zDimS, proshade_signed fakeMapKernel );
    void getCorrelationMapMask                        ( proshade_double*& map, proshade_double*& fakeHalfMap, proshade_double*& correlationMask,
                                                        proshade_unsign xDimS, proshade_unsign yDimS,
                                                       proshade_unsign zDimS, proshade_signed corrMaskKernel );
    proshade_single getIndicesFromAngstroms           ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim,
                                                        proshade_single xAngs, proshade_single yAngs, proshade_single zAngs,
                                                        proshade_single dist );
    void connectMaskBlobs                             ( proshade_double*& mask, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim,
                                                        proshade_single xAngs, proshade_single yAngs,
                                                        proshade_single zAngs, proshade_single maskThres );
    void beautifyBoundaries                           ( proshade_signed*& bounds, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim,
                                                        proshade_signed boundsDiffThres, proshade_signed verbose );
    proshade_signed betterClosePrimeFactors           ( proshade_signed fromRange, proshade_signed toRange );
    void distributeSpaceToBoundaries                  ( proshade_signed& minBound, proshade_signed& maxBound, proshade_signed oldBoundRange, proshade_signed newBoundRange );
    void copyMapByBounds                              ( proshade_signed xFrom, proshade_signed xTo, proshade_signed yFrom, proshade_signed yTo, proshade_signed zFrom,
                                                        proshade_signed zTo, proshade_signed origXFrom, proshade_signed origYFrom, proshade_signed origZFrom,
                                                        proshade_signed yDimIndices, proshade_signed zDimIndices, proshade_signed origXDimIndices,
                                                        proshade_signed origYDimIndices, proshade_signed origZDimIndices, proshade_double*& newMap,
                                                        proshade_double* origMap );
}

#endif
