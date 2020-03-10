/*! \file ProSHADE_mapManip.hpp
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
#include "ProSHADE_io.hpp"

#ifndef __PROSHADE_MAPMANIP__
#define __PROSHADE_MAPMANIP__

//============================================ ProSHADE_internal_mapManip namespace
/*! \namespace ProSHADE_internal_mapManip
 \brief This namespace contains the internal functions for manipulating maps already present in the internal structures.
 
 The ProSHADE_internal_mapManip namespace contains helper functions for map manipulation and processing. However, these functions do make
 minimum assumptions about the map internal organisation and variables and simply perform tasks on generic variables. None of these
 functions should be used directly be the user.
 */
namespace ProSHADE_internal_mapManip
{
    proshade_signed myRound                   ( proshade_double x );
    void determinePDBRanges                   ( clipper::mmdb::CMMDBManager* pdbFile, proshade_single* xFrom, proshade_single* xTo, proshade_single* yFrom,
                                                proshade_single* yTo, proshade_single* zFrom, proshade_single* zTo, int* noAt );
    void changePDBBFactors                    ( clipper::mmdb::CMMDBManager* pdbFile, proshade_double newBFactorValue );
    void movePDBForClipper                    ( clipper::mmdb::CMMDBManager* pdbFile, proshade_single xMov, proshade_single yMov, proshade_single zMov );
    void generateMapFromPDB                   ( clipper::mmdb::CMMDBManager* pdbFile, proshade_double*& map, proshade_single requestedResolution,
                                                proshade_single xCell, proshade_single yCell, proshade_single zCell, int noAtoms, proshade_signed* xTo,
                                                proshade_signed* yTo, proshade_signed* zTo );
    void moveMapByIndices                     ( proshade_single* xMov, proshade_single* yMov, proshade_single* zMov, proshade_single xAngs, proshade_single yAngs,
                                                proshade_single zAngs, proshade_signed* xFrom, proshade_signed* xTo, proshade_signed* yFrom, proshade_signed* yTo,
                                                proshade_signed* zFrom, proshade_signed* zTo, proshade_signed* xOrigin, proshade_signed* yOrigin,
                                                proshade_signed* zOrigin );
    void moveMapByFourier                     ( proshade_double*& map, proshade_single xMov, proshade_single yMov, proshade_single zMov, proshade_single xAngs,
                                                proshade_single yAngs, proshade_single zAngs, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim );
    void blurSharpenMap                       ( proshade_double*& map, proshade_double*& maskedMap, proshade_unsign xDimS, proshade_unsign yDimS,
                                                proshade_unsign zDimS, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_single blurringFactor );
    void getMaskFromBlurr                     ( proshade_double*& blurMap, proshade_double*& outMap, proshade_unsign xDimS, proshade_unsign yDimS,
                                                proshade_unsign zDimS, proshade_single noIQRs );
    void getNonZeroBounds                     ( proshade_double* map, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim, proshade_single xAngs,
                                                proshade_single yAngs, proshade_single zAngs, proshade_signed*& ret );
    void addExtraBoundSpace                   ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_single xAngs, proshade_single yAngs,
                                                proshade_single zAngs, proshade_signed*& bounds, proshade_single extraSpace );
    void reSampleMapToResolutionTrilinear     ( proshade_double*& map, proshade_single resolution, proshade_unsign xDimS, proshade_unsign yDimS,
                                                proshade_unsign zDimS, proshade_single xAngs, proshade_single yAngs,  proshade_single zAngs,
                                                proshade_single*& corrs );
    void removeMapPhase                       ( fftw_complex*& mapCoeffs, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim );
    void getFakeHalfMap                       ( proshade_double*& map, proshade_double*& fakeHalfMap, proshade_unsign xDimS, proshade_unsign yDimS,
                                               proshade_unsign zDimS, proshade_signed fakeMapKernel );
    void getCorrelationMapMask                ( proshade_double*& map, proshade_double*& fakeHalfMap, proshade_double*& correlationMask, proshade_unsign xDimS, proshade_unsign yDimS,
                                               proshade_unsign zDimS, proshade_signed corrMaskKernel );
    proshade_single getIndicesFromAngstroms   ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs,
                                                proshade_single dist );
    void connectMaskBlobs                     ( proshade_double*& mask, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim, proshade_single xAngs, proshade_single yAngs,
                                                proshade_single zAngs, proshade_single maskThres );
    void beautifyBoundaries                   ( proshade_signed*& bounds, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_signed boundsDiffThres, proshade_signed verbose );
    proshade_signed betterClosePrimeFactors   ( proshade_signed fromRange, proshade_signed toRange );
    void distributeSpaceToBoundaries          ( proshade_signed& minBound, proshade_signed& maxBound, proshade_signed oldBoundRange, proshade_signed newBoundRange );
    void copyMapByBounds                      ( proshade_signed xFrom, proshade_signed xTo, proshade_signed yFrom, proshade_signed yTo, proshade_signed zFrom, proshade_signed zTo,
                                                proshade_signed origXFrom, proshade_signed origYFrom, proshade_signed origZFrom, proshade_signed yDimIndices, proshade_signed zDimIndices,
                                                proshade_signed origXDimIndices, proshade_signed origYDimIndices, proshade_signed origZDimIndices, proshade_double*& newMap,
                                                proshade_double* origMap );
}

#endif
