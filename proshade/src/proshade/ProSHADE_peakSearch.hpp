/*! \file ProSHADE_peakSearch.hpp
    \brief This header file declares functions required for 3D map peak searching.
 
    This header file declares the ProSHADE_internal_peakSearch namespace, which groups all the ProSHADE functions used to search for 3D density
    map peaks and processing them, including peak position optimisation.
 
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
#include "ProSHADE_spheres.hpp"

//==================================================== Overinclusion protection
#ifndef PROSHADE_PEAKSEARCH
#define PROSHADE_PEAKSEARCH

//==================================================== ProSHADE_internal_peakSearch Namespace
/*! \namespace ProSHADE_internal_peakSearch
    \brief This namespace contains all the functions required for peak searching in 3D maps.
 
    Peak searching may be one of the bottlenecks of the "old" ProSHADE and therefore a dedicated header with
    functions related to 3D peak searching is created to allow more testing and work to be done on this task.
 */
namespace ProSHADE_internal_peakSearch
{
    std::vector< proshade_double* > findAllPointsAboveNeighbours ( proshade_complex* map, proshade_unsign dim,
                                                                   proshade_signed peakSize, proshade_double* medianIQR );
    void pointsAboveNeighboursRemoveSmallHeight       ( std::vector< proshade_double* >* pointVec, proshade_double* medianIQR,
                                                        proshade_double noIQRs );
    void allocatePeakOptimisationMemory               ( proshade_double*& avgMat, proshade_double*& hlpMap, proshade_double*& eA,
                                                        proshade_double*& eB, proshade_double*& eG, proshade_double*& uAV );
    void releasePeakOptimisationMemory                ( proshade_double*& avgMat, proshade_double*& hlpMap, proshade_double*& eA,
                                                       proshade_double*& eB, proshade_double*& eG, proshade_double*& uAV );
    void optimisePeakPositions                        ( std::vector< proshade_double* >* pointVec, proshade_signed peakSize, proshade_signed band );
    std::vector< proshade_double* > getAllPeaksNaive  ( proshade_complex* map, proshade_unsign dim, proshade_signed peakSize, proshade_double noIQRs );
    void getBestPeakEulerAngsNaive                    ( proshade_complex* map, proshade_unsign dim, proshade_double* eulA, proshade_double* eulB,
                                                        proshade_double* eulG, ProSHADE_settings* settings );
    void allocateSmoothingZScoreMemory                ( proshade_unsign dim, proshade_double*& scoreOverVals, proshade_signed*& signals,
                                                        proshade_double*& filteredY, proshade_double*& avgFilter, proshade_double*& stdFilter,
                                                        proshade_double*& subVec, proshade_double*& medianIQR, proshade_double*& YZMap,
                                                        proshade_double*& XZMap, proshade_double*& XYMap, proshade_unsign smLag );
    void releaseSmoothingZScoreMemory                 ( proshade_double*& scoreOverVals, proshade_signed*& signals, proshade_double*& filteredY,
                                                        proshade_double*& avgFilter, proshade_double*& stdFilter, proshade_double*& subVec,
                                                        proshade_double*& medianIQR, proshade_double*& YZMap, proshade_double*& XZMap,
                                                        proshade_double*& XYMap );
    void getSmoothedZScore1D                          ( proshade_unsign dim, proshade_unsign smoothingLag, proshade_double ZScoreThreshold,
                                                        proshade_signed* signals, proshade_double* filteredY, proshade_double* avgFilter,
                                                        proshade_double* stdFilter, proshade_double* subVec, proshade_double* medianIQR,
                                                        proshade_double* scoreOverVals );
    void getXAxisArraysSmoothedZScorePeaks            ( proshade_unsign dim, proshade_unsign smoothingLag, proshade_double ZScoreThreshold,
                                                        proshade_signed* signals, proshade_double* filteredY, proshade_double* avgFilter,
                                                        proshade_double* stdFilter, proshade_double* subVec, proshade_double* medianIQR,
                                                        proshade_double* scoreOverVals, proshade_complex* map, proshade_double* YZMap );
    void getYAxisArraysSmoothedZScorePeaks            ( proshade_unsign dim, proshade_unsign smoothingLag, proshade_double ZScoreThreshold,
                                                        proshade_signed* signals, proshade_double* filteredY, proshade_double* avgFilter,
                                                        proshade_double* stdFilter, proshade_double* subVec, proshade_double* medianIQR,
                                                        proshade_double* scoreOverVals, proshade_complex* map, proshade_double* XZMap );
    void getZAxisArraysSmoothedZScorePeaks            ( proshade_unsign dim, proshade_unsign smoothingLag, proshade_double ZScoreThreshold,
                                                        proshade_signed* signals, proshade_double* filteredY, proshade_double* avgFilter,
                                                        proshade_double* stdFilter, proshade_double* subVec, proshade_double* medianIQR,
                                                        proshade_double* scoreOverVals, proshade_complex* map, proshade_double* XYMap );
    void findAllPointNeighbours                       ( proshade_double* YZMap, proshade_double* XZMap, proshade_double* XYMap, proshade_unsign* visitedMap,
                                                        proshade_signed dim, proshade_signed x, proshade_signed y, proshade_signed z,
                                                        std::vector< proshade_unsign >* retVals );
    void findAllDisconnectedIslands                   ( proshade_complex* map, proshade_double* YZMap, proshade_double* XZMap, proshade_double* XYMap,
                                                        proshade_unsign dim, std::vector< proshade_unsign >* allIslandBests );
    void findAllSmoothedZScorePeaksWithNeighbours     ( proshade_complex* map, proshade_double* YZMap, proshade_double* XZMap, proshade_double* XYMap,
                                                        proshade_signed dim, proshade_signed peakSize,
                                                        std::vector< proshade_double* >* allPeaksWithNeighbours );
    std::vector< proshade_double* > getAllPeaksSmoothedZ ( proshade_complex* map, proshade_unsign dim, proshade_double smoothingFraction,
                                                           proshade_double noIQRs, proshade_signed peakSize );
    void getBestPeakEulerAngsSmoothedZ                ( proshade_complex* map, proshade_unsign dim, proshade_double smoothingFraction,
                                                        proshade_double noIQRs, proshade_signed peakSize, proshade_double* eulA, proshade_double* eulB,
                                                        proshade_double* eulG );
}

#endif
