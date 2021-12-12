/*! \file ProSHADE_peakSearch.cpp
    \brief This source file declares functions required for peak searching and peak position optimisation.
 
    The functions declared in this source file all deal with 3D peak detection in density maps or with the optimisation of the exact peak position in such a maps.
    There are actually two sets of functions for this purpose, the naive approach and the Smoothed Z Score approach, however, the SZA is not full working at the
    moment (and as it is right now, it is bested by the neive approach is quite a few cases).
 
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
#include "ProSHADE_peakSearch.hpp"

/*! \brief This function finds all indices with higher value then all neighbours.
 
    This function finds all map indices with higher height than all their neighbours (the size of neighbours in terms of surrounding points can be changed)
    and saves these as well as the neighbours into the output vector of double*. It also keeps track of the non-peak values and saves their median height
    and its IQR, so that the output vector could be cleared from these small (background) values.
 
    \param[in] map Pointer to 1D array holding the 3D map value in which the peaks are to be found. Map must be cube!
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] peakSize The number of neighbouring points in single direction which should be considered as neighbours.
    \param[in] medianIQR Pointer to double[2] array where median and IQR of non-peak values will be saved.
    \param[out] X Vector of pointers to double numbers, first is the 'peak' (x,y,z and height) and then all neighbours follow in sets of 4 numbers.
 */
std::vector< proshade_double* > ProSHADE_internal_peakSearch::findAllPointsAboveNeighbours ( proshade_complex* map, proshade_unsign dim, proshade_signed peakSize, proshade_double* medianIQR )
{
    //================================================ Initialise local variables
    std::vector< proshade_double* > ret;
    std::vector< proshade_double > nonPeakVals;
    proshade_double* retHlp                           = nullptr;
    proshade_double pointHeight                       = 0.0;
    proshade_signed x, y, z, peakX, peakY, peakZ, newIter, ptrIter;
    proshade_signed xDim                              = static_cast< proshade_signed > ( std::pow ( dim, 2 ) );
    proshade_signed yDim                              = static_cast< proshade_signed > ( dim );
    bool breakPeak;
    
    //================================================ Check all points
    for ( proshade_signed iter = 0; iter < static_cast< proshade_signed > ( pow( static_cast<proshade_unsign> ( dim ), 3 ) ); iter++ )
    {
        //============================================ Find point height
        pointHeight                                   = pow( map[iter][0], 2.0 ) + pow( map[iter][1], 2.0 );
        
        //============================================ Find the x, y and z
        x                                             = static_cast< proshade_signed > ( std::floor ( iter / xDim ) );
        y                                             = static_cast< proshade_signed > ( std::floor ( ( iter - x * xDim ) / yDim ) );
        z                                             = static_cast< proshade_signed > ( iter - x * xDim - y * yDim );
        
        //==================================== Initialise the output pointer
        if ( retHlp == nullptr )
        {
            retHlp                                    = new proshade_double[static_cast<proshade_unsign>( pow( ((peakSize*2)+1), 3) * 4 )];
            ProSHADE_internal_misc::checkMemoryAllocation ( retHlp, __FILE__, __LINE__, __func__ );
        }
        
        //============================================ Get the X surrounding values and check for all being lower than point in question
        breakPeak                                     = false;
        ptrIter                                       = 4;
        for ( proshade_signed xCh = -peakSize; xCh <= +peakSize; xCh++ )
        {
            if ( breakPeak ) { break; }
            for ( proshade_signed yCh = -peakSize; yCh <= +peakSize; yCh++ )
            {
                if ( breakPeak ) { break; }
                for ( proshade_signed zCh = -peakSize; zCh <= +peakSize; zCh++ )
                {
                    if ( breakPeak ) { break; }
                    if ( ( xCh == 0 ) && ( yCh == 0 ) && ( zCh == 0 ) ) { continue; }
                    
                    //================================ Find the nieghbout peak indices (with periodicity)
                    peakX                             = x + xCh; if ( peakX >= yDim ) { peakX = yDim - 1; }; if ( peakX < 0 ) { peakX = 0; }
                    peakY                             = y + yCh; if ( peakY >= yDim ) { peakY = yDim - 1; }; if ( peakY < 0 ) { peakY = 0; }
                    peakZ                             = z + zCh; if ( peakZ >= yDim ) { peakZ = yDim - 1; }; if ( peakZ < 0 ) { peakZ = 0; }
                    newIter                           = peakX * xDim + peakY * yDim + peakZ;
                    
                    //================================ If neighbour has higher value than index in question, break out
                    if ( ( pow ( map[newIter][0], 2.0 ) + pow( map[newIter][1], 2.0 ) ) > pointHeight ) { breakPeak = true; break; }
                    
                    //================================ Save neighbour values for optimisation of peaks later
                    retHlp[ptrIter]                   = static_cast<proshade_double> ( peakX );
                    retHlp[ptrIter+1]                 = static_cast<proshade_double> ( peakY );
                    retHlp[ptrIter+2]                 = static_cast<proshade_double> ( peakZ );
                    retHlp[ptrIter+3]                 = pow ( map[newIter][0], 2.0 ) + pow( map[newIter][1], 2.0 );
                    ptrIter                          += 4;
                }
            }
        }
        if ( breakPeak ) { ProSHADE_internal_misc::addToDoubleVector ( &nonPeakVals, pointHeight ); continue; }
        
        //============================================ If passing, save for returning
        retHlp[0]                                     = static_cast< proshade_double > ( x );
        retHlp[1]                                     = static_cast< proshade_double > ( y );
        retHlp[2]                                     = static_cast< proshade_double > ( z );
        retHlp[3]                                     = pointHeight;
        ProSHADE_internal_misc::addToDblPtrVector ( &ret, retHlp );
        retHlp                                        = nullptr;
    }
    
    //================================================ Save non-peak median and IQR
    ProSHADE_internal_maths::vectorMedianAndIQR       ( &nonPeakVals, medianIQR );
    
    //================================================ Release memory if need be
    if ( retHlp != nullptr ) { delete[] retHlp; }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function clears the 'higher than neighbour' vector from background values.
 
    This function simply takes all the 'higher than all neighbours' indices vector and a threshold defining values and proceeds to
    compute the threshold for calling the indice a peak. Then, it applies this threshold to the vector, keeping only the indices
    which have higher value than the threshold.
 
    \param[in] pointVec Vector of double pointers as returned by findAllPointsAboveNeighbours().
    \param[in] medianIQR Array of two numbers, the median and IQR for which should be used to remove 'background' points in the vector.
    \param[in] noIQRs The number of IQRs from median to be used to create the cut-off threshold.
 */
void ProSHADE_internal_peakSearch::pointsAboveNeighboursRemoveSmallHeight ( std::vector< proshade_double* >* pointVec, proshade_double* medianIQR, proshade_double noIQRs )
{
    //================================================ Determine the threshold
    proshade_double backgroundThreshold               = std::min ( std::max ( medianIQR[0] + ( medianIQR[1] * noIQRs ), 0.05 ), 0.3 );
    
    //================================================ Check for passing the threshold
    for ( proshade_signed iter = static_cast<proshade_signed> ( pointVec->size()-1 ); iter >= 0  ; iter-- )
    {
        if ( pointVec->at( static_cast< size_t > ( iter ) )[3] <= backgroundThreshold )
        {
            delete[] pointVec->at( static_cast< size_t > ( iter ) );
            pointVec->erase                           ( pointVec->begin() + iter );
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allocates and checks all the peak optimisation memory.
 
    \param[in] avgMat The matrix which will take the weighted sum of the neighbour values.
    \param[in] hlpMat The helper matrix to store temporary results.
    \param[in] eA The Euler angle alpha value.
    \param[in] eB The Euler beta alpha value.
    \param[in] eG The Euler gamma alpha value.
    \param[in] uAV The U and V^T matrices resulting from the SVD.
 */
void ProSHADE_internal_peakSearch::allocatePeakOptimisationMemory ( proshade_double*& avgMat, proshade_double*& hlpMat, proshade_double*& eA, proshade_double*& eB, proshade_double*& eG, proshade_double*& uAV )
{
    //================================================ Allocate the memory
    avgMat                                            = new proshade_double [9];
    hlpMat                                            = new proshade_double [9];
    eA                                                = new proshade_double;
    eB                                                = new proshade_double;
    eG                                                = new proshade_double;
    uAV                                               = new proshade_double [18];
    
    //================================================ Check the memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( avgMat, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( hlpMat, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( eA,     __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( eB,     __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( eG,     __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( uAV,    __FILE__, __LINE__, __func__ );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function deletes all the peak optimisation memory.
 
    \param[in] avgMat The matrix which will take the weighted sum of the neighbour values.
    \param[in] hlpMat The helper matrix to store temporary results.
    \param[in] eA The Euler angle alpha value.
    \param[in] eB The Euler beta alpha value.
    \param[in] eG The Euler gamma alpha value.
    \param[in] uAV The U and V^T matrices resulting from the SVD.
 */
void ProSHADE_internal_peakSearch::releasePeakOptimisationMemory ( proshade_double*& avgMat, proshade_double*& hlpMat, proshade_double*& eA, proshade_double*& eB, proshade_double*& eG, proshade_double*& uAV )
{
    //================================================ Release the memory
    delete[] avgMat;
    delete[] hlpMat;
    delete[] eA;
    delete[] eB;
    delete[] eG;
    delete[] uAV;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function optimises all the peaks in the input vector using the values of their neighbours.
 
    This function takes the vector of arrays containing the peak and all its neighbours and proceeds to obtain the
    rotation matrices for the peak and all its neighbours. It then computes the weighted average matrix from the combination
    of all the matrices. By subjecting this weighted average matrix to SVD (using LAPACK) and combining the two rotation
    matrices (U and V^T), the optimised rotation matrix is obtained. This rotation matrix is then converted to Euler angles
    representation and these values are finally saved over the old all neighbours containing pointers. Thus, the number
    of peaks does not change, but each peak now only has 4 values - alpha, beta, gamma and peak height.
 
    \param[in] pointVec Vector of double pointers as returned by pointsAboveNeighboursRemoveSmallHeight().
    \param[in] peakSize The number of neighbouring points in single direction which should be considered as neighbours.
    \param[in] band The maximum bandwidth of the computation.
 */
void ProSHADE_internal_peakSearch::optimisePeakPositions ( std::vector< proshade_double* >* pointVec, proshade_signed peakSize, proshade_signed band )
{
    //================================================ Initialise local variables
    proshade_double *avgMat, *hlpMat, *eulA, *eulB, *eulG, *uAndV;
    proshade_unsign noPoints                          = static_cast< proshade_unsign > ( std::pow( ( ( peakSize * 2 ) + 1 ), 3 ) * 4 );
    proshade_double matWeight;
    
    //================================================ Allocate the required memory
    allocatePeakOptimisationMemory ( avgMat, hlpMat, eulA, eulB, eulG, uAndV );
    
    //================================================ For each peak
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( pointVec->size() ); iter++ )
    {
        //============================================ Re-set loop
        for ( proshade_unsign i = 0; i < 9; i++ ) { avgMat[i] = 0.0; }
        matWeight                                     = 0.0;
        
        //============================================ For each neighbout and the point itself
        for ( proshade_unsign it = 0; it < noPoints; it += 4 )
        {
            //======================================== Get the Euler angles
            ProSHADE_internal_maths::getEulerZXZFromSOFTPosition ( band, static_cast< proshade_signed > ( pointVec->at(iter)[it+0] ), static_cast< proshade_signed > ( pointVec->at(iter)[it+1] ), static_cast< proshade_signed > ( pointVec->at(iter)[it+2] ), eulA, eulB, eulG );
            
            //======================================== Convert Euler angles to rotation matrix
            ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( *eulA, *eulB, *eulG, hlpMat );
            
            //======================================== Add the matrix to the sum
            for ( proshade_unsign i = 0; i < 9; i++ ) { avgMat[i] += hlpMat[i] * pointVec->at(iter)[it+3]; }
            
            //======================================== Find matrix weight
            matWeight                                += pointVec->at(iter)[it+3];
        }
        
        //============================================ Normalise weighted sum matrix by sum of weights
        for ( proshade_unsign i = 0; i < 9; i++ ) { avgMat[i] /= matWeight; }
        
        //============================================  Decompose the average matrix using SVD
        ProSHADE_internal_maths::realMatrixSVDUandVOnly ( avgMat, 3, uAndV, false );
        const FloatingPoint< proshade_double > lhs ( uAndV[0] ), rhs ( -777.7 );
        if ( lhs.AlmostEquals ( rhs ) )
        {
            //======================================== SVD Failed. Just use the central value
            ProSHADE_internal_maths::getEulerZXZFromSOFTPosition ( band, static_cast< proshade_signed > ( pointVec->at(iter)[0] ), static_cast< proshade_signed > ( pointVec->at(iter)[1] ), static_cast< proshade_signed > ( pointVec->at(iter)[2] ), eulA, eulB, eulG );
            matWeight                                 = pointVec->at(iter)[3];
            pointVec->at(iter)                        = new proshade_double [4];
            ProSHADE_internal_misc::checkMemoryAllocation ( pointVec->at(iter), __FILE__, __LINE__, __func__ );
            pointVec->at(iter)[0]                     = *eulA;
            pointVec->at(iter)[1]                     = *eulB;
            pointVec->at(iter)[2]                     = *eulG;
            pointVec->at(iter)[3]                     = matWeight;
            
            continue ;
        }
        
        //============================================ SVD Succeeded. Compute U * V^T
        for ( proshade_unsign i = 0; i < 9; i++ ) { avgMat[i] = 0.0; }
        ProSHADE_internal_maths::multiplyTwoSquareMatrices ( uAndV, uAndV+9, avgMat, 3 );
        
        //============================================ Convert to Euler
        ProSHADE_internal_maths::getEulerZXZFromRotMatrix ( avgMat, eulA, eulB, eulG );
        
        //============================================ Save over input
        matWeight                                     = pointVec->at(iter)[3];
        delete[] pointVec->at(iter);
        pointVec->at(iter)                            = new proshade_double [4];
        ProSHADE_internal_misc::checkMemoryAllocation ( pointVec->at(iter), __FILE__, __LINE__, __func__ );
        pointVec->at(iter)[0]                         = *eulA;
        pointVec->at(iter)[1]                         = *eulB;
        pointVec->at(iter)[2]                         = *eulG;
        pointVec->at(iter)[3]                         = matWeight;
    }
    
    //================================================ Release memory
    releasePeakOptimisationMemory                     ( avgMat, hlpMat, eulA, eulB, eulG, uAndV );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds peaks in the 3D map using the "naive" approach.
 
    This function uses the "naive" approach (used in early versions of ProSHADE) to find all significant peaks in a map. To do this,
    it firstly locates all map points which have higher value than all their neighbours in all directions. It also computes the median
    and IQR of all non-higher points and it then uses the median + x * IQR threshold (x is the noIQRs parameter) to remove all map points
    with value under this thereshold. Finally, it optimises all the remaining values using the weighted average of all the neighbour
    points. The final output then is a vector with a single entry for each passing peak; this array has 4 values, the alpha, beta and gamma
    Euler angle values and the maximum peak heigh.
 
    \param[in] map Pointer to 1D array holding the 3D map value in which the peaks are to be found. Map must be cube!
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] peakSize The number of neighbouring points in single direction which should be considered as neighbours.
    \param[in] noIQRs The number of IQRs from the median to determine minimal peak height.
    \param[out] X Vector of located peaks with pointers to arrays of 4 values: alpha, beta, gamma and peak heighht.
 */
std::vector< proshade_double* > ProSHADE_internal_peakSearch::getAllPeaksNaive ( proshade_complex* map, proshade_unsign dim, proshade_signed peakSize, proshade_double noIQRs )
{
    //================================================ Find all indices with higher value than all neighbours
    std::vector< proshade_double* > allHigherIndices;
    proshade_double* nonPeakMedianIQR                 = new  proshade_double[2];
    ProSHADE_internal_misc::checkMemoryAllocation ( nonPeakMedianIQR, __FILE__, __LINE__, __func__ );
    allHigherIndices                                  = findAllPointsAboveNeighbours ( map, dim, peakSize, nonPeakMedianIQR );
    
    //================================================ Remove all indices with too small height
    pointsAboveNeighboursRemoveSmallHeight            ( &allHigherIndices, nonPeakMedianIQR, noIQRs );
    
    //================================================ Optimise the peaks using the neighbour values
    optimisePeakPositions                             ( &allHigherIndices, peakSize, dim/2 );
    
    //================================================ Release memory
    delete[]  nonPeakMedianIQR;
    
    //================================================ Done
    return                                            ( allHigherIndices );
    
}

/*! \brief This function finds the highest peaks optimised Euler angles using the "naive" approach.
 
    This function uses the same "naive" approach as discussed in the previous function. It firstly gets the full list of all
    "naive" peaks, it then optimises all such peaks and finds the largest (i.e. highest height) peak. For this single best
    peak, it finally returns the Euler angle values.
 
    \param[in] map Pointer to 1D array holding the 3D map value in which the peaks are to be found. Map must be cube!
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] noIQRs The number of IQRs from the median to determine minimal peak height.
    \param[in] eulA Pointer to where the Euler alpha angle value will be saved.
    \param[in] eulB Pointer to where the Euler beta angle value will be saved.
    \param[in] eulG Pointer to where the Euler gamma angle value will be saved.
    \param[in] settings The ProSHADE_settings object containing all the values required for making decisions.
 */
void ProSHADE_internal_peakSearch::getBestPeakEulerAngsNaive ( proshade_complex* map, proshade_unsign dim, proshade_double* eulA, proshade_double* eulB, proshade_double* eulG, ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Looking for Euler angles of highest peak.", settings->messageShift );
    
    //================================================ Get all peaks
    std::vector< proshade_double* > allPeaks          = getAllPeaksNaive ( map, dim, static_cast< proshade_signed > ( settings->peakNeighbours ), settings->noIQRsFromMedianNaivePeak );
    
    //================================================ Report progress
    std::stringstream hlpSSP;
    hlpSSP << "Found " << allPeaks.size() << " possible peaks.";
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, hlpSSP.str(), settings->messageShift );
    
    //================================================ Sanity check
    if ( allPeaks.size() == 0 )
    {
        *eulA                                         = 0.0;
        *eulB                                         = 0.0;
        *eulG                                         = 0.0;
        return ;
    }
    
    //================================================ Find the highest peak from the list
    proshade_double highestPeak                       = 0.0;
    proshade_unsign highestPeakIndex                  = 0;
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign>( allPeaks.size() ); iter++ )
    {
        if ( allPeaks.at(iter)[3] > highestPeak ) { highestPeak = allPeaks.at(iter)[3]; highestPeakIndex = iter; }
    }

    //================================================ Get Euler ZXZ for the highest peak
   *eulA                                              = allPeaks.at(highestPeakIndex)[0];
   *eulB                                              = allPeaks.at(highestPeakIndex)[1];
   *eulG                                              = allPeaks.at(highestPeakIndex)[2];
    
    //================================================ Release memory
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( allPeaks.size() ); iter++ )
    {
        delete[] allPeaks.at(iter);
    }
    
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "Optimal Euler angles are " << *eulA << " ; " << *eulB << " ; " << *eulG << " with peak height " << highestPeak;
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, hlpSS.str(), settings->messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allocates the memory required for smoothed Z score computation.
 
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] scoreOverVals Pointer to where the 1D cuts from the 3D map will be saved.
    \param[in] signals Pointer to where the results of the 1D peak searches will be saved.
    \param[in] filteredY Pointer to where the weighted smoothed scores will be saved.
    \param[in] avgFilter Pointer to where the median values will be saved.
    \param[in] stdFilter Pointer to where the IQR values will be saved.
    \param[in] subVec Pointer to where temporary  computations will be done.
    \param[in] medianIQR Pointer to simple array of 2 for returning results.
    \param[in] YZMap Pointer to where the X-axis peaks will be saved.
    \param[in] XZMap Pointer to where the Y-axis peaks will be saved.
    \param[in] XYMap Pointer to where the Z-axis peaks will be saved.
    \param[in] PeakMap Pointer to where the axis peak combinations will be saved.
    \param[in] smLag The size of the smoothing window.
 */
void ProSHADE_internal_peakSearch::allocateSmoothingZScoreMemory ( proshade_unsign dim, proshade_double*& scoreOverVals, proshade_signed*& signals, proshade_double*& filteredY, proshade_double*& avgFilter, proshade_double*& stdFilter, proshade_double*& subVec, proshade_double*& medianIQR, proshade_double*& YZMap, proshade_double*& XZMap, proshade_double*& XYMap, proshade_unsign smLag )
{
    //================================================ Allocate the memory
    signals                                           = new proshade_signed[dim];
    scoreOverVals                                     = new proshade_double[dim];
    filteredY                                         = new proshade_double[dim];
    avgFilter                                         = new proshade_double[dim];
    stdFilter                                         = new proshade_double[dim];
    subVec                                            = new proshade_double[smLag];
    medianIQR                                         = new proshade_double[2];
    YZMap                                             = new proshade_double[dim * dim * dim];
    XZMap                                             = new proshade_double[dim * dim * dim];
    XYMap                                             = new proshade_double[dim * dim * dim];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( signals,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( scoreOverVals, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( filteredY,     __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( avgFilter,     __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( stdFilter,     __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( subVec,        __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( medianIQR,     __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( YZMap,         __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( XZMap,         __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( XYMap,         __FILE__, __LINE__, __func__ );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function releases the memory required for smoothed Z score computation.
 
    \param[in] scoreOverVals Pointer to where the 1D cuts from the 3D map will be saved.
    \param[in] signals Pointer to where the results of the 1D peak searches will be saved.
    \param[in] filteredY Pointer to where the weighted smoothed scores will be saved.
    \param[in] avgFilter Pointer to where the median values will be saved.
    \param[in] stdFilter Pointer to where the IQR values will be saved.
    \param[in] subVec Pointer to where temporary  computations will be done.
    \param[in] medianIQR Pointer to simple array of 2 for returning results.
    \param[in] YZMap Pointer to where the X-axis peaks will be saved.
    \param[in] XZMap Pointer to where the Y-axis peaks will be saved.
    \param[in] XYMap Pointer to where the Z-axis peaks will be saved.
 */
void ProSHADE_internal_peakSearch::releaseSmoothingZScoreMemory ( proshade_double*& scoreOverVals, proshade_signed*& signals, proshade_double*& filteredY, proshade_double*& avgFilter, proshade_double*& stdFilter, proshade_double*& subVec, proshade_double*& medianIQR, proshade_double*& YZMap, proshade_double*& XZMap, proshade_double*& XYMap )
{
    //================================================ Release the memory
    delete[] scoreOverVals;
    delete[] signals;
    delete[] filteredY;
    delete[] avgFilter;
    delete[] stdFilter;
    delete[] subVec;
    delete[] medianIQR;
    delete[] YZMap;
    delete[] XZMap;
    delete[] XYMap;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the 1D peaks for a 1D input array and returns array of int's as signal.
 
    This function implements the smoothed Z score peak searching algorithm. It takes a number of previous values (assuming
    periodicity) and computes the median and IQR. It then checks if the current value is X IQR's from the median and reports
    peak if so. In this case, it also saves the value with decreased weight to make sure the sliding window for following indices
    will not have too high values due to previous peaks. If no peak is found, it saves the value as is. Finally, the function pre-
    computes the median and IQR for the next index and iterates.
 
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] smoothingLag The size of the smoothing window.
    \param[in] ZScoreThreshold The number of IQRs from median forming a peak.
    \param[in] signals Pointer to where the results of the 1D peak searches will be saved.
    \param[in] filteredY Pointer to where the weighted smoothed scores will be saved.
    \param[in] avgFilter Pointer to where the median values will be saved.
    \param[in] stdFilter Pointer to where the IQR values will be saved.
    \param[in] subVec Pointer to where temporary computations will be done.
    \param[in] medianIQR Pointer to simple array of 2 for returning results.
    \param[in] scoreOverVals Pointer to where the 1D cuts from the 3D map live.
 */
void ProSHADE_internal_peakSearch::getSmoothedZScore1D ( proshade_unsign dim, proshade_unsign smoothingLag, proshade_double ZScoreThreshold, proshade_signed* signals, proshade_double* filteredY, proshade_double* avgFilter, proshade_double* stdFilter, proshade_double* subVec, proshade_double* medianIQR, proshade_double* scoreOverVals )
{
    //================================================ Re-set the run
    for ( proshade_unsign i = 0; i < dim; i++ ) { signals[i] = 0; avgFilter[i] = 0.0; stdFilter[i] = 0.0; filteredY[i] = 0.0; }
    for ( proshade_unsign i = 0; i < smoothingLag; i++ ) { subVec[i] = scoreOverVals[i-smoothingLag+dim]; }
    ProSHADE_internal_maths::arrayMedianAndIQR ( subVec, smoothingLag, medianIQR );
    avgFilter[0]                                      = medianIQR[0];
    stdFilter[0]                                      = medianIQR[1];
    
    //================================================ Find peaks for 1D array
    for ( proshade_unsign i = 0; i < dim; i++)
    {
        //============================================ Is this peak?
        if ( std::abs ( scoreOverVals[i] - avgFilter[i]) > ZScoreThreshold * stdFilter[i] )
        {
            if ( scoreOverVals[i] > avgFilter[i] )
            {
                signals[i]                            = 1;
            }
            else
            {
                signals[i]                            = -1;
            }
            
            //======================================== Decrease influence for this window
            if ( i != 0 ) { filteredY[i] = 0.5 * scoreOverVals[i] + (1 - 0.5) * filteredY[i - 1]; }
            else          { filteredY[i] = 0.5 * scoreOverVals[i]; }
        }
        else
        {
            signals[i]                                = 0;
            filteredY[i]                              = scoreOverVals[i];
        }
        
        //============================================ Filters adjustments
        for ( proshade_signed subIt = 0; subIt < static_cast<proshade_signed> ( smoothingLag ); subIt++ )
        {
            if ( ( static_cast<proshade_signed>(i) + static_cast<proshade_signed>(subIt) - static_cast<proshade_signed>(smoothingLag) + 1 ) < 0 )
            {
                subVec[subIt]                         = scoreOverVals[( i + static_cast< proshade_unsign > ( subIt ) - smoothingLag + 1 ) + dim];
            }
            else
            {
                subVec[subIt]                         = filteredY[( i + static_cast< proshade_unsign > ( subIt ) - smoothingLag + 1 )];
            }
        }
        ProSHADE_internal_maths::arrayMedianAndIQR    ( subVec, smoothingLag, medianIQR );
        avgFilter[i+1]                                = medianIQR[0];
        stdFilter[i+1]                                = medianIQR[1];
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function runs the 1D smoothed Z score algorithm on all X-axis arrays as its inputs.
 
    This function iterates through all YZ map positions and at each, takes all the X values and subjects them to the
    1D smoothed Z score peak search. It also saves the resulting signals into a map, which therefore contains all the
    X-axis peaks of the map.
 
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] smoothingLag The size of the smoothing window.
    \param[in] ZScoreThreshold The number of IQRs from median forming a peak.
    \param[in] signals Pointer to where the results of the 1D peak searches will be saved.
    \param[in] filteredY Pointer to where the weighted smoothed scores will be saved.
    \param[in] avgFilter Pointer to where the median values will be saved.
    \param[in] stdFilter Pointer to where the IQR values will be saved.
    \param[in] subVec Pointer to where temporary  computations will be done.
    \param[in] medianIQR Pointer to simple array of 2 for returning results.
    \param[in] scoreOverVals Pointer to where the 1D cuts from the 3D map live.
    \param[in] map The map in which the peaks are to be found.
    \param[in] YZMap The map where the results will be saved.
 */
void ProSHADE_internal_peakSearch::getXAxisArraysSmoothedZScorePeaks ( proshade_unsign dim, proshade_unsign smoothingLag, proshade_double ZScoreThreshold, proshade_signed* signals, proshade_double* filteredY, proshade_double* avgFilter, proshade_double* stdFilter, proshade_double* subVec, proshade_double* medianIQR, proshade_double* scoreOverVals, proshade_complex* map, proshade_double* YZMap )
{
    //================================================ For each YZ point
    for ( proshade_unsign yIt = 0; yIt < dim; yIt++ )
    {
        for ( proshade_unsign zIt = 0; zIt < dim; zIt++ )
        {
            //======================================== Get all X values for this YZ position
            for ( proshade_unsign xIt = 0; xIt < dim; xIt++ )
            {
                scoreOverVals[xIt]                    = pow ( map[xIt * static_cast<proshade_unsign>(4 * pow ( ( dim / 2 ), 2 )) + yIt * dim + zIt][0], 2.0 ) +
                                                        pow ( map[xIt * static_cast<proshade_unsign>(4 * pow ( ( dim / 2 ), 2 )) + yIt * dim + zIt][1], 2.0 );
            }
            
            //======================================== Run 1D smoothed Z score algorithm
            getSmoothedZScore1D                       ( dim, smoothingLag, ZScoreThreshold, signals, filteredY, avgFilter, stdFilter, subVec, medianIQR, scoreOverVals );
            
            //======================================== Save signals to YZMap
            for ( proshade_unsign xIt = 0; xIt < dim; xIt++ )
            {
                YZMap[xIt * static_cast< proshade_unsign > (4 * pow ( ( dim / 2 ), 2 )) + yIt * dim + zIt] = static_cast< proshade_double > ( signals[xIt] );
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function runs the 1D smoothed Z score algorithm on all Y-axis arrays as its inputs.
 
    This function iterates through all YZ map positions and at each, takes all the Y values and subjects them to the
    1D smoothed Z score peak search. It also saves the resulting signals into a map, which therefore contains all the
    Y-axis peaks of the map.
 
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] smoothingLag The size of the smoothing window.
    \param[in] ZScoreThreshold The number of IQRs from median forming a peak.
    \param[in] signals Pointer to where the results of the 1D peak searches will be saved.
    \param[in] filteredY Pointer to where the weighted smoothed scores will be saved.
    \param[in] avgFilter Pointer to where the median values will be saved.
    \param[in] stdFilter Pointer to where the IQR values will be saved.
    \param[in] subVec Pointer to where temporary  computations will be done.
    \param[in] medianIQR Pointer to simple array of 2 for returning results.
    \param[in] scoreOverVals Pointer to where the 1D cuts from the 3D map live.
    \param[in] map The map in which the peaks are to be found.
    \param[in] XZMap The map where the results will be saved.
 */
void ProSHADE_internal_peakSearch::getYAxisArraysSmoothedZScorePeaks ( proshade_unsign dim, proshade_unsign smoothingLag, proshade_double ZScoreThreshold, proshade_signed* signals, proshade_double* filteredY, proshade_double* avgFilter, proshade_double* stdFilter, proshade_double* subVec, proshade_double* medianIQR, proshade_double* scoreOverVals, proshade_complex* map, proshade_double* XZMap )
{
    //================================================ For each XZ point
    for ( proshade_unsign xIt = 0; xIt < dim; xIt++ )
    {
        for ( proshade_unsign zIt = 0; zIt < dim; zIt++ )
        {
            //======================================== Get all Y values for this XZ position
            for ( proshade_unsign yIt = 0; yIt < dim; yIt++ )
            {
                scoreOverVals[yIt]                    = pow ( map[xIt * static_cast<proshade_unsign>(4 * pow ( ( dim / 2 ), 2 )) + yIt * dim + zIt][0], 2.0 ) +
                                                        pow ( map[xIt * static_cast<proshade_unsign>(4 * pow ( ( dim / 2 ), 2 )) + yIt * dim + zIt][1], 2.0 );
            }
            
            //======================================== Run 1D smoothed Z score algorithm
            getSmoothedZScore1D                       ( dim, smoothingLag, ZScoreThreshold, signals, filteredY, avgFilter, stdFilter, subVec, medianIQR, scoreOverVals );
            
            //======================================== Save signals to YZMap
            for ( proshade_unsign yIt = 0; yIt < dim; yIt++ )
            {
                XZMap[xIt * static_cast< proshade_unsign > (4 * pow ( ( dim / 2 ), 2 )) + yIt * dim + zIt] = static_cast< proshade_double > ( signals[xIt] );
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function runs the 1D smoothed Z score algorithm on all Z-axis arrays as its inputs.
 
    This function iterates through all XY map positions and at each, takes all the Z values and subjects them to the
    1D smoothed Z score peak search. It also saves the resulting signals into a map, which therefore contains all the
    Z-axis peaks of the map.
 
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] smoothingLag The size of the smoothing window.
    \param[in] ZScoreThreshold The number of IQRs from median forming a peak.
    \param[in] signals Pointer to where the results of the 1D peak searches will be saved.
    \param[in] filteredY Pointer to where the weighted smoothed scores will be saved.
    \param[in] avgFilter Pointer to where the median values will be saved.
    \param[in] stdFilter Pointer to where the IQR values will be saved.
    \param[in] subVec Pointer to where temporary  computations will be done.
    \param[in] medianIQR Pointer to simple array of 2 for returning results.
    \param[in] scoreOverVals Pointer to where the 1D cuts from the 3D map live.
    \param[in] map The map in which the peaks are to be found.
    \param[in] XYMap The map where the results will be saved.
 */
void ProSHADE_internal_peakSearch::getZAxisArraysSmoothedZScorePeaks ( proshade_unsign dim, proshade_unsign smoothingLag, proshade_double ZScoreThreshold, proshade_signed* signals, proshade_double* filteredY, proshade_double* avgFilter, proshade_double* stdFilter, proshade_double* subVec, proshade_double* medianIQR, proshade_double* scoreOverVals, proshade_complex* map, proshade_double* XYMap )
{
    //================================================ For each XZ point
    for ( proshade_unsign xIt = 0; xIt < dim; xIt++ )
    {
        for ( proshade_unsign yIt = 0; yIt < dim; yIt++ )
        {
            //======================================== Get all Y values for this XZ position
            for ( proshade_unsign zIt = 0; zIt < dim; zIt++ )
            {
                scoreOverVals[yIt]                    = pow ( map[xIt * static_cast<proshade_unsign>(4 * pow ( ( dim / 2 ), 2 )) + yIt * dim + zIt][0], 2.0 ) +
                                                        pow ( map[xIt * static_cast<proshade_unsign>(4 * pow ( ( dim / 2 ), 2 )) + yIt * dim + zIt][1], 2.0 );
            }
            
            //======================================== Run 1D smoothed Z score algorithm
            getSmoothedZScore1D                       ( dim, smoothingLag, ZScoreThreshold, signals, filteredY, avgFilter, stdFilter, subVec, medianIQR, scoreOverVals );
            
            //======================================== Save signals to YZMap
            for ( proshade_unsign zIt = 0; zIt < dim; zIt++ )
            {
                XYMap[xIt * static_cast< proshade_unsign > (4 * pow ( ( dim / 2 ), 2 )) + yIt * dim + zIt] = static_cast< proshade_double > ( signals[xIt] );
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This is a support function for the Z-score peak detection. It is currently not being used.
 
    This function should not really be present here, it is a support function for smoothed Z-score peak detection, which is currently not being used.
 
    \param[in] YZMap The map containing peaks detected along the X-axis dimension using the smoothed Z score method.
    \param[in] XZMap The map containing peaks detected along the Y-axis dimension using the smoothed Z score method.
    \param[in] XYMap The map containing peaks detected along the z-axis dimension using the smoothed Z score method.
    \param[in] visitedMap The map with 0 if point not yet considered and 1 if already done.
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] x The X-axis position of the currently considered point.
    \param[in] y The Y-axis position of the currently considered point.
    \param[in] z The Z-axis position of the currently considered point.
    \param[in] retVals A vector to which all neighbour indices will be added.
 */

void ProSHADE_internal_peakSearch::findAllPointNeighbours ( proshade_double* YZMap, proshade_double* XZMap, proshade_double* XYMap, proshade_unsign* visitedMap, proshade_signed dim, proshade_signed x, proshade_signed y, proshade_signed z, std::vector< proshade_unsign >* retVals )
{
    //================================================ Initialise local variables
    proshade_signed newIter                           = 0;
    proshade_signed peakX, peakY, peakZ;
    proshade_signed xDim                              = static_cast<proshade_signed>(4 * pow ( ( dim / 2 ), 2 ));
    
    //================================================ Iterate through all neighbours
    for ( proshade_signed xCh = -1; xCh <= +1; xCh++ )
    {
        for ( proshade_signed yCh = -1; yCh <= +1; yCh++ )
        {
            for ( proshade_signed zCh = -1; zCh <= +1; zCh++ )
            {
                //==================================== Do not use this point
                if ( ( xCh == 0 ) && ( yCh == 0 ) && ( zCh == 0 ) ) { continue; }
                
                //==================================== Find the nieghbout peak indices (with periodicity)
                peakX                                 = x + xCh; if ( peakX >= dim ) { peakX -= dim; }; if ( peakX < 0 ) { peakX += dim; }
                peakY                                 = y + yCh; if ( peakY >= dim ) { peakY -= dim; }; if ( peakY < 0 ) { peakY += dim; }
                peakZ                                 = z + zCh; if ( peakZ >= dim ) { peakZ -= dim; }; if ( peakZ < 0 ) { peakZ += dim; }
                newIter                               = peakX * xDim + peakY * dim + peakZ;
                
                //==================================== If already visited, next point
                if ( visitedMap[newIter] == 1 ) { continue; }
                else                            { visitedMap[newIter] = 1; }
                
                //==================================== If not peak, next point
                const FloatingPoint< proshade_double > lhs ( YZMap[newIter] + XZMap[newIter] + XYMap[newIter] ), rhs ( 3.0 );
                if ( lhs.AlmostEquals ( rhs ) ) { continue; }
                
                //==================================== This is a valid neighbour! Save
                ProSHADE_internal_misc::addToUnsignVector ( retVals, static_cast< proshade_unsign > ( newIter ) );
                
                //==================================== ... and recurse!
                findAllPointNeighbours                ( YZMap, XZMap, XYMap, visitedMap, dim, peakX, peakY, peakZ, retVals );
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function combines the three Z score maps, locates individual islands and returns a vector of highest point indices for each such island.
 
    This function considers each point common to all the three input maps and checks if this point has a signal along all three axes
    (i.e. is a peak in terms of the 3D smoothed Z score). If so, it finds all its neighbours and from this set (island), it finds the
    one point with the highest height. It then saves the index of this point and proceeds to another common point in the three input
    maps.
 
    \param[in] map Pointer to 1D array holding the 3D map value in which the peaks are to be found. Map must be cube!
    \param[in] YZMap The map containing peaks detected along the X-axis dimension using the smoothed Z score method.
    \param[in] XZMap The map containing peaks detected along the Y-axis dimension using the smoothed Z score method.
    \param[in] XYMap The map containing peaks detected along the z-axis dimension using the smoothed Z score method.
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] allIslandBests A vector which will be filled with the indices of the highest height point indices for each island.
 */
void ProSHADE_internal_peakSearch::findAllDisconnectedIslands ( proshade_complex* map, proshade_double* YZMap, proshade_double* XZMap, proshade_double* XYMap, proshade_unsign dim, std::vector< proshade_unsign >* allIslandBests )
{
    //================================================ Initialise local variables
    std::vector< proshade_unsign > ret;
    std::vector< proshade_unsign > thisIsland;
    proshade_unsign* visitedMap                       = new proshade_unsign[dim*dim*dim];
    ProSHADE_internal_misc::checkMemoryAllocation     ( visitedMap, __FILE__, __LINE__, __func__ );
    for ( proshade_unsign i = 0; i < ( dim * dim * dim ); i++ ) { visitedMap[i] = 0; }
    proshade_unsign index, maxIndex;
    proshade_double maxHeight;
    
    //================================================ For each point
    for ( proshade_unsign xIt = 0; xIt < dim; xIt++ )
    {
        for ( proshade_unsign yIt = 0; yIt < dim; yIt++ )
        {
            for ( proshade_unsign zIt = 0; zIt < dim; zIt++ )
            {
                //==================================== Find index
                index                                 = xIt * static_cast<proshade_unsign>(4 * pow ( ( dim / 2 ), 2 )) + yIt * dim + zIt;
                
                //==================================== If already visited, next point
                if ( visitedMap[index] == 1 ) { continue; }
                else                          { visitedMap[index] = 1; }
                
                //==================================== If not peak, next point
                const FloatingPoint< proshade_double > lhs ( YZMap[index] + XZMap[index] + XYMap[index] ), rhs ( 3.0 );
                if ( !lhs.AlmostEquals ( rhs ) ) { continue; }
                
                //==================================== This is a new island! Save this point
                thisIsland.clear                      ( );
                ProSHADE_internal_misc::addToUnsignVector ( &thisIsland, index );
                
                //==================================== ... and find all neighbours
                findAllPointNeighbours                ( YZMap, XZMap, XYMap, visitedMap, static_cast<proshade_signed>(dim), static_cast<proshade_signed>(xIt), static_cast<proshade_signed>(yIt), static_cast<proshade_signed>(zIt), &thisIsland );
                
                //==================================== Find the largest of the island peaks
                maxHeight                             = 0.0;
                maxIndex                              = 0;
                for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( thisIsland.size() ); iter++ )
                {
                    if ( ( pow( map[thisIsland.at(iter)][0], 2.0 ) + pow( map[thisIsland.at(iter)][1], 2.0) ) > maxHeight )
                    {
                        maxHeight                     = pow( map[thisIsland.at(iter)][0], 2.0 ) + pow( map[thisIsland.at(iter)][1], 2.0 );
                        maxIndex                      = thisIsland.at(iter);
                    }
                }
                
                //==================================== Save the results
                ProSHADE_internal_misc::addToUnsignVector ( allIslandBests, maxIndex );
            }
        }
    }
    
    //================================================ Release memory
    delete[] visitedMap;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function firstly determines the highest peak of all smoothed Z score islands and then returns this point as well as all its neighbours.
 
    This function firstly uses the findAllDisconnectedIslands function to combine the smoothed Z score maps and detect all discontinuous inslands
    in the combined map. Then, it uses the highest point in each island, finds all its neighbours as well as its X, Y and Z position and returns
    a vector containing all these values in a proshade_double pointer, which it allocates.
 
    \param[in] map Pointer to 1D array holding the 3D map value in which the peaks are to be found. Map must be cube!
    \param[in] YZMap The map containing peaks detected along the X-axis dimension using the smoothed Z score method.
    \param[in] XZMap The map containing peaks detected along the Y-axis dimension using the smoothed Z score method.
    \param[in] XYMap The map containing peaks detected along the z-axis dimension using the smoothed Z score method.
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] peakSize The number of neighbouring points in single direction which should be considered as neighbours.
    \param[in] allPeaksWithNeighbours A vector which will have the x, y, and z positions and heights of all peaks and all their neighbours.
 */
void ProSHADE_internal_peakSearch::findAllSmoothedZScorePeaksWithNeighbours ( proshade_complex* map, proshade_double* YZMap, proshade_double* XZMap, proshade_double* XYMap, proshade_signed dim, proshade_signed peakSize, std::vector< proshade_double* >* allPeaksWithNeighbours )
{
    //================================================ Initialise local variables
    proshade_signed x, y, z, peakX, peakY, peakZ, newIter, ptrIter;
    proshade_unsign noNeighbours                      = static_cast< proshade_unsign > ( std::pow( ( ( peakSize * 2 ) + 1 ), 3 ) * 4 );
    
    //================================================ Find all islands and their best representative
    std::vector < proshade_unsign > allPeaksRepresentatives;
    findAllDisconnectedIslands                        ( map, YZMap, XZMap, XYMap, static_cast< proshade_unsign > ( dim ), &allPeaksRepresentatives );
    
    //================================================ For each peak
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( allPeaksRepresentatives.size() ); iter++ )
    {
        //============================================ Reset loop
        ptrIter                                       = 0;
        
        //============================================ Determine x, y and z
        z                                             = ( static_cast< proshade_signed > ( allPeaksRepresentatives.at(iter) ) % (dim*dim) ) % dim;
        y                                             = ( static_cast< proshade_signed > ( allPeaksRepresentatives.at(iter) ) - z ) % (dim*dim) / dim;
        x                                             = ( static_cast< proshade_signed > ( allPeaksRepresentatives.at(iter) ) - z - ( y * dim ) ) / (dim*dim);
        
        //============================================ Allocate the memory
        proshade_double* neighbours                   = new proshade_double[noNeighbours];
        ProSHADE_internal_misc::checkMemoryAllocation ( neighbours, __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::addToDblPtrVector     ( allPeaksWithNeighbours, neighbours );
        
        //============================================ Find all neighbours
        for ( proshade_signed xCh = -peakSize; xCh <= +peakSize; xCh++ )
        {
            for ( proshade_signed yCh = -peakSize; yCh <= +peakSize; yCh++ )
            {
                for ( proshade_signed zCh = -peakSize; zCh <= +peakSize; zCh++ )
                {
                    //================================ Find the nieghbout peak indices (with periodicity)
                    peakX                             = x + xCh; if ( peakX >= dim ) { peakX -= dim; }; if ( peakX < 0 ) { peakX += dim; }
                    peakY                             = y + yCh; if ( peakY >= dim ) { peakY -= dim; }; if ( peakY < 0 ) { peakY += dim; }
                    peakZ                             = z + zCh; if ( peakZ >= dim ) { peakZ -= dim; }; if ( peakZ < 0 ) { peakZ += dim; }
                    newIter                           = peakX * (dim*dim) + peakY * dim + peakZ;
                    
                    //================================ Save neighbour values for optimisation of peaks later
                    allPeaksWithNeighbours->at(iter)[ptrIter]   = static_cast<proshade_double> ( peakX );
                    allPeaksWithNeighbours->at(iter)[ptrIter+1] = static_cast<proshade_double> ( peakY );
                    allPeaksWithNeighbours->at(iter)[ptrIter+2] = static_cast<proshade_double> ( peakZ );
                    allPeaksWithNeighbours->at(iter)[ptrIter+3] = pow ( map[newIter][0], 2.0 ) + pow( map[newIter][1], 2.0 );
                    ptrIter                          += 4;
                }
            }
        }
        
        //============================================ Release memory
        delete[] neighbours;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds peaks in the 3D map using the smoothed Z score approach.
 
    This function implements the smoothed Z-score peak detection. It is currently not being used and is most likely severly bugged.
 
    \param[in] map Pointer to 1D array holding the 3D map values in which the peaks are to be found. Map must be cube!
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] peakSize The number of neighbouring points in single direction which should be considered as neighbours.
    \param[in] noIQRs The number of IQRs from the median to determine minimal peak height.
    \param[in] peakSize The number of neighbouring points in single direction which should be considered as neighbours.
    \param[out] X Vector of located peaks with pointers to arrays of 5 values: x, y, z, angle and peak heighht.
 */
std::vector< proshade_double* > ProSHADE_internal_peakSearch::getAllPeaksSmoothedZ ( proshade_complex* map, proshade_unsign dim, proshade_double smoothingFraction, proshade_double noIQRs, proshade_signed peakSize )
{
    //================================================ Initialise local variables
    std::vector< proshade_double* > allHigherIndices;
    proshade_unsign smoothingLag                      = static_cast< proshade_unsign > ( std::floor ( smoothingFraction * static_cast< proshade_double > ( dim ) ) );
    proshade_double ZScoreThreshold                   = noIQRs;
    proshade_signed *signals;
    proshade_double *scoreOverVals, *filteredY, *avgFilter, *stdFilter, *subVec, *medianIQR, *YZMap, *XZMap, *XYMap;
    
    //================================================ Sanity check
    if ( dim <= smoothingLag + 2)
    {
        // throw error
    }
    
    //================================================ Allocate required memory
    allocateSmoothingZScoreMemory                     ( dim, scoreOverVals, signals, filteredY, avgFilter, stdFilter, subVec, medianIQR, YZMap, XZMap, XYMap, smoothingLag );
    
    //================================================ Get smoothed Z score peaks for X-axis arrays
    getXAxisArraysSmoothedZScorePeaks                 ( dim, smoothingLag, ZScoreThreshold, signals, filteredY, avgFilter, stdFilter, subVec, medianIQR, scoreOverVals, map, YZMap );
    
    //================================================ Get smoothed Z score peaks for Y-axis arrays
    getYAxisArraysSmoothedZScorePeaks                 ( dim, smoothingLag, ZScoreThreshold, signals, filteredY, avgFilter, stdFilter, subVec, medianIQR, scoreOverVals, map, XZMap );
    
    //================================================ Get smoothed Z score peaks for Y-axis arrays
    getZAxisArraysSmoothedZScorePeaks                 ( dim, smoothingLag, ZScoreThreshold, signals, filteredY, avgFilter, stdFilter, subVec, medianIQR, scoreOverVals, map, XYMap );
    
    //================================================ Put all dimension peaks together and detect disconnected islands
    findAllSmoothedZScorePeaksWithNeighbours          ( map, YZMap, XZMap, XYMap, static_cast< proshade_signed > ( dim ), peakSize, &allHigherIndices );
    
    //================================================ Optimise the peaks using the neighbour values
    optimisePeakPositions                             ( &allHigherIndices, peakSize, dim/2 );
    
    //================================================ Release the required memory
    releaseSmoothingZScoreMemory                      ( scoreOverVals, signals, filteredY, avgFilter, stdFilter, subVec, medianIQR, YZMap, XZMap, XYMap );
    
    //================================================ Done
    return                                            ( allHigherIndices );
    
}

/*! \brief This function finds the highest peaks optimised Euler angles using the smoothed Z score approach.
 
    This function uses the best peak detected using the smoothed Z-score peak detection to compute its Euler angles. It is currently not being used and is most likely severly bugged.
 
    \param[in] map Pointer to 1D array holding the 3D map value in which the peaks are to be found. Map must be cube!
    \param[in] dim The size of one dimension of the map (assuming cube map).
    \param[in] peakSize The number of neighbouring points in single direction which should be considered as neighbours.
    \param[in] noIQRs The number of IQRs from the median to determine minimal peak height.
    \param[in] eulA Pointer to where the Euler alpha angle value will be saved.
    \param[in] eulB Pointer to where the Euler beta angle value will be saved.
    \param[in] eulG Pointer to where the Euler gamma angle value will be saved.
 */
void ProSHADE_internal_peakSearch::getBestPeakEulerAngsSmoothedZ ( proshade_complex* map, proshade_unsign dim, proshade_double smoothingFraction, proshade_double noIQRs, proshade_signed peakSize, proshade_double* eulA, proshade_double* eulB, proshade_double* eulG )
{
    //================================================ Get all peaks
    std::vector< proshade_double* > allPeaks          = getAllPeaksSmoothedZ ( map, dim, smoothingFraction, noIQRs, peakSize );
    
    //================================================ Find the highest peak from the list
    proshade_double highestPeak                       = 0.0;
    proshade_unsign highestPeakIndex                  = 0;
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign>( allPeaks.size() ); iter++ )
    {
        if ( allPeaks.at(iter)[4] > highestPeak ) { highestPeak = allPeaks.at(iter)[4]; highestPeakIndex = iter; }
    }
    
    //================================================ Get Euler ZXZ from Angle-axis
    proshade_double* rotMat                           = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( rotMat, allPeaks.at(highestPeakIndex)[0], allPeaks.at(highestPeakIndex)[1], allPeaks.at(highestPeakIndex)[2], allPeaks.at(highestPeakIndex)[3] );
    ProSHADE_internal_maths::getEulerZXZFromRotMatrix ( rotMat, eulA, eulB, eulG );
    
    //================================================ Release memory
    delete[] rotMat;
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( allPeaks.size() ); iter++ )
    {
        delete[] allPeaks.at(iter);
    }
    
    //================================================ Done
    return ;
    
}
