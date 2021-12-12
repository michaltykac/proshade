/*! \file ProSHADE_symmetry.cpp
    \brief This source file contains all the functions required to detect symmetry axes and types from the inverse SOFT map.
 
    The functions in this source file are required to allow detection of symmetry axes and symmetries from the inverse SOFT map. The currect functionality can detect C, D, T, O and I symmetries
    with the C and D symmetries having their fold automatically detected as well.
 
    Copyright by Michal Tykac and individual contributors. All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
 
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.6.1
    \date      AUG 2021
 */

//==================================================== ProSHADE
#include "ProSHADE_symmetry.hpp"

//==================================================== Local functions prototypes
proshade_double                               determinePeakThreshold      ( std::vector < proshade_double > inArr, proshade_double noIQRsFromMedian );
bool                                          sortProSHADESymmetryByPeak  ( proshade_double* a, proshade_double* b );
std::vector < std::pair< proshade_unsign, proshade_unsign > > findBestIcosDihedralPair    ( std::vector< proshade_double* >* CSymList, proshade_double minPeakHeight, proshade_double axErr );
std::pair< proshade_unsign, proshade_unsign > findBestOctaDihedralPair    ( std::vector< proshade_double* >* CSymList, proshade_double minPeakHeight, proshade_double axErr );
std::pair< proshade_unsign, proshade_unsign > findBestTetraDihedralPair   ( std::vector< proshade_double* >* CSymList, proshade_double minPeakHeight, proshade_double axErr );

/*! \brief This function computes the self-rotation function for this structure.
 
    This function assumes that the spherical harmonics have been computed for a data object. It can be then called on this
    object and it will proceed to compute the E matrices for this object against itself. From these "self E matrices", the
    function will generate the SO(3) transform coefficients and finally it will invert transform these coefficients back,
    thus getting the self-rotation function.
 
    \param[in] settings A pointer to settings class containing all the information required for map self-rotation function computation.
 */
void ProSHADE_internal_data::ProSHADE_data::computeRotationFunction ( ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting self-rotation function computation.", settings->messageShift );
    
    //================================================ Compute un-weighted E matrices and their weights
    ProSHADE_internal_distances::computeEMatrices     ( this, this, settings );
    
    //================================================ Normalise E matrices by the magnitudes
    ProSHADE_internal_distances::normaliseEMatrices   ( this, this, settings );
    
    //================================================ Generate SO(3) coefficients
    ProSHADE_internal_distances::generateSO3CoeffsFromEMatrices ( this, this, settings );
    
    //================================================ Compute the inverse SO(3) Fourier Transform (SOFT) on the newly computed coefficients
    ProSHADE_internal_distances::computeInverseSOFTTransform ( this, this, settings );
    
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Self-rotation function obtained.", settings->messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes a vector of values and determines the threshold for removing noise from it.

    \param[in] inArr A vector of values for which the threshold is to be determined.
    \param[out] ret The threshold.
 */
proshade_double determinePeakThreshold ( std::vector < proshade_double > inArr, proshade_double noIQRsFromMedian )
{
    //================================================ Initialise variables
    proshade_double ret                               = 0.0;
    proshade_double rmsd                              = 0.0;
    size_t vecSize                                    = inArr.size();
    
    //================================================ Deal with low number of input cases
    if ( vecSize == 0 )                               { return ( ret ); }                                                                                                           // Return 0
    if ( vecSize <= 4 )                               { ret = std::accumulate ( inArr.begin(), inArr.end(), 0.0 ) / static_cast< proshade_double > ( vecSize ); return ( ret ); }   // Return mean
    
    //================================================ Deal with reasonable number in input cases
    else
    {
        //============================================ Find mean
        ret                                           = std::accumulate ( inArr.begin(), inArr.end(), 0.0 ) / static_cast< proshade_double > ( vecSize );
        
        //============================================ Get the RMS distance
        for ( size_t i = 0; i < vecSize; i++ )
        {
            rmsd                                     += std::pow ( ret - inArr.at(i), 2.0 );
        }
        rmsd                                          = std::sqrt ( rmsd );
        
        //============================================ Get the threshold
        ret                                           = ret + ( noIQRsFromMedian * rmsd );
    }
    
    //================================================ Sanity checks
    if ( ret > *( std::max_element ( inArr.begin(), inArr.end() ) ) )
    {
        ret                                           = *( std::max_element ( inArr.begin(), inArr.end() ) );
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function converts the self-rotation function of this structure to angle-axis representation.
 
    This function creates a set of concentric spheres in a spherical co-ordinates space, where the radius is the angle-axis representation angle and
    the lattitude and longitude angles are the angle-axis representation axis vector. I.e. each of the spheres contains all angle-axis representation
    axes for a single given angle.
 
    Then, it proceeds to interpolate the rotation function for each point in this space, thus effectivelly re-sampling the rotation function onto the required
    space.
 
    \param[in] settings A pointer to settings class containing all the information required for map self-rotation function computation.
 */
void ProSHADE_internal_data::ProSHADE_data::convertRotationFunction ( ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting self-rotation function conversion to angle-axis representation.", settings->messageShift );
    
    //================================================ Initialise variables
    proshade_double shellSpacing                      = ( 2.0 * M_PI ) / static_cast<proshade_double> ( this->maxShellBand ) * 2.0;
    std::vector< proshade_double > allPeakHeights;
    
    //================================================ Initialise the spheres
    for ( proshade_unsign spIt = 1; spIt < ( this->maxShellBand * 2 ); spIt++ )
    {
        this->sphereMappedRotFun.emplace_back         ( new ProSHADE_internal_spheres::ProSHADE_rotFun_sphere( static_cast<proshade_double> ( spIt ) * shellSpacing,
                                                                                                               shellSpacing,
                                                                                                               this->maxShellBand * 2,
                                                                                                               static_cast<proshade_double> ( spIt ) * shellSpacing,
                                                                                                               spIt - 1 ) );
    }

    //================================================ Interpolate the rotation function onto the spheres
    for ( proshade_unsign shIt = 0; shIt < static_cast<proshade_unsign> ( sphereMappedRotFun.size() ); shIt++ )
    {
        //============================================ Report progress
        std::stringstream hlpSS;
        hlpSS << "Interpolating sphere " << shIt << " ( radius: " << this->sphereMappedRotFun.at(shIt)->getRadius() << " ).";
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 3, hlpSS.str(), settings->messageShift );
        
        //============================================ Interpolate onto spheres
        this->sphereMappedRotFun.at(shIt)->interpolateSphereValues ( this->getInvSO3Coeffs ( ) );
    }
    
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Self-rotation function converted to spherical angle-axis space.", settings->messageShift );
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Started peak detection on the angle-axis spheres.", settings->messageShift );
    
    //================================================ Find all peaks in the sphere grids
    for ( proshade_unsign shIt = 0; shIt < static_cast<proshade_unsign> ( this->sphereMappedRotFun.size() ); shIt++ )
    {
        this->sphereMappedRotFun.at(shIt)->findAllPeaks  ( static_cast< proshade_signed > ( settings->peakNeighbours ), &allPeakHeights );
    }
    
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "Detected " << allPeakHeights.size() << " peaks with any height.";
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, hlpSS.str(), settings->messageShift );
    
    //================================================ Compute threshold for small peaks
    proshade_double peakThres                         = std::max ( settings->minSymPeak, determinePeakThreshold ( allPeakHeights, settings->noIQRsFromMedianNaivePeak ) );
    
    //================================================ Report progress
    std::stringstream hlpSS2;
    hlpSS2 << "From these peaks, decided the threshold will be " << peakThres << " peak height.";
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 4, hlpSS2.str(), settings->messageShift );

    //================================================ Remove too small peaks
    for ( proshade_unsign shIt = 0; shIt < static_cast<proshade_unsign> ( this->sphereMappedRotFun.size() ); shIt++ )
    {
        this->sphereMappedRotFun.at(shIt)->removeSmallPeaks  ( peakThres );
    }

    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 3, "Peaks detected for all spheres.", settings->messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function checks if a very similar symmetry is not already saved.
 
    This is a simple function comparing a single double array of 6 to a vector of these, returning whether the vector already
    contains a very similar entry to the rested one. If the new has better height, replacement will take place.
 
    \param[in] ret This is the variable where the tested array will be saved if passed. It is a vector of double[6] arrays with the following meaning: [0] = fold, [1] = x-axis, [2] = y-axis, [3] = z-axis, [4] = angle, [5] = average peak height.
    \param[in] sym This is a double array of 6 which is to be compared to all the vector entries.
    \param[in] simThres The threshold for dot product comparison similarity.
    \param[in] matchedPos Pointer to variable where the matched position (if any axis is matched) is saved, or -1 is written.
    \param[out] X Boolean value stating whether a similar entry has been found (true = it was, false = it was not).
 */
bool ProSHADE_internal_symmetry::isSymmetrySame ( std::vector< proshade_double* >* ret, proshade_double* sym, proshade_double simThres, proshade_signed* matchedPos )
{
    //================================================ Initialise variables
    proshade_double dotProduct                        = 0.0;
   *matchedPos                                        = -1;
    
    //================================================ Check
    for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( ret->size() ); symIt++ )
    {
        //============================================ Minor speed-up => only test for same folds
        const FloatingPoint< proshade_double > lhs ( ret->at(symIt)[0] ), rhs ( sym[0] );
        if ( lhs.AlmostEquals ( rhs ) )
        {
            //======================================== Is axis the same?
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &ret->at(symIt)[1], &ret->at(symIt)[2],
                                                                                                     &ret->at(symIt)[3], &sym[1], &sym[2], &sym[3] );
            if ( ( ( 1.0 > ( dotProduct - simThres ) ) && ( 1.0 < ( dotProduct + simThres ) ) ) || ( ( -1.0 > ( dotProduct - simThres ) ) && ( -1.0 < ( dotProduct + simThres ) ) ) )
            {
                //==================================== Matched. Save the index
               *matchedPos                            = static_cast< proshade_signed > ( symIt );
                
                //==================================== Does the already saved have higher height?
                if ( ret->at(symIt)[5] >= sym[5] ) { return ( true ); }
                
                //==================================== In this case, new is better than old - sort it out
                ret->at(symIt)[1]                     = sym[1];
                ret->at(symIt)[2]                     = sym[2];
                ret->at(symIt)[3]                     = sym[3];
                ret->at(symIt)[5]                     = sym[5];
                return                                ( true );
            }
        }
    }
    
    //================================================ Done - no matches found
    return                                            ( false );
    
}

/*! \brief This function checks if a very similar symmetry is not already saved.
 
    This is a simple function comparing a single double array of 6 to a vector of these, returning whether the vector already
    contains a very similar entry to the rested one. If the new has better height, replacement will take place.
 
    \param[in] ret This is the variable where the tested array will be saved if passed. It is a vector of double[6] arrays with the following meaning: [0] = fold, [1] = x-axis, [2] = y-axis, [3] = z-axis, [4] = angle, [5] = average peak height.
    \param[in] sym This is a double array of 6 which is to be compared to all the vector entries.
    \param[in] simThres The threshold for dot product comparison similarity.
    \param[in] matchedPos Pointer to variable where the matched position (if any axis is matched) is saved, or -1 is written.
    \param[in] fscVal Value to be used as FSC in case of a match.
    \param[out] X Boolean value stating whether a similar entry has been found (true = it was, false = it was not).
 */
bool ProSHADE_internal_symmetry::isSymmetrySame ( std::vector< proshade_double* >* ret, proshade_double* sym, proshade_double simThres, proshade_signed* matchedPos, proshade_double fscVal )
{
    //================================================ Initialise variables
    proshade_double dotProduct                        = 0.0;
   *matchedPos                                        = -1;
    
    //================================================ Check
    for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( ret->size() ); symIt++ )
    {
        //============================================ Minor speed-up => only test for same folds
        const FloatingPoint< proshade_double > lhs ( ret->at(symIt)[0] ), rhs ( sym[0] );
        if ( lhs.AlmostEquals ( rhs ) )
        {
            //======================================== Is axis the same?
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &ret->at(symIt)[1], &ret->at(symIt)[2],
                                                                                                     &ret->at(symIt)[3], &sym[1], &sym[2], &sym[3] );
            if ( ( ( 1.0 > ( dotProduct - simThres ) ) && ( 1.0 < ( dotProduct + simThres ) ) ) || ( ( -1.0 > ( dotProduct - simThres ) ) && ( -1.0 < ( dotProduct + simThres ) ) ) )
            {
                //==================================== Matched. Save the index
               *matchedPos                            = static_cast< proshade_signed > ( symIt );
                
                //==================================== Does the already saved have higher height?
                if ( ret->at(symIt)[5] >= sym[5] ) { return ( true ); }
                
                //==================================== In this case, new is better than old - sort it out
                ret->at(symIt)[1]                     = sym[1];
                ret->at(symIt)[2]                     = sym[2];
                ret->at(symIt)[3]                     = sym[3];
                ret->at(symIt)[5]                     = sym[5];
                ret->at(symIt)[6]                     = fscVal;
                return                                ( true );
            }
        }
    }
    
    //================================================ Done - no matches found
    return                                            ( false );
    
}

/*! \brief This function obtains a list of all D symmetries from already computed C  symmetries list.
 
    This function simply returns a vector of C symmetry pairs which are perpendicular to each other (and therefore form dihedral symmetry).
    The vector contains arrays of 12 double numbers with the following format: [0] = Fold of axis 1; [1] = X-axis of axis 1; [2] Y-axis of
    axis 1; [3] = Z-axis of axis 1; [4] = angle of axis 1; [5] = average peak height of axis 1; [6] = Fold of axis 2; [7] = X-axis of axis 2;
    [8] Y-axis of axis 2; [9] = Z-axis of axis 2; [10] = angle of axis 2; [11] = average peak height of axis 2. Note that the larger fold axis
    is listed first in this format.
 
    \param[in] settings A pointer to settings class containing all the information required for symmetry detection.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
 */
std::vector< proshade_double* > ProSHADE_internal_data::ProSHADE_data::getDihedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList )
{
    //================================================ Initialise variables
    std::vector< proshade_double* > ret;
    proshade_double dotProduct;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Starting D symmetry detection.", settings->messageShift );
    
    //================================================If not enough axes, just end here
    if ( CSymList->size() < 2 ) { return ( ret ); }
    
    //================================================ For each unique pair of axes
    for ( proshade_unsign ax1 = 0; ax1 < static_cast<proshade_unsign> ( CSymList->size() ); ax1++ )
    {
        //============================================ Ignore small axes
        const FloatingPoint< proshade_double > lhs1 ( CSymList->at(ax1)[5] ), rhs1 ( -999.9 );
        if ( ( CSymList->at(ax1)[5] < settings->minSymPeak ) && !( lhs1.AlmostEquals ( rhs1 ) ) ) { continue; }
        
        for ( proshade_unsign ax2 = 1; ax2 < static_cast<proshade_unsign> ( CSymList->size() ); ax2++ )
        {
            //======================================= Use unique pairs only
            if ( ax1 >= ax2 ) { continue; }
            
            //======================================== Ignore small axes
            const FloatingPoint< proshade_double > lhs2 ( CSymList->at(ax2)[5] ), rhs2 ( -999.9 );
            if ( ( CSymList->at(ax2)[5] < settings->minSymPeak ) && !( lhs2.AlmostEquals ( rhs2 ) ) ) { continue; }
            
            //======================================= Compute the dot product
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &CSymList->at(ax1)[1], &CSymList->at(ax1)[2],
                                                                                                     &CSymList->at(ax1)[3], &CSymList->at(ax2)[1],
                                                                                                     &CSymList->at(ax2)[2], &CSymList->at(ax2)[3] );
            
            //======================================== If close to zero, these two axes are perpendicular
            if ( std::abs( dotProduct ) < settings->axisErrTolerance )
            {
                //==================================== Save
                if ( CSymList->at(ax1)[0] >= CSymList->at(ax2)[0] )
                {
                    ProSHADE_internal_symmetry::saveDSymmetry ( &ret, CSymList, ax1, ax2 );
                    
                    std::vector< proshade_unsign > DSymInd;
                    ProSHADE_internal_misc::addToUnsignVector ( &DSymInd, ax1 );
                    ProSHADE_internal_misc::addToUnsignVector ( &DSymInd, ax2 );
                    ProSHADE_internal_misc::addToUnsignVectorVector ( &settings->allDetectedDAxes, DSymInd );
                    
                }
                else
                {
                    ProSHADE_internal_symmetry::saveDSymmetry ( &ret, CSymList, ax2, ax1 );
                    
                    std::vector< proshade_unsign > DSymInd;
                    ProSHADE_internal_misc::addToUnsignVector ( &DSymInd, ax2 );
                    ProSHADE_internal_misc::addToUnsignVector ( &DSymInd, ax1 );
                    ProSHADE_internal_misc::addToUnsignVectorVector ( &settings->allDetectedDAxes, DSymInd );
                }
            }
        }
    }
    
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "Detected " << ret.size() << " D symmetries.";
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, hlpSS.str(), settings->messageShift );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function saves a detected dihedral symmetry to the dihedral symmetries list.
 
    This function takes two C symmetry axes as supplied by the calling function and the list of the detected C symmetries. It then
    produces the saving structure for a dihedral symmetry formed by the two supplied axes and saves this structure to the supplied
    dihedral symmetry list vector - ret.
 
    \param[in] ret The vector of double pointers to which the symmetry is to be saved to.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] axisOne The index of the first C symmetry forming the dihedral symmetry.
    \param[in] axisTwo The index of the second C symmetry forming the dihedral symmetry.
 */
void ProSHADE_internal_symmetry::saveDSymmetry ( std::vector< proshade_double* >* ret, std::vector< proshade_double* >* CSymList, proshade_unsign axisOne, proshade_unsign axisTwo )
{
    //================================================ Allocate the memory
    proshade_double* hlpP                             = new proshade_double [14];
    ProSHADE_internal_misc::checkMemoryAllocation     ( hlpP, __FILE__, __LINE__, __func__ );
    
    //================================================ Set the axis and heights
    hlpP[0]                                           = CSymList->at(axisOne)[0];
    hlpP[1]                                           = CSymList->at(axisOne)[1];
    hlpP[2]                                           = CSymList->at(axisOne)[2];
    hlpP[3]                                           = CSymList->at(axisOne)[3];
    hlpP[4]                                           = CSymList->at(axisOne)[4];
    hlpP[5]                                           = CSymList->at(axisOne)[5];
    hlpP[6]                                           = CSymList->at(axisOne)[6];
    hlpP[7]                                           = CSymList->at(axisTwo)[0];
    hlpP[8]                                           = CSymList->at(axisTwo)[1];
    hlpP[9]                                           = CSymList->at(axisTwo)[2];
    hlpP[10]                                          = CSymList->at(axisTwo)[3];
    hlpP[11]                                          = CSymList->at(axisTwo)[4];
    hlpP[12]                                          = CSymList->at(axisTwo)[5];
    hlpP[13]                                          = CSymList->at(axisTwo)[6];
    
    //================================================ Save to ret
    ProSHADE_internal_misc::addToDblPtrVector         ( ret, hlpP );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function obtains a list of all T symmetry axes from the already computed C symmetries list.
 
    This function starts by checking if there are two C3 symmetries with the tetrahedral dihedral angle. If so, it proceeds to search for all seven symmetry axes
    expected to form a full tetrahedral symmetry. It then returns the list of found symmetries; if full tetrahedral symmetry was found, seven axes (four C3s and
    three C2s) are returned. If less than seven symmetries are returned, the procedure has failed and no tetrahedral symmetry was found.
 
    \param[in] settings A pointer to settings class containing all the information required for symmetry detection.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
 */
std::vector< proshade_double* > ProSHADE_internal_data::ProSHADE_data::getTetrahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList )
{
    //================================================ Initialise variables
    std::vector< proshade_double* > ret;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting T symmetry detection.", settings->messageShift );
    
    //================================================ Are the basic requirements for tetrahedral symmetry met?
    if ( ProSHADE_internal_symmetry::detectTetrahedralSymmetry ( CSymList, settings->axisErrTolerance, settings->minSymPeak ) )
    {
        //============================================ Search for all the symmetry axes
        ProSHADE_internal_symmetry::findTetra4C3s     ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->messageShift, settings->minSymPeak );
        if ( ret.size() != 4 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        
        ProSHADE_internal_symmetry::findTetra3C2s     ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->messageShift, settings->minSymPeak );
        if ( ret.size() != 7 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        else
        {
            for ( proshade_unsign csIt = 0; csIt < static_cast<proshade_unsign> ( CSymList->size() ); csIt++ )
            {
                for ( proshade_unsign retIt = 0; retIt < static_cast<proshade_unsign> ( ret.size() ); retIt++ )
                {
                    //======================================== Sort ret by fold
                    std::sort                                 ( ret.begin(), ret.end(), ProSHADE_internal_misc::sortSymInvFoldHlp );
                    
                    //======================================== Save indices
                    const FloatingPoint< proshade_double > lhs1 ( CSymList->at(csIt)[0] ), rhs1 ( ret.at(retIt)[0] );
                    const FloatingPoint< proshade_double > lhs2 ( CSymList->at(csIt)[1] ), rhs2 ( ret.at(retIt)[1] );
                    const FloatingPoint< proshade_double > lhs3 ( CSymList->at(csIt)[2] ), rhs3 ( ret.at(retIt)[2] );
                    const FloatingPoint< proshade_double > lhs4 ( CSymList->at(csIt)[3] ), rhs4 ( ret.at(retIt)[3] );
                    const FloatingPoint< proshade_double > lhs5 ( CSymList->at(csIt)[4] ), rhs5 ( ret.at(retIt)[4] );
                    const FloatingPoint< proshade_double > lhs6 ( CSymList->at(csIt)[5] ), rhs6 ( ret.at(retIt)[5] );
                    if ( ( lhs1.AlmostEquals ( rhs1 ) ) &&
                         ( lhs2.AlmostEquals ( rhs2 ) ) &&
                         ( lhs3.AlmostEquals ( rhs3 ) ) &&
                         ( lhs4.AlmostEquals ( rhs4 ) ) &&
                         ( lhs5.AlmostEquals ( rhs5 ) ) &&
                         ( lhs6.AlmostEquals ( rhs6 ) ) )
                    {
                        ProSHADE_internal_misc::addToUnsignVector ( &settings->allDetectedTAxes, csIt );
                    }
                }
            }
        }
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "T symmetry detection complete.", settings->messageShift );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function takes the list of C symmetries and decides whether basic requirements for tetrahedral symmetry are there.
 
    This function first finds all the C3 symmetries in the C symmetries list and then it checks all pais of such present C3s for have the angle
    between the pair equal to the dihedral angle of a tetrahedron ( acos(1/3) ). If a single such pair is detected, this is likely a tetrahedral
    symmetry and all other axes need to be located. Otherwise, false is returned.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
    \param[out] X Boolean value telling whether there are two C3 symmetries with tetrahedral dihhedral angle.
 */
bool ProSHADE_internal_symmetry::detectTetrahedralSymmetry ( std::vector< proshade_double* >* CSymList, proshade_double axErr, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C3List;
    proshade_double dotProduct;
    
    //================================================ Find all C3 symmetries
    for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ )
    {
        const FloatingPoint< proshade_double > lhs ( CSymList->at(cSym)[0] ), rhs ( 3.0 );
        if ( lhs.AlmostEquals ( rhs ) && CSymList->at(cSym)[5] >= minPeakHeight ) { ProSHADE_internal_misc::addToUnsignVector ( &C3List, cSym ); }
    }
    
    //================================================ For each unique pair of C3s
    for ( proshade_unsign c31 = 0; c31 < static_cast<proshade_unsign> ( C3List.size() ); c31++ )
    {
        for ( proshade_unsign c32 = 1; c32 < static_cast<proshade_unsign> ( C3List.size() ); c32++ )
        {
            //================================ Unique pairs only
            if ( c31 >= c32 ) { continue; }
            
            //========================================  Check the angle between the C3 axes
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &CSymList->at(C3List.at(c31))[1], &CSymList->at(C3List.at(c31))[2], &CSymList->at(C3List.at(c31))[3], &CSymList->at(C3List.at(c32))[1], &CSymList->at(C3List.at(c32))[2], &CSymList->at(C3List.at(c32))[3] );
            
            //================================ Is the angle approximately the dihedral angle
            if ( ( ( 1.0 / 3.0 ) > ( dotProduct - axErr ) ) && ( ( 1.0 / 3.0 ) < ( dotProduct + axErr ) ) )
            {
                return                                ( true );
            }
        }
    }
    
    //================================================ Done
    return                                            ( false );
    
}

/*! \brief This function takes the list of C symmetries and finds the 4 C3 symmetries with correct angles required for full tetrahedral symmetry.
 
    This function is specific to detecting the tetrahedral symmetry. It should be called once tetrahedral symmetry is suspected (by detecting its dihedral
    angles) and it needs to be fully described. This function specifically searches for the four C3 symmetries which must all be detected in order to fully
    describe tetrahedral symmetry. If all four are found, the ret vector will contain these as its only four entries, while it will be empty if some of the
    C3 symmetries are not found. The missing symmetry axis detection is implemented as part of this function as well.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret The vector .
    \param[in] axErr The error tolerance on angle matching.
    \param[in] verbose How loud the announcments should be?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
 */
void ProSHADE_internal_symmetry::findTetra4C3s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C3PossibilitiesHlp;
    std::vector< std::vector< proshade_unsign > > C3Possibilities;
    bool groupMatched;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage ( verbose, 2, "Starting detection of four C3 axes.", messageShift );
    
    //================================================ For all symmetries in the C symmetries list
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
    {
        //============================================ Search only using C3s
        if ( CSymList->at(cIt)[0] != 3.0 || CSymList->at(cIt)[0] < minPeakHeight ) { continue; }
        
        //============================================ If this is the first C3, then just save it to the first group of the temporary holder
        if ( C3Possibilities.size() == 0 ) { ProSHADE_internal_misc::addToUnsignVector ( &C3PossibilitiesHlp, cIt ); ProSHADE_internal_misc::addToUnsignVectorVector ( &C3Possibilities, C3PossibilitiesHlp ); continue; }
        
        //============================================ If second or more C3, check if it has the correct angle to all other already found C3s for each group
        groupMatched                                  = false;
        for ( proshade_unsign gIt = 0; gIt < static_cast<proshade_unsign> ( C3Possibilities.size() ); gIt++ )
        {
            if ( ProSHADE_internal_symmetry::testGroupAgainstSymmetry ( CSymList, &C3Possibilities.at(gIt), CSymList->at(cIt), axErr, 1.0/3.0, true, cIt ) ) { ProSHADE_internal_misc::addToUnsignVector ( &C3Possibilities.at(gIt), cIt ); groupMatched = true; break; }
        }
        
        //============================================ If no group matched, create a new group
        if ( !groupMatched ) { C3PossibilitiesHlp.clear(); ProSHADE_internal_misc::addToUnsignVector ( &C3PossibilitiesHlp, cIt ); ProSHADE_internal_misc::addToUnsignVectorVector ( &C3Possibilities, C3PossibilitiesHlp ); continue; }
    }
    
    //================================================ Test for missing symmetry axes, if need be
    ProSHADE_internal_symmetry::findMissingAxes       ( &C3Possibilities, CSymList, 4, axErr, 1.0/3.0, 3, dataObj, minPeakHeight );
    
    //================================================ Any group has 4 entries? If more such groups, take the one with highest average height.
    proshade_double maxHeight = 0.0; proshade_unsign maxGrp = 0;
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( C3Possibilities.size() ); iter++ ) { if ( C3Possibilities.at(iter).size() == 4 ) { if ( ( ( CSymList->at(C3Possibilities.at(iter).at(0))[5] + CSymList->at(C3Possibilities.at(iter).at(1))[5] + CSymList->at(C3Possibilities.at(iter).at(2))[5] + CSymList->at(C3Possibilities.at(iter).at(3))[5] ) / 4.0 ) > maxHeight ) { maxHeight = ( ( CSymList->at(C3Possibilities.at(iter).at(0))[5] + CSymList->at(C3Possibilities.at(iter).at(1))[5] + CSymList->at(C3Possibilities.at(iter).at(2))[5] + CSymList->at(C3Possibilities.at(iter).at(3))[5] ) / 4.0 ); maxGrp = iter; } } }
    
    if ( C3Possibilities.at(maxGrp).size() == 4 )
    {
        //============================================ Success! Save and exit
        for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( C3Possibilities.at(maxGrp).size() ); it++ ) { ProSHADE_internal_misc::addToDblPtrVector ( ret, CSymList->at(C3Possibilities.at(maxGrp).at(it)) ); }
        
        //============================================ Report progress
        ProSHADE_internal_messages::printProgressMessage ( verbose, 3, "Detection of four C3 axes successfull.", messageShift );
        
        //============================================ Done
        return ;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function tests whether a  symmetry has particular angle to all members of a group.
 
    This utility function tests if a sinlge symmetry axis has a given angle to all member of a particular symmetry group as given by the
    vector of indices and a vector of all symmetries. If the improve parameter is true, that it will also check for the tested axis for
    being parallel to any of the group axes while having higher average peak height - and in such cases, the function will replace the existing
    axis with the tested axis index as given in the pos argument. This utility is useful when searching for all axes of polyhedral symmetry groups.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] grp A vector of indices (relating to CSymList) of the group members.
    \param[in] sym A double pointer to array containing the symmetry to be tested against the group.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] angle The angle that each group member is required to have against the symmetry.
    \param[in] improve Boolead value stating whether an axis with higher average height should be used instead of equal axis with lower average height, if such axis is found.
    \param[in] pos This is the CSymList index of the axis tested against the group. It will be used if improve = true to change the grp entry which is identical, but has lower height.
    \param[out] X Boolean value speciying whether all group members have the angle to the symmetry or not.
 */
bool ProSHADE_internal_symmetry::testGroupAgainstSymmetry ( std::vector< proshade_double* >* CSymList, std::vector< proshade_unsign >* grp, proshade_double* sym, proshade_double axErr, proshade_double angle, bool improve, proshade_unsign pos )
{
    //================================================ Initialise variables
    bool allAnglesMet                                 = true;
    proshade_double dotProduct;
    
    //================================================ Improve if required
    if ( improve )
    {
        for ( proshade_unsign mIt = 0; mIt < static_cast<proshade_unsign> ( grp->size() ); mIt++ )
        {
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &CSymList->at(grp->at(mIt))[1], &CSymList->at(grp->at(mIt))[2], &CSymList->at(grp->at(mIt))[3], &sym[1], &sym[2], &sym[3] );
            
            if ( ( ( 1.0 > ( dotProduct - axErr ) ) && ( 1.0 < ( dotProduct + axErr ) ) ) || ( ( -1.0 > ( dotProduct - axErr ) ) && ( -1.0 < ( dotProduct + axErr ) ) ) )
            {
                if ( sym[5] > CSymList->at(grp->at(mIt))[5] )
                {
                    grp->at(mIt)                      = pos;
                }
                else
                {
                    allAnglesMet                      = false;
                    return                            ( allAnglesMet );
                }
            }
        }
    }
    
    //================================================ For all group members
    for ( proshade_unsign mIt = 0; mIt < static_cast<proshade_unsign> ( grp->size() ); mIt++ )
    {
        dotProduct                                    = ProSHADE_internal_maths::computeDotProduct ( &CSymList->at(grp->at(mIt))[1], &CSymList->at(grp->at(mIt))[2], &CSymList->at(grp->at(mIt))[3], &sym[1], &sym[2], &sym[3] );
        
        if ( ( angle > ( std::abs ( dotProduct ) - axErr ) ) &&
             ( angle < ( std::abs ( dotProduct ) + axErr ) ) )
        {
            //======================================== Matching group memner - try next one
        }
        else
        {
            //======================================== Group member not matched - try next group
            allAnglesMet                              = false;
            break;
        }
    }
    
    //================================================ Done
    return                                           ( allAnglesMet );
    
}

/*! \brief This function tries to find an axis which would complete a particular group of axes for polyhedral symmetry detection.
 
    This function assumes that there is a set of already detected axes and that for a polyhedral symmetry, another axis with known fold and angle to some of the already detected axis needs
    to be found. It uses algebraic solution to try to find such an axis (or a given number of them) and also tests for these newly detected axes being unique and having at least minPeakHeight
    average peak height. If such axes are found, they are added to the CSymList vector and their indices are also added to the possibilities vector.
 
    \param[in] possibilities A vector of vectors of indices to the cyclic symmetries list with all the already determined axes.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] requiredNoAxes Number of axes required for positive result.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] angle The angle that each group member is required to have against the symmetry.
    \param[in] fold The fold of the searched for axis.
    \param[in] dataObj The full data holding object pointer - this is to get access to self-rotation function values.
    \param[in] minPeakHeight The minimum new axis average peak height in order for the axis to be added.
    \param[out] atLeastOne Boolean value speciying whether at least the minimum required number of axes was found.
 */
bool ProSHADE_internal_symmetry::findMissingAxes ( std::vector< std::vector< proshade_unsign > >* possibilities, std::vector< proshade_double* >* CSymList, proshade_unsign requiredNoAxes, proshade_double axErr, proshade_double angle, proshade_unsign fold, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_double* > hlpVec;
    bool atLeastOne                                   = false;
    
    //================================================ Proceed only if need be
    for ( proshade_unsign gIt = 0; gIt < static_cast<proshade_unsign> ( possibilities->size() ); gIt++ )
    {
        if ( static_cast<proshade_unsign> ( possibilities->at(gIt).size() ) == requiredNoAxes ) { atLeastOne = true; return ( atLeastOne ); }
    }
    
    //================================================ For each possible group
    for ( proshade_unsign gIt = 0; gIt < static_cast<proshade_unsign> ( possibilities->size() ); gIt++ )
    {
        //============================================ This will not work for less than two axes in group
        if ( possibilities->at(gIt).size() < 2 ) { continue; }
        
        //============================================ Prepare iteration
        hlpVec.clear                                  ( );
        
        //============================================ Search for missing axes
        ProSHADE_internal_symmetry::searchMissingSymmetrySpace ( dataObj, CSymList, &possibilities->at(gIt), &hlpVec, axErr, angle, fold, minPeakHeight );
        
        //============================================ Add missing axes
        if ( hlpVec.size() > 0 )
        {
            //======================================== Start adding by highest first
            std::sort                                 ( hlpVec.begin(), hlpVec.end(), ProSHADE_internal_misc::sortSymHlpInv );
            
            //======================================== For each missing axis
            for ( proshade_unsign axIt = 0; axIt < static_cast<proshade_unsign> ( hlpVec.size() ); axIt++ )
            {
                if ( ProSHADE_internal_symmetry::testGroupAgainstSymmetry ( CSymList, &possibilities->at(gIt), hlpVec.at(axIt), axErr, angle, false ) )
                {
                    //================================ Check for uniqueness
                    if ( ProSHADE_internal_maths::isAxisUnique ( CSymList, hlpVec.at(axIt), axErr ) )
                    {
                        //============================ Add
                        ProSHADE_internal_misc::addToDblPtrVector ( CSymList, hlpVec.at(axIt) );
                        ProSHADE_internal_misc::addToUnsignVector ( &possibilities->at(gIt), static_cast<proshade_unsign> ( CSymList->size()-1 ) );
                    }
                }
            }
        }
        
        if ( possibilities->at(gIt).size() == requiredNoAxes ) { atLeastOne = true; }
    }
    
    //================================================ Done
    return                                            ( atLeastOne );
    
}

/*! \brief This function compares two arrays of two based on the first number.
 
 \param[in] a The first array to compare.
 \param[in] b The second array to compare.
 \param[out] X Boolean whether the first is smaller than the second.
 */
bool ProSHADE_internal_symmetry::sortArrVecHlp ( const proshade_double* a, const proshade_double* b )
{
    //================================================ Compare
    return                                            ( a[0] < b[0] );
    
}

/*! \brief This function searches for the highest peaks average that would produce the required axis and fold.
 
    This function starts by finding all self-rotation map points with corresponding axis and recording the angle and map heights of these points. It then
    sorts these and searches for a combination of fold points separated by the 2pi/fold distance with the highest average map height. In this way, the highest
    average symmetry height is determined for any axis. This does not, however, check if such symmetry does indeed exist!
 
    \param[in] xVal The x-axis element of the axis to have the height detected.
    \param[in] yVal The y-axis element of the axis to have the height detected.
    \param[in] zVal The z-axis element of the axis to have the height detected.
    \param[in] dataObj The full data holding object pointer - this is to get access to self-rotation function values.
    \param[in] fold The fold of the searched for axis.
    \param[in] axErr The error tolerance on angle matching.
    \param[out] X The highest height value found for the axis with the given fold.
 */
proshade_double ProSHADE_internal_symmetry::missingAxisHeight ( proshade_double xVal, proshade_double yVal, proshade_double zVal, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_unsign fold, proshade_double axErr )
{
    //================================================ Initialise variables
    proshade_double ret                               = 0.0;
    proshade_double curSum                            = 0.0;
    proshade_double maxVal                            = 0.0;
    proshade_double angStep                           = std::acos ( 1.0 - axErr ) / 2;
    std::vector< proshade_double* > angVec;
    
    //================================================ Find map points conforming to the axis
    angVec                                            = ProSHADE_internal_symmetry::findMissingAxisPoints ( xVal, yVal, zVal, dataObj, axErr );
    
    //================================================ Sort points by angle
    std::sort                                         ( angVec.begin(), angVec.end(), ProSHADE_internal_symmetry::sortArrVecHlp );
    
    //================================================ Find the best X peaks with correct distances
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( std::floor ( ( 2.0 * M_PI / angStep ) / static_cast< proshade_double > ( fold ) ) ); iter++ )
    {
        //============================================ Initialise new ang group iteration
        curSum                                        = 0.0;
        
        //============================================ For each of the fold times
        for ( proshade_unsign angCmb = 0; angCmb < static_cast<proshade_unsign> ( fold ); angCmb++ )
        {
            //======================================== Initialise
            maxVal                                    = 0.0;
            
            //======================================== Search
            for ( proshade_unsign angIt = 0; angIt < static_cast<proshade_unsign> ( angVec.size() ); angIt++ )
            {
                if ( angVec.at(angIt)[0] < ( ( static_cast< proshade_double > ( iter ) * angStep ) +
                                             ( ( 2.0 * M_PI / static_cast< proshade_double > ( fold ) ) * static_cast< proshade_double > ( angCmb ) ) ) ) { continue; }
                if ( angVec.at(angIt)[0] > ( ( ( static_cast< proshade_double > ( iter ) + 1.0 ) * angStep ) +
                                             ( ( 2.0 * M_PI / static_cast< proshade_double > ( fold ) ) * static_cast< proshade_double > ( angCmb ) ) ) ) { break; }

                if ( angVec.at(angIt)[1] > maxVal ) { maxVal = angVec.at(angIt)[1]; }
            }
            curSum                                   += maxVal;
        }
        curSum                                      /= static_cast<proshade_double> ( fold );
        if ( ret < curSum ) { ret = curSum; }
    }
    
    //================================================ Release memory
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( angVec.size() ); iter++ ) { delete[] angVec.at(iter); }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function searches for all the self-rotation map points conforming to the axis, returning their angles and heights.
 
    This helper function searches the self-rotation map point by point for all points which represent the same rotation axis as required
    by the input parameters. For all such points, it records the angle they represent and the map height associated with them. Finally, it
    returns a vector of all detected points.
 
    \param[in] xVal The x-axis element of the axis to have the height detected.
    \param[in] yVal The y-axis element of the axis to have the height detected.
    \param[in] zVal The z-axis element of the axis to have the height detected.
    \param[in] dataObj The full data holding object pointer - this is to get access to self-rotation function values.
    \param[in] axErr The error tolerance on angle matching.
    \param[out] angVec Vector containing all map points which conform to the required axis along with their heights.
 */
std::vector < proshade_double* > ProSHADE_internal_symmetry::findMissingAxisPoints ( proshade_double xVal, proshade_double yVal, proshade_double zVal, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_double axErr )
{
    //================================================ Initialise variables
    proshade_double euA, euB, euG, xPk, yPk, zPk, anglPk;
    proshade_double* rotMat                           = new proshade_double [9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    proshade_unsign arrIndex;
    std::vector< proshade_double* > angVec;
    
    //================================================ Search the self-rotation map
    for ( proshade_unsign xIt = 0; xIt < ( dataObj->getMaxBand() * 2 ); xIt++ )
    {
        for ( proshade_unsign yIt = 0; yIt < ( dataObj->getMaxBand() * 2 ); yIt++ )
        {
            for ( proshade_unsign zIt = 0; zIt < ( dataObj->getMaxBand() * 2 ); zIt++ )
            {
                //==================================== Get height and check against threshold
                arrIndex                              = zIt  + ( dataObj->getMaxBand() * 2 ) * ( yIt  + ( dataObj->getMaxBand() * 2 ) * xIt );
                
                //==================================== Get angle-axis values
                ProSHADE_internal_maths::getEulerZXZFromSOFTPosition ( static_cast< proshade_signed > ( dataObj->getMaxBand() ), static_cast< proshade_signed > ( xIt ),
                                                                       static_cast< proshade_signed > ( yIt ), static_cast< proshade_signed > ( zIt ),
                                                                       &euA, &euB, &euG );
                ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( euA, euB, euG, rotMat );
                ProSHADE_internal_maths::getAxisAngleFromRotationMatrix ( rotMat, &xPk, &yPk, &zPk, &anglPk );
                
                //==================================== Set largest axis element to positive
                const FloatingPoint< proshade_double > lhs1 ( std::max ( std::abs ( xPk ), std::max( std::abs ( yPk ), std::abs ( zPk ) ) ) );
                const FloatingPoint< proshade_double > rhs1 ( std::abs ( xPk ));
                const FloatingPoint< proshade_double > rhs2 ( std::abs ( yPk ) );
                const FloatingPoint< proshade_double > rhs3 ( std::abs ( zPk ) );
                if ( ( lhs1.AlmostEquals ( rhs1 ) && ( xPk < 0.0 ) ) ||
                     ( lhs1.AlmostEquals ( rhs2 ) && ( yPk < 0.0 ) ) ||
                     ( lhs1.AlmostEquals ( rhs3 ) && ( zPk < 0.0 ) ) )
                {
                    xPk                              *= -1.0;
                    yPk                              *= -1.0;
                    zPk                              *= -1.0;
                    anglPk                           *= -1.0;
                }
                
                //==================================== Does the peak match the required axis?
                if ( ProSHADE_internal_maths::vectorOrientationSimilarity ( xPk, yPk, zPk, xVal, yVal, zVal, axErr ) )
                {
                    //================================ Matching map point - save it
                    proshade_double* hlpArr           = new proshade_double [2];
                    ProSHADE_internal_misc::checkMemoryAllocation ( hlpArr, __FILE__, __LINE__, __func__ );
                    hlpArr[0]                         = anglPk + M_PI;
                    hlpArr[1]                         = pow( dataObj->getInvSO3Coeffs()[arrIndex][0], 2.0 ) +
                                                        pow( dataObj->getInvSO3Coeffs()[arrIndex][1], 2.0 );
                    ProSHADE_internal_misc::addToDblPtrVector ( &angVec, hlpArr );
                }
            }
        }
    }
    
    //================================================ Release memory
    delete[] rotMat;
    
    //================================================ Done
    return                                            ( angVec );
    
}

/*! \brief This function saves the recovered information about missing axis into a full symmetry, making sure no duplicates are created.
 
    This function takes the information about the missing symmetry and proceeds to create a full symmetry description out of it. It then checks whether
    the vector already contains similar symmetry, either replacing the old or ignoring the new symmetry based on which has hiher height. If the symmetry
    does not match anything in the vector, it will be copied as a new vector entry.
 
    \param[in] axVec Vector containing all already detected missing axes.
    \param[in] axX The x-axis element of the missing axis.
    \param[in] axY The y-axis element of the missing axis.
    \param[in] axZ The z-axis element of the missing axis.
    \param[in] height The average map height for this new axis.
    \param[in] fold The fold of the searched for axis.
    \param[in] axErr The error tolerance on angle matching.
 */
void ProSHADE_internal_symmetry::saveMissingAxisNewOnly ( std::vector< proshade_double* >* axVec, proshade_double axX, proshade_double axY, proshade_double axZ, proshade_double height, proshade_unsign fold, proshade_double axErr )
{
    //================================================  Create symmetry array from the inputs
    proshade_double* hlpSym                           = new proshade_double [6];
    ProSHADE_internal_misc::checkMemoryAllocation     ( hlpSym, __FILE__, __LINE__, __func__ );
    
    //================================================ Fill it in
    hlpSym[0]                                         = static_cast<proshade_double> ( fold );
    hlpSym[1]                                         = axX;
    hlpSym[2]                                         = axY;
    hlpSym[3]                                         = axZ;
    hlpSym[4]                                         = ( 2.0 * M_PI ) / static_cast<proshade_double> ( fold );
    hlpSym[5]                                         = height;
    
    //================================================ Check if similar symmetry does not exist already
    for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( axVec->size() ); symIt++ )
    {
        //============================================ Minor speed-up => only test for same folds
        const FloatingPoint< proshade_double > lhs1 ( axVec->at(symIt)[0] ), rhs1 ( hlpSym[0] );
        if ( lhs1.AlmostEquals ( rhs1 ) )
        {
            if ( ProSHADE_internal_maths::vectorOrientationSimilarity ( axVec->at(symIt)[1],
                                                                        axVec->at(symIt)[2],
                                                                        axVec->at(symIt)[3],
                                                                        hlpSym[1],
                                                                        hlpSym[2],
                                                                        hlpSym[3],
                                                                        axErr ) )
            {
                //==================================== Almost identical entry
                if ( axVec->at(symIt)[5] < hlpSym[5] )
                {
                    //================================ If higher, save
                    delete[] axVec->at(symIt);
                    axVec->at(symIt)                  = hlpSym;
                    return ;
                }
                else
                {
                    //================================ or just terminate if better is already saved
                    delete[] hlpSym;
                    return ;
                }
            }
        }
    }
    
    //================================================ Not matched to anything
    ProSHADE_internal_misc::addToDblPtrVector         ( axVec, hlpSym );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function tests feasible axes against the missing axis criteria, returning a set of matching axes.
 
    This function does the real missing axis searching. It starts by taking all supplied axes and algebraically computing the vector which has the required angle to two of the supplied axes. This computed axis is
    then tested against the group for being unique and having an average height at least as high as required. If such axis is found, it is added to the axes list.
 
    \param[in] dataObj The full data holding object pointer - this is to get access to self-rotation function values.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] grp A vector of indices (relating to CSymList) of the group members.
    \param[in] hlpVec A vector which will hold the detected, but not verified axes to be returned to the caller function.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] angle The angle that each group member is required to have against the symmetry.
    \param[in] fold The fold of the searched for axis.
    \param[in] minPeakHeight The minimum new axis average peak height in order for the axis to be added.
 */
void ProSHADE_internal_symmetry::searchMissingSymmetrySpace ( ProSHADE_internal_data::ProSHADE_data* dataObj, std::vector< proshade_double* >* CSymList, std::vector< proshade_unsign >* grp, std::vector< proshade_double* >* hlpVec, proshade_double axErr, proshade_double angle, proshade_unsign fold, proshade_double minPeakHeight )
{
    //================================================ Sanity check
    if ( grp->size() < 2 ) { return; }
    
    //================================================ Initialise variables
    proshade_double axHeight                          = 0.0;
    proshade_double* symHlp                           = new proshade_double[7];
    ProSHADE_internal_misc::checkMemoryAllocation     ( symHlp, __FILE__, __LINE__, __func__ );
    
    //================================================ For each axis pair in the group, find the possible solutions
    for ( proshade_unsign fAx = 0; fAx < static_cast<proshade_unsign> ( grp->size() ); fAx++ )
    {
        for ( proshade_unsign sAx = 1; sAx < static_cast<proshade_unsign> ( grp->size() ); sAx++ )
        {
            //======================================== Only unique pairs
            if ( fAx >= sAx ) { continue; }
            
            //======================================== Find possible axis having the required angle to this pair ( solution 1 )
            std::vector< proshade_double > solVec     = ProSHADE_internal_maths::findVectorFromTwoVAndTwoD ( CSymList->at(grp->at(fAx))[1],
                                                                                                             CSymList->at(grp->at(fAx))[2],
                                                                                                             CSymList->at(grp->at(fAx))[3],
                                                                                                             CSymList->at(grp->at(sAx))[1],
                                                                                                             CSymList->at(grp->at(sAx))[2],
                                                                                                             CSymList->at(grp->at(sAx))[3], angle, angle );
            
            //======================================== Set largest axis element to positive
            const FloatingPoint< proshade_double > lhs1 ( std::max ( std::abs ( solVec.at(0) ), std::max( std::abs ( solVec.at(1) ), std::abs ( solVec.at(2) ) ) ) );
            const FloatingPoint< proshade_double > rhs1 ( std::abs ( solVec.at(0) ) );
            const FloatingPoint< proshade_double > rhs2 ( std::abs ( solVec.at(1) ) );
            const FloatingPoint< proshade_double > rhs3 ( std::abs ( solVec.at(2) ) );
            if ( ( lhs1.AlmostEquals ( rhs1 ) && ( solVec.at(0) < 0.0 ) ) ||
                 ( lhs1.AlmostEquals ( rhs2 ) && ( solVec.at(1) < 0.0 ) ) ||
                 ( lhs1.AlmostEquals ( rhs3 ) && ( solVec.at(2) < 0.0 ) ) )
            {
                solVec.at(0)                         *= -1.0;
                solVec.at(1)                         *= -1.0;
                solVec.at(2)                         *= -1.0;
            }
            
            //======================================== Does the solution fit the whole group?
            symHlp[1] = solVec.at(0); symHlp[2] = solVec.at(1); symHlp[3] = solVec.at(2); symHlp[6] = -std::numeric_limits < proshade_double >::infinity();
            if ( ProSHADE_internal_symmetry::testGroupAgainstSymmetry ( CSymList, grp, symHlp, axErr, angle, true ) )
            {
                //==================================== Find the height for the axis
                axHeight                              = ProSHADE_internal_symmetry::missingAxisHeight ( solVec.at(0), solVec.at(1), solVec.at(2), dataObj, fold, axErr );
                
                //================================ Save max height result
                if ( axHeight >= minPeakHeight ) { ProSHADE_internal_symmetry::saveMissingAxisNewOnly ( hlpVec, solVec.at(0), solVec.at(1), solVec.at(2), axHeight, fold, axErr ); }
            }
            
            //======================================== Find possible axis having the required angle to this pair ( solution 2 )
            solVec                                    = ProSHADE_internal_maths::findVectorFromTwoVAndTwoD ( CSymList->at(grp->at(fAx))[1],
                                                                                                             CSymList->at(grp->at(fAx))[2],
                                                                                                             CSymList->at(grp->at(fAx))[3],
                                                                                                             CSymList->at(grp->at(sAx))[1],
                                                                                                             CSymList->at(grp->at(sAx))[2],
                                                                                                             CSymList->at(grp->at(sAx))[3], -angle, -angle );
            
            //======================================== Set largest axis element to positive
            const FloatingPoint< proshade_double > lhs2 ( std::max ( std::abs ( solVec.at(0) ), std::max( std::abs ( solVec.at(1) ), std::abs ( solVec.at(2) ) ) ) );
            const FloatingPoint< proshade_double > rhs4 ( std::abs ( solVec.at(0) ) );
            const FloatingPoint< proshade_double > rhs5 ( std::abs ( solVec.at(1) ) );
            const FloatingPoint< proshade_double > rhs6 ( std::abs ( solVec.at(2) ) );
            if ( ( lhs2.AlmostEquals ( rhs4 ) && ( solVec.at(0) < 0.0 ) ) ||
                 ( lhs2.AlmostEquals ( rhs5 ) && ( solVec.at(1) < 0.0 ) ) ||
                 ( lhs2.AlmostEquals ( rhs6 ) && ( solVec.at(2) < 0.0 ) ) )
            {
                solVec.at(0)                         *= -1.0;
                solVec.at(1)                         *= -1.0;
                solVec.at(2)                         *= -1.0;
            }
            
            //======================================== Does the solution fit the whole group?
            symHlp[1] = solVec.at(0); symHlp[2] = solVec.at(1); symHlp[3] = solVec.at(2); symHlp[6] = -std::numeric_limits < proshade_double >::infinity();
            if ( ProSHADE_internal_symmetry::testGroupAgainstSymmetry ( CSymList, grp, symHlp, axErr, angle, true ) )
            {
                //==================================== Find the height for the axis
                axHeight                              = ProSHADE_internal_symmetry::missingAxisHeight ( solVec.at(0), solVec.at(1), solVec.at(2), dataObj, fold, axErr );
                
                //================================ Save max height result
                if ( axHeight >= minPeakHeight ) { ProSHADE_internal_symmetry::saveMissingAxisNewOnly ( hlpVec, solVec.at(0), solVec.at(1), solVec.at(2), axHeight, fold, axErr ); }
            }
        }
    }
    
    //================================================ Release memory
    delete[] symHlp;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes the list of C symmetries and finds the 3 C2 symmetries with correct angles required for full tetrahedral symmetry.
 
    This is a specific helper function for detecting three C2 symmetries perpendicular to each other hand having a specific angle ( acos(0.5) ) to one of the already
    detected C3 symmetries of the sought after tetrahedral symmetry. It firstly finds all C2s and tests these for having the acos(0.5) angle to the already found C3s.
    From this list of passing C2s, it then tries to find three mutually perpendicular axes, including searching for missing axes. If no such axes are found, the ret
    array will still have 4 entries, while if they are found, the ret array will have these added to the total of 7 entries.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret The vector .
    \param[in] axErr The error tolerance on angle matching.
    \param[in] verobse How loud the announcments should be?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
 */
void ProSHADE_internal_symmetry::findTetra3C2s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C3s, prospectiveC2s, C2PossibilitiesHlp;
    std::vector< std::vector< proshade_unsign > > C2Possibilities;
    proshade_double dotProd;
    bool groupMatched;
    for ( proshade_unsign iter = 0; iter < 4; iter++ ) { ProSHADE_internal_misc::addToUnsignVector ( &C3s, iter ); }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of three C2 axes.", messageShift );
    
    //================================================ For each C3
    for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( ret->size() ); rIt++ )
    {
        //============================================ For each C2, check it has angle ( acos(0.5) ) to the tested C3
        for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
        {
            //======================================== Search only using C2s
            const FloatingPoint< proshade_double > lhs999 ( CSymList->at(cIt)[5] ), rhs999 ( static_cast< proshade_double > ( -999.9 ) );
            if ( CSymList->at(cIt)[0] != 2.0 || ( ( CSymList->at(cIt)[5] < minPeakHeight ) && !( lhs999.AlmostEquals( rhs999 ) ) ) ) { continue; }
            
            //======================================== Check the C2 axis to the C3 ( acos ( 0.5 ) )
            dotProd = ProSHADE_internal_maths::computeDotProduct ( &ret->at(rIt)[1], &ret->at(rIt)[2], &ret->at(rIt)[3],
                                                                   &CSymList->at(cIt)[1], &CSymList->at(cIt)[2], &CSymList->at(cIt)[3] );
            
            if ( ( std::abs ( dotProd ) > ( 0.5 - axErr ) ) && ( std::abs ( dotProd ) < ( 0.5 + axErr ) ) ) { ProSHADE_internal_misc::addToUnsignVector ( &prospectiveC2s, cIt ); }
        }
    }
        
    //================================================ Group the prospective C2s
    C2Possibilities.clear(); C2PossibilitiesHlp.clear();
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( prospectiveC2s.size() ); cIt++ )
    {
        //============================================ If second or more C2, check if it can be placed in any group with being perpendicular to all its members
        groupMatched                                  = false;
        for ( proshade_unsign gIt = 0; gIt < static_cast<proshade_unsign> ( C2Possibilities.size() ); gIt++ )
        {
            if ( ProSHADE_internal_symmetry::testGroupAgainstSymmetry ( CSymList, &C2Possibilities.at(gIt), CSymList->at(prospectiveC2s.at(cIt)), axErr, 0.0, true, prospectiveC2s.at(cIt) ) ) { ProSHADE_internal_misc::addToUnsignVector ( &C2Possibilities.at(gIt), prospectiveC2s.at(cIt) ); groupMatched = true; break; }
        }
        
        //============================================ If no group matched, create a new group
        if ( !groupMatched ) { C2PossibilitiesHlp.clear(); ProSHADE_internal_misc::addToUnsignVector ( &C2PossibilitiesHlp, prospectiveC2s.at(cIt) ); ProSHADE_internal_misc::addToUnsignVectorVector ( &C2Possibilities, C2PossibilitiesHlp ); continue; }
    }
    
    //================================================ Find the best group or return empty
    while ( C2Possibilities.size() != 0 )
    {
        //============================================ Test for missing symmetry axes, if need be
        ProSHADE_internal_symmetry::findMissingAxes   ( &C2Possibilities, CSymList, 3, axErr, 0.0, 2, dataObj, minPeakHeight );
        
        //============================================ Found 3 C2s?
        if ( C2Possibilities.at(0).size() == 3 )
        {
            //======================================== Success! Save and exit
            for ( proshade_unsign it = 0; it < 3; it++ ) { ProSHADE_internal_misc::addToDblPtrVector ( ret, CSymList->at(C2Possibilities.at(0).at(it)) ); }
            
            //======================================== Report progress
            ProSHADE_internal_messages::printProgressMessage ( verbose, 3, "Detection of three C2 axes successfull.", messageShift );
            
            //======================================== Done
            return ;
        }
        else { C2Possibilities.erase ( C2Possibilities.begin() ); }
    }
  
    //================================================ Done
    return ;
    
}

/*! \brief This function compares two groups of axes for a single pair having the required angle.
 
    This simple helper function takes two sets of symmetry axes and two vectors of indices, each relating to one of the two sets. It then proceeds to
    check each of the indexed axes in each set against all the indexed axes in the other set, searching for a particular angle. If this angle is found
    for at least one pair, true is returned, while otherwise false is returned.
 
    \param[in] GrList1 A vector containing the symmetries for the group 1.
    \param[in] grp1 The indices respective to GrList1 which form group 1.
    \param[in] GrList2 A vector containing the symmetries for the group 2.
    \param[in] grp2 The indices respective to GrList1 which form group 2.
    \param[in] angle The angle which needs to be found between any pair of axes in group 1 and 2.
    \param[in] axErr The error tolerance on angle matching.
    \param[out] ret True if succeeded, false otherwise.
 */
bool ProSHADE_internal_symmetry::testGroupAgainstGroup ( std::vector< proshade_double* >* GrList1, std::vector< proshade_unsign >* grp1, std::vector< proshade_double* >* GrList2, std::vector< proshade_unsign >* grp2, proshade_double angle, proshade_double axErr )
{
    //================================================ Initialise variables
    bool ret                                          = false;
    proshade_double dotProduct;
    
    //================================================ For all pairs of axes
    for ( proshade_unsign g1It = 0; g1It < static_cast<proshade_unsign> ( grp1->size() ); g1It++ )
    {
        for ( proshade_unsign g2It = 0; g2It < static_cast<proshade_unsign> ( grp2->size() ); g2It++ )
        {
            //======================================== Find the angle
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &GrList1->at(grp1->at(g1It))[1],
                                                                                                     &GrList1->at(grp1->at(g1It))[2],
                                                                                                     &GrList1->at(grp1->at(g1It))[3],
                                                                                                     &GrList2->at(grp2->at(g2It))[1],
                                                                                                     &GrList2->at(grp2->at(g2It))[2],
                                                                                                     &GrList2->at(grp2->at(g2It))[3] );
            
            //======================================== Check the angle
            if ( ( angle > ( dotProduct - axErr ) ) && ( angle < ( dotProduct + axErr ) ) )
            {
                ret                                   = true;
                return                                ( ret );
            }
        }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function obtains a list of all O symmetry axes from the already computed C symmetries list.
 
    This function starts by checking if there is a pair of C3 and C4 symmetries with the octahedron dihedral angle ( acos ( 1/sqrt(3) ) ). If so, it will
    then assume existence of octahedral symmetry and it will search for three C4 axes, four C3 axes and six C2 axes with the correct angle to each other
    and within the group. If all required axes are detected, it will return a list of 13 axes, otherwise it will return empty or shorter list. Automated
    missing symmetry axis detection is also included.
 
    \param[in] settings A pointer to settings class containing all the information required for symmetry detection.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
 */
std::vector< proshade_double* > ProSHADE_internal_data::ProSHADE_data::getOctahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList )
{
    //================================================ Initialise variables
    std::vector< proshade_double* > ret;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting O symmetry detection.", settings->messageShift );
    
    //================================================ Are the basic requirements for tetrahedral symmetry met?
    if ( ProSHADE_internal_symmetry::detectOctahedralSymmetry ( CSymList, settings->axisErrTolerance, settings->minSymPeak ) )
    {
        //============================================ Search for all the symmetry axes
        ProSHADE_internal_symmetry::findOcta3C4s ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->messageShift, settings->minSymPeak );
        if ( ret.size() != 3 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        
        ProSHADE_internal_symmetry::findOcta4C3s ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->messageShift, settings->minSymPeak );
        if ( ret.size() != 7 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        
        ProSHADE_internal_symmetry::findOcta6C2s ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->messageShift, settings->minSymPeak );
        if ( ret.size() != 13 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        else
        {
            for ( proshade_unsign csIt = 0; csIt < static_cast<proshade_unsign> ( CSymList->size() ); csIt++ )
            {
                for ( proshade_unsign retIt = 0; retIt < static_cast<proshade_unsign> ( ret.size() ); retIt++ )
                {
                    //======================================== Sort ret by fold
                    std::sort                                 ( ret.begin(), ret.end(), ProSHADE_internal_misc::sortSymInvFoldHlp );
                    
                    //======================================== Save indices
                    const FloatingPoint< proshade_double > lhs1 ( CSymList->at(csIt)[0] ), rhs1 ( ret.at(retIt)[0] );
                    const FloatingPoint< proshade_double > lhs2 ( CSymList->at(csIt)[1] ), rhs2 ( ret.at(retIt)[1] );
                    const FloatingPoint< proshade_double > lhs3 ( CSymList->at(csIt)[2] ), rhs3 ( ret.at(retIt)[2] );
                    const FloatingPoint< proshade_double > lhs4 ( CSymList->at(csIt)[3] ), rhs4 ( ret.at(retIt)[3] );
                    const FloatingPoint< proshade_double > lhs5 ( CSymList->at(csIt)[4] ), rhs5 ( ret.at(retIt)[4] );
                    const FloatingPoint< proshade_double > lhs6 ( CSymList->at(csIt)[5] ), rhs6 ( ret.at(retIt)[5] );
                    if ( lhs1.AlmostEquals ( rhs1 ) &&
                         lhs2.AlmostEquals ( rhs2 ) &&
                         lhs3.AlmostEquals ( rhs3 ) &&
                         lhs4.AlmostEquals ( rhs4 ) &&
                         lhs5.AlmostEquals ( rhs5 ) &&
                         lhs6.AlmostEquals ( rhs6 ) )
                    {
                        ProSHADE_internal_misc::addToUnsignVector ( &settings->allDetectedOAxes, csIt );
                    }
                }
            }
        }
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "O symmetry detection complete.", settings->messageShift );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function takes the list of C symmetries and decides whether basic requirements for octahhedral symmetry are there.
 
    This function first finds all the C4 symmetries in the C symmetries list and then it checks each present C4 against all C3 symmetries for having
    the angle between the pair equal to the dihedral angle of an octahedron ( acos(1/sqrt(3)) ). If a single such pair is detected, this is likely an
    octahedral symmetry and all other axes need to be located. Otherwise, false is returned.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
    \param[out] X Boolean value telling whether there are C4 and C3 symmetries with octahedral dihhedral angle.
 */
bool ProSHADE_internal_symmetry::detectOctahedralSymmetry ( std::vector< proshade_double* >* CSymList, proshade_double axErr, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C4List;
    proshade_double dotProduct;
    
    //================================================ Find all C4 symmetries
    for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ )
    {
        const FloatingPoint< proshade_double > lhs1 ( CSymList->at(cSym)[0] ), rhs1 ( 4.0 );
        if ( lhs1.AlmostEquals ( rhs1 ) && CSymList->at(cSym)[5] >= minPeakHeight ) { ProSHADE_internal_misc::addToUnsignVector ( &C4List, cSym ); }
    }
    
    //================================================ For each unique pair of C3s
    for ( proshade_unsign c4 = 0; c4 < static_cast<proshade_unsign> ( C4List.size() ); c4++ )
    {
        for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ )
        {
            //======================================== Compare only C3s to the C3List
            const FloatingPoint< proshade_double > lhs1 ( CSymList->at(cSym)[0] ), rhs1 ( 3.0 );
            if ( !lhs1.AlmostEquals ( rhs1 ) ) { continue; }
            
            //========================================  Check the angle between the C4 and C3 axes
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &CSymList->at(C4List.at(c4))[1],
                                                                                                     &CSymList->at(C4List.at(c4))[2],
                                                                                                     &CSymList->at(C4List.at(c4))[3],
                                                                                                     &CSymList->at(cSym)[1],
                                                                                                     &CSymList->at(cSym)[2],
                                                                                                     &CSymList->at(cSym)[3] );
            
            //======================================== Is the angle approximately the dihedral angle
            if ( ( ( 1.0 / sqrt ( 3.0 ) ) > ( dotProduct - axErr ) ) && ( ( 1.0 / sqrt ( 3.0 ) ) < ( dotProduct + axErr ) ) )
            {
                return                                ( true );
            }
        }
    }
    
    //================================================ Done
    return                                            ( false );
    
}

/*! \brief This function takes the list of C symmetries and finds the 3 C4 symmetries with perpendicular angles required for full octahedral symmetry.
 
    This function is specific to detecting the octahedral symmetry. It should be called once octahedral symmetry is suspected (by detecting its dihedral
    angles) and it needs to be fully described. This function specifically searches for the three C4 symmetries which must all be detected in order to fully
    describe octahedral symmetry. If all three are found, the ret vector will contain these as its only four entries, while it will be empty if some of the
    C4 symmetries are not found. The missing symmetry axis detection is implemented as part of this function as well.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret The vector containing all axes required for the octahedral symmetry detected so far.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
    \param[in] verbose How loud the announcments should be?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
    \param[in] minPeakHeight The threshold for peak height.
 */
void ProSHADE_internal_symmetry::findOcta3C4s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C4PossibilitiesHlp;
    std::vector< std::vector< proshade_unsign > > C4Possibilities;
    bool groupMatched;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of three C4 axes.", messageShift );
    
    //================================================ For all symmetries in the C symmetries list
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
    {
        //============================================ Search only using C4s
        if ( CSymList->at(cIt)[0] != 4.0 || CSymList->at(cIt)[5] < minPeakHeight ) { continue; }

        //============================================ If second or more C4, check if it has the correct angle to all other already found C4s for each group
        groupMatched                                  = false;
        for ( proshade_unsign gIt = 0; gIt < static_cast<proshade_unsign> ( C4Possibilities.size() ); gIt++ )
        {
            if ( ProSHADE_internal_symmetry::testGroupAgainstSymmetry ( CSymList, &C4Possibilities.at(gIt), CSymList->at(cIt), axErr, 0.0, true, cIt ) ) { ProSHADE_internal_misc::addToUnsignVector ( &C4Possibilities.at(gIt), cIt ); groupMatched = true; break; }
        }

        //=========================================== If no group matched, create a new group
        if ( !groupMatched ) { C4PossibilitiesHlp.clear(); ProSHADE_internal_misc::addToUnsignVector ( &C4PossibilitiesHlp, cIt ); ProSHADE_internal_misc::addToUnsignVectorVector ( &C4Possibilities, C4PossibilitiesHlp ); continue; }
    }

    //================================================ Test for missing symmetry axes, if need be
    ProSHADE_internal_symmetry::findMissingAxes       ( &C4Possibilities, CSymList, 3, axErr, 0.0, 4, dataObj, minPeakHeight );

    //================================================ Any group has 3 entries? If more such groups, take the one with highest average height.
    proshade_double maxHeight = 0.0; proshade_unsign maxGrp = 0;
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( C4Possibilities.size() ); iter++ ) { if ( C4Possibilities.at(iter).size() == 3 ) { if ( ( ( CSymList->at(C4Possibilities.at(iter).at(0))[5] + CSymList->at(C4Possibilities.at(iter).at(1))[5] + CSymList->at(C4Possibilities.at(iter).at(2))[5] ) / 3.0 ) > maxHeight ) { maxHeight = ( ( CSymList->at(C4Possibilities.at(iter).at(0))[5] + CSymList->at(C4Possibilities.at(iter).at(1))[5] + CSymList->at(C4Possibilities.at(iter).at(2))[5] ) / 3.0 ); maxGrp = iter; } } }
    
    if ( C4Possibilities.at(maxGrp).size() == 3 )
    {
        //============================================ Success! Save and exit
        for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( C4Possibilities.at(maxGrp).size() ); it++ ) { ProSHADE_internal_misc::addToDblPtrVector ( ret, CSymList->at(C4Possibilities.at(maxGrp).at(it)) ); }
        
        //============================================ Report progress
        ProSHADE_internal_messages::printProgressMessage ( verbose, 3, "Detection of three C4 axes successfull.", messageShift );
        
        //============================================ Done
        return ;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes the list of C symmetries and finds the four C3 symmetries with correct angles required for full octahedral symmetry.
 
    This function is specific to detecting the tetrahedral symmetry. It should be called once tetrahedral symmetry is suspected (by detecting its dihedral
    angles) and it needs to be fully described. This function specifically searches for the four C3 symmetries which must all be detected in order to fully
    describe octahedral symmetry. If all four are found, the ret vector will have these four axes added to the already present three C4 axes; alternatively,
    the ret array size will not change. In order not to replicate computations, if tetrahedral symmetry has already been detected, the four axes sought here
    are the same as the first four axes detected there, so simple copying is used instead of re-computing the results anew.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret The vector .
    \param[in] axErr The error tolerance on angle matching.
    \param[in] verbose How loud the announcments should be?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
    \param[in] TetraSymList A vector containing the already detected tetrahedral symmetries - this is to avoid the same search for four C3 symmetry axes.
 */
void ProSHADE_internal_symmetry::findOcta4C3s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C4s, prospectiveC3s, C3PossibilitiesHlp;
    std::vector< std::vector< proshade_unsign > > C3Possibilities;
    proshade_double dotProd;
    bool groupMatched;
    for ( proshade_unsign iter = 0; iter < 3; iter++ ) { ProSHADE_internal_misc::addToUnsignVector ( &C4s, iter ); }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of four C3 axes.", messageShift );
    
    //================================================ For each C4
    for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( ret->size() ); rIt++ )
    {
        //============================================ For each C3, check it has angle ( acos( 1/sqrt(3) ) ) to the tested C4
        for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
        {
            //======================================== Search only using C3s
            if ( CSymList->at(cIt)[0] != 3.0 || CSymList->at(cIt)[5] < minPeakHeight ) { continue; }
            
            //======================================== Check the C3 axis to the C4 ( acos ( 1/sqrt(3) ) )
            dotProd = ProSHADE_internal_maths::computeDotProduct ( &ret->at(rIt)[1], &ret->at(rIt)[2], &ret->at(rIt)[3], &CSymList->at(cIt)[1], &CSymList->at(cIt)[2], &CSymList->at(cIt)[3] );
            
            if ( ( std::abs ( dotProd ) > ( ( 1.0 / sqrt(3.0) ) - axErr ) ) && ( std::abs ( dotProd ) < ( ( 1.0 / sqrt(3.0) ) + axErr ) ) ) { ProSHADE_internal_misc::addToUnsignVector ( &prospectiveC3s, cIt ); }
        }
    }
    
    //================================================ Group the prospective C3s
    C3Possibilities.clear(); C3PossibilitiesHlp.clear();
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( prospectiveC3s.size() ); cIt++ )
    {
        //============================================ If second or more C3, check if it can be placed in any group with having acos (1/3) to all its members
        groupMatched                                  = false;
        for ( proshade_unsign gIt = 0; gIt < static_cast<proshade_unsign> ( C3Possibilities.size() ); gIt++ )
        {
            if ( ProSHADE_internal_symmetry::testGroupAgainstSymmetry ( CSymList, &C3Possibilities.at(gIt), CSymList->at(prospectiveC3s.at(cIt)), axErr, 1.0/3.0, true, prospectiveC3s.at(cIt) ) ) { ProSHADE_internal_misc::addToUnsignVector ( &C3Possibilities.at(gIt), prospectiveC3s.at(cIt) ); groupMatched = true; break; }
        }

        //============================================ If no group matched, create a new group
        if ( !groupMatched ) { C3PossibilitiesHlp.clear(); ProSHADE_internal_misc::addToUnsignVector ( &C3PossibilitiesHlp, prospectiveC3s.at(cIt) ); ProSHADE_internal_misc::addToUnsignVectorVector ( &C3Possibilities, C3PossibilitiesHlp ); continue; }
    }

    //================================================ Find the best group or return empty
    while ( C3Possibilities.size() != 0 )
    {
        //============================================ Test for missing symmetry axes, if need be
        ProSHADE_internal_symmetry::findMissingAxes ( &C3Possibilities, CSymList, 4, axErr, 1.0/3.0, 3, dataObj, minPeakHeight );

        //============================================ Found four C3s?
        if ( C3Possibilities.at(0).size() == 4 )
        {
            //======================================== Success! Save and exit
            for ( proshade_unsign it = 0; it < 4; it++ ) { ProSHADE_internal_misc::addToDblPtrVector ( ret, CSymList->at(C3Possibilities.at(0).at(it)) ); }

            //======================================== Report progress
            ProSHADE_internal_messages::printProgressMessage ( verbose, 3, "Detection of four C3 axes successfull.", messageShift );

            //======================================== Done
            return ;
        }
        else { C3Possibilities.erase                  ( C3Possibilities.begin() ); }
    }

    //================================================ Done
    return ;
    
}

/*! \brief This function takes the list of C symmetries and finds the six C2 symmetries with correct angles required for full octahedral symmetry.
 
    This function is specific to detecting the octahedral symmetry. It should be called once octahedral symmetry is suspected (by detecting its dihedral
    angles) and it needs to be fully described. This function specifically searches for the six C2 symmetries which must all be detected in order to fully
    describe octahedral symmetry. If all six are found, the ret vector will have these six axes added to the already present three C4 axes and the four C3 axes;
    alternatively, the ret array size will not change.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret The vector containing the already detected axes to which newly detected axes (if any) will be added.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] verbose How loud the announcments should be?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
 */
void ProSHADE_internal_symmetry::findOcta6C2s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > prospectiveC2s, retGrp;
    proshade_double dotProd;
    proshade_unsign noPerpendicular, noSqrtTwo;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of six C2 axes.", messageShift );
    
    //================================================ For each C2
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
    {
        //============================================ Use only C2s
        const FloatingPoint< proshade_double > lhs999 ( CSymList->at(cIt)[5] ), rhs999 ( static_cast< proshade_double > ( -999.9 ) );
        if ( CSymList->at(cIt)[0] != 2.0 || ( ( CSymList->at(cIt)[5] < minPeakHeight ) && ! ( lhs999.AlmostEquals( rhs999 ) ) ) ) { continue; }
        
        //============================================ Check the C2 has acos ( 1/sqrt(2) ) to 2 C4s and acos ( 0.0 ) to the third C4
        noPerpendicular = 0; noSqrtTwo = 0;
        for ( proshade_unsign rIt = 0; rIt < 3; rIt++ )
        {
            dotProd                                   = ProSHADE_internal_maths::computeDotProduct ( &ret->at(rIt)[1],
                                                                                                     &ret->at(rIt)[2],
                                                                                                     &ret->at(rIt)[3],
                                                                                                     &CSymList->at(cIt)[1],
                                                                                                     &CSymList->at(cIt)[2],
                                                                                                     &CSymList->at(cIt)[3] );
            
            if ( ( std::abs ( dotProd ) > ( ( 1.0 / sqrt(2.0) ) - axErr ) ) && ( std::abs ( dotProd ) < ( ( 1.0 / sqrt(2.0) ) + axErr ) ) ) { noSqrtTwo       += 1; continue; }
            if ( ( std::abs ( dotProd ) > (   0.0               - axErr ) ) && ( std::abs ( dotProd ) < (   0.0               + axErr ) ) ) { noPerpendicular += 1; continue; }
        }
        
        //============================================ If correct angles distribution is found, save the axis
        if ( ( noSqrtTwo == 2 ) && ( noPerpendicular == 1 ) )
        {
            ProSHADE_internal_misc::addToUnsignVector ( &prospectiveC2s, cIt );
        }
    }
    
    //================================================ Search for missing axes
    for ( proshade_unsign iter = 0; iter < 3; iter++ ) { ProSHADE_internal_misc::addToUnsignVector ( &retGrp, iter ); }
    if ( !ProSHADE_internal_symmetry::findMissingAxesDual ( &prospectiveC2s, CSymList, ret, &retGrp, 6, axErr, 1, 0.0, 2, 1/sqrt(2.0), 2, dataObj ) )
    {
        return ;
    }
    
    //================================================ Found correct number of axes! Now save the
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( prospectiveC2s.size() ); iter++ )
    {
        ProSHADE_internal_misc::addToDblPtrVector     ( ret, CSymList->at(prospectiveC2s.at(iter)) );
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "Detection of six C2 axes successfull.", messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function tries to find a particular symmetry axes which would complete a group of symmetries with two different angle requirement to another group.

    This function takes a list of axes to which a new axis should have two particular angles (to two different group members, that is). It then uses algebraic solution finding approach
    to compute possible solutions which would satisfy this condition, testing whether such solutions comply with the appropriate number of angles to number of members and for
    the new solutions being unique. If the required number of solutions is found, it will add the newly detected solutions to the CSymList vector and update the possibilities indices list,
    otherwise it will leave both alone.

    \param[in] possibilities A vector of already detected axis indices which should be extended.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret A list of already detected octahedral axes.
    \param[in] retGroup A vector of indices in the ret list which form the group to which new axes are compared to.
    \param[in] requiredNoAxes Number of axes required for positive result.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] noMatchesG1 The number of axes from ret that need to be matched with angle1.
    \param[in] angle1 The angle with which noMatchesG1 axes need to be matched with the retGroup axes.
    \param[in] noMatchesG2 The number of axes from ret that need to be matched with angle2.
    \param[in] angle2 The angle with which noMatchesG2 axes need to be matched with the retGroup axes.
    \param[in] fold The fold of the searched for axis.
    \param[in] dataObj The full data holding object pointer - this is to get access to self-rotation function values.
    \param[out] atLeastOne Boolean value speciying whether at least the minimum required number of axes was found.
*/
bool ProSHADE_internal_symmetry::findMissingAxesDual ( std::vector< proshade_unsign >* possibilities, std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, std::vector< proshade_unsign >* retGroup, proshade_unsign requiredNoAxes, proshade_double axErr, proshade_unsign noMatchesG1, proshade_double angle1, proshade_unsign noMatchesG2, proshade_double angle2, proshade_unsign fold, ProSHADE_internal_data::ProSHADE_data* dataObj )
{
    //================================================ Initialise variables
    bool atLeastOne                                   = false;
    std::vector< proshade_double* > prosp;
    std::vector< proshade_double > sol;
    
    //================================================ Proceed only if need be
    if ( static_cast<proshade_unsign> ( possibilities->size() ) == requiredNoAxes ) { atLeastOne = true; return ( atLeastOne ); }
    
    //================================================ Copy already found to prospective
    for ( proshade_unsign prIt = 0; prIt < static_cast<proshade_unsign> ( possibilities->size() ); prIt++ )
    {
        ProSHADE_internal_symmetry::addAxisUnlessSame ( static_cast< proshade_unsign > ( CSymList->at(possibilities->at(prIt))[0] ),
                                                        CSymList->at(possibilities->at(prIt))[1],
                                                        CSymList->at(possibilities->at(prIt))[2],
                                                        CSymList->at(possibilities->at(prIt))[3],
                                                        CSymList->at(possibilities->at(prIt))[5], &prosp, axErr );
    }
    
    //================================================ Start generating possible solutions
    for ( proshade_unsign rgIt1 = 0; rgIt1 < static_cast<proshade_unsign> ( retGroup->size() ); rgIt1++ )
    {
        for ( proshade_unsign rgIt2 = 0; rgIt2 < static_cast<proshade_unsign> ( retGroup->size() ); rgIt2++ )
        {
            //======================================== Use unique combinations (order matters here!)
            if ( rgIt1 == rgIt2 ) { continue; }
            
            //======================================== Generate possible solution (1)
            sol                                       = ProSHADE_internal_maths::findVectorFromTwoVAndTwoD ( ret->at(rgIt1)[1], ret->at(rgIt1)[2], ret->at(rgIt1)[3],
                                                                                                             ret->at(rgIt2)[1], ret->at(rgIt2)[2], ret->at(rgIt2)[3], angle1, angle2 );
            
            //======================================== Check if solution fits the group completely
            ProSHADE_internal_symmetry::checkFittingAxisDualAndSave ( retGroup, ret, fold, sol.at(0), sol.at(1), sol.at(2), &prosp, axErr, noMatchesG1, angle1, noMatchesG2, angle2, dataObj );
            if ( prosp.size() == requiredNoAxes ) { break; }
            
            //======================================== Generate possible solution (2)
            sol                                       = ProSHADE_internal_maths::findVectorFromTwoVAndTwoD ( ret->at(rgIt1)[1], ret->at(rgIt1)[2], ret->at(rgIt1)[3],
                                                                                                             ret->at(rgIt2)[1], ret->at(rgIt2)[2], ret->at(rgIt2)[3], -angle1, -angle2 );
            
            //======================================== Check if solution fits the group completely
            ProSHADE_internal_symmetry::checkFittingAxisDualAndSave ( retGroup, ret, fold, sol.at(0), sol.at(1), sol.at(2), &prosp, axErr, noMatchesG1, angle1, noMatchesG2, angle2, dataObj );
            if ( prosp.size() == requiredNoAxes ) { break; }
        }
        
        if ( prosp.size() == requiredNoAxes ) { break; }
    }

    //================================================ Found all required axes!
    if ( static_cast<proshade_unsign> ( prosp.size() ) == requiredNoAxes )
    {
        //============================================ Copy the detected axes
        for ( proshade_unsign iter = static_cast<proshade_unsign> ( possibilities->size() ); iter < static_cast<proshade_unsign> ( prosp.size() ); iter++ )
        {
            if ( ProSHADE_internal_maths::isAxisUnique ( CSymList, prosp.at(iter), axErr ) )
            {
                //==================================== Add
                ProSHADE_internal_misc::addToUnsignVector ( possibilities, static_cast< proshade_unsign > ( CSymList->size() ) );
                ProSHADE_internal_misc::addToDblPtrVector ( CSymList, prosp.at(iter) );
            }
        }
        
        //============================================ Done
        atLeastOne                                    = true;
        return                                        ( atLeastOne );
    }
    else
    {
        //============================================ Delete the created, but not used axes
        for ( proshade_unsign iter = static_cast<proshade_unsign> ( possibilities->size() ); iter < static_cast<proshade_unsign> ( prosp.size() ); iter++ )
        {
            delete[] prosp.at(iter);
        }
    }
    
    //================================================ Done
    return                                            ( atLeastOne );
    
}

/*! \brief This function simply creates a new axis from information in aruments and tests if no such axis already exists, saving it if need be.
 
    This is a simple helper function, which takes all the new axis information and creates the ProSHADE axis representation from these. It then proceeds to check
    if such axis does not already exist in the supplied vector, if not, it saves the new axis; alternatively, it just discards the created axis and terminates.
 
    \param[in] fold The fold of the searched for axis.
    \param[in] axX The x-axis element of the new axis.
    \param[in] axY The y-axis element of the new axis.
    \param[in] axZ The z-axis element of the new axis.
    \param[in] axHeight The average peak height of the new axis.
    \param[in] averageFSC The value of average FSC, if computed - otherwise enter -1.0.
    \param[in] prosp The vector to which the axis is to be saved.
    \param[in] axErr The error tolerance on angle matching.
    \param[out] addedNo Position at which the symmetry either already is, or was added to.
 */
proshade_signed ProSHADE_internal_symmetry::addAxisUnlessSame ( proshade_unsign fold, proshade_double axX, proshade_double axY, proshade_double axZ, proshade_double axHeight, proshade_double averageFSC, std::vector< proshade_double* >* prosp, proshade_double axErr )
{
    //================================================ Initialise variables
    proshade_double* symHlp                           = new proshade_double[7];
    ProSHADE_internal_misc::checkMemoryAllocation     ( symHlp, __FILE__, __LINE__, __func__ );
    proshade_signed ret                               = -1;
    
    //================================================ Fill in the prospective axis
    symHlp[0]                                         = static_cast<proshade_double> ( fold );
    symHlp[1]                                         = axX;
    symHlp[2]                                         = axY;
    symHlp[3]                                         = axZ;
    symHlp[4]                                         = 2.0 * M_PI / symHlp[0];
    symHlp[5]                                         = axHeight;
    symHlp[6]                                         = averageFSC;
    
    //================================================ If it is not the same as already saved axes
    if ( !ProSHADE_internal_symmetry::isSymmetrySame ( prosp, symHlp, axErr, &ret, averageFSC ) )
    {
        ProSHADE_internal_misc::addToDblPtrVector     ( prosp, symHlp );
        return                                        ( static_cast< proshade_signed > ( prosp->size() - 1 ) );
    }
    else
    {
        delete[] symHlp;
        return                                        ( ret );
    }
    
    //================================================ Done
    
}

/*! \brief This function simply creates a new axis from information in aruments and tests if no such axis already exists, saving it if need be.
 
    This is a simple helper function, which takes all the new axis information and creates the ProSHADE axis representation from these. It then proceeds to check
    if such axis does not already exist in the supplied vector, if not, it saves the new axis; alternatively, it just discards the created axis and terminates.
 
    \param[in] fold The fold of the searched for axis.
    \param[in] axX The x-axis element of the new axis.
    \param[in] axY The y-axis element of the new axis.
    \param[in] axZ The z-axis element of the new axis.
    \param[in] axHeight The average peak height of the new axis.
    \param[in] prosp The vector to which the axis is to be saved.
    \param[in] axErr The error tolerance on angle matching.
    \param[out] addedNo Position at which the symmetry either already is, or was added to.
 */
proshade_signed ProSHADE_internal_symmetry::addAxisUnlessSame ( proshade_unsign fold, proshade_double axX, proshade_double axY, proshade_double axZ, proshade_double axHeight, std::vector< proshade_double* >* prosp, proshade_double axErr )
{
    //================================================ Initialise variables
    proshade_double* symHlp                           = new proshade_double[7];
    ProSHADE_internal_misc::checkMemoryAllocation     ( symHlp, __FILE__, __LINE__, __func__ );
    proshade_signed ret                               = -1;
    
    //================================================ Fill in the prospective axis
    symHlp[0]                                         = static_cast<proshade_double> ( fold );
    symHlp[1]                                         = axX;
    symHlp[2]                                         = axY;
    symHlp[3]                                         = axZ;
    symHlp[4]                                         = 2.0 * M_PI / symHlp[0];
    symHlp[5]                                         = axHeight;
    symHlp[6]                                         = -std::numeric_limits < proshade_double >::infinity();
    
    //================================================ If it is not the same as already saved axes
    if ( !ProSHADE_internal_symmetry::isSymmetrySame ( prosp, symHlp, axErr, &ret ) )
    {
        ProSHADE_internal_misc::addToDblPtrVector     ( prosp, symHlp );
        return                                        ( static_cast< proshade_signed > ( prosp->size() - 1 ) );
    }
    else
    {
        delete[] symHlp;
        return                                        ( ret );
    }
    
    //================================================ Done
    
}

/*! \brief This function takes a newly detected "missing" axis and tests it for belonging to the group, checking the height and replacing lower height members with better members.
 
    This function takes the list of already detected axes, information about the tested new axis and the conditions for belonging. It then proceeds to check if the
    new axis conforms to the conditions of belonging. If so, it then checks if the axis height is high enough to be considered as part of the group. Again, if so,
    it will save this new axis to the old set, replacing any old axis with this new one, if it is the same and has better height.
 
    \param[in] retGroup A vector of indices in the ret list which form the group to which new axes are compared to.
    \param[in] ret A list of already detected axes.
    \param[in] fold The fold of the searched for axis.
    \param[in] axX The x-axis element of the new axis.
    \param[in] axY The y-axis element of the new axis.
    \param[in] axZ The z-axis element of the new axis.
    \param[in] prosp The vector to which the axis is to be saved.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] noMatchesG1 The number of axes from ret that need to be matched with angle1.
    \param[in] angle1 The angle with which noMatchesG1 axes need to be matched with the retGroup axes.
    \param[in] noMatchesG2 The number of axes from ret that need to be matched with angle2.
    \param[in] angle2 The angle with which noMatchesG2 axes need to be matched with the retGroup axes.
    \param[in] dataObj The full data holding object pointer - this is to get access to self-rotation function values.
    \param[out] Bool True if the axis was added to the group, false otherwise.
 */
bool ProSHADE_internal_symmetry::checkFittingAxisDualAndSave ( std::vector< proshade_unsign >* retGroup, std::vector< proshade_double* >* ret, proshade_unsign fold, proshade_double axX, proshade_double axY, proshade_double axZ, std::vector< proshade_double* >* prosp, proshade_double axErr, proshade_unsign noMatchesG1, proshade_double angle1, proshade_unsign noMatchesG2, proshade_double angle2, ProSHADE_internal_data::ProSHADE_data* dataObj )
{
    //================================================ Initialise variables
    proshade_unsign noG1                              = 0;
    proshade_unsign noG2                              = 0;
    proshade_double dotProd                           = 0.0;
    proshade_double axHeight                          = 0.0;
    
    //================================================ Find the angle and count dual matching frequencies
    for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( retGroup->size() ); rIt++ )
    {
        dotProd                                       = ProSHADE_internal_maths::computeDotProduct ( &ret->at(retGroup->at(rIt))[1],
                                                                                                     &ret->at(retGroup->at(rIt))[2],
                                                                                                     &ret->at(retGroup->at(rIt))[3],
                                                                                                     &axX, &axY, &axZ );
        
        if ( ( std::abs ( dotProd ) > ( angle1 - axErr ) ) && ( std::abs ( dotProd ) < ( angle1 + axErr ) ) ) { noG1 += 1; continue; }
        if ( ( std::abs ( dotProd ) > ( angle2 - axErr ) ) && ( std::abs ( dotProd ) < ( angle2 + axErr ) ) ) { noG2 += 1; continue; }
    }
    
    //================================================ If correct frequencies are matched, check height.
    if ( ( noG1 == noMatchesG1 ) && ( noG2 == noMatchesG2 ) )
    {
        //============================================ Is the height good enough?
        axHeight                                      = ProSHADE_internal_symmetry::missingAxisHeight ( axX, axY, axZ, dataObj, fold, axErr );
        
        //============================================ If so, save
        if ( axHeight > 0.1 )
        {
            proshade_unsign prevProsp                 = static_cast<proshade_unsign> ( prosp->size() );
            ProSHADE_internal_symmetry::addAxisUnlessSame ( fold, axX, axY, axZ, axHeight, prosp, axErr );
            
            if ( static_cast<proshade_unsign> ( prosp->size() ) > prevProsp ) { return ( true ); }
            else                                                              { return ( false ); }
        }
    }
    
    //================================================ Done
    return                                            ( false );
    
}

/*! \brief This function obtains a list of all I symmetry axes from the already computed C symmetries list.
 
    This function starts by checking if there is a pair of C3 and C4 symmetries with the octahedron dihedral angle ( acos ( 1/sqrt(3) ) ). If so, it will
    then assume existence of octahedral symmetry and it will search for three C4 axes, four C3 axes and six C2 axes with the correct angle to each other
    and within the group. If all required axes are detected, it will return a list of 13 axes, otherwise it will return empty or shorter list. Automated
    missing symmetry axis detection is also included.
 
    \param[in] settings A pointer to settings class containing all the information required for symmetry detection.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
 */
std::vector< proshade_double* > ProSHADE_internal_data::ProSHADE_data::getIcosahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList )
{
    //================================================ Initialise variables
    std::vector< proshade_double* > ret;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting I symmetry detection.", settings->messageShift );
    
    //================================================ Are the basic requirements for icosahedral symmetry met?
    if ( ProSHADE_internal_symmetry::detectIcosahedralSymmetry ( CSymList, settings->axisErrTolerance, settings->minSymPeak ) )
    {
        //============================================ Search for all the symmetry axes
        ProSHADE_internal_symmetry::findIcos6C5s      ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->messageShift, settings->minSymPeak );
        if ( ret.size() != 6 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }

        ProSHADE_internal_symmetry::findIcos10C3s     ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->messageShift, settings->minSymPeak );
        if ( ret.size() != 16 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }

        ProSHADE_internal_symmetry::findIcos15C2s     ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->messageShift, settings->minSymPeak );
        if ( ret.size() != 31 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        else
        {
            //======================================== Sort ret by fold
            std::sort                                 ( ret.begin(), ret.end(), ProSHADE_internal_misc::sortSymInvFoldHlp );
            
            //======================================== Save indices
            for ( proshade_unsign csIt = 0; csIt < static_cast<proshade_unsign> ( CSymList->size() ); csIt++ )
            {
                for ( proshade_unsign retIt = 0; retIt < static_cast<proshade_unsign> ( ret.size() ); retIt++ )
                {
                    const FloatingPoint< proshade_double > lhs1 ( CSymList->at(csIt)[0] ), rhs1 ( ret.at(retIt)[0] );
                    const FloatingPoint< proshade_double > lhs2 ( CSymList->at(csIt)[1] ), rhs2 ( ret.at(retIt)[1] );
                    const FloatingPoint< proshade_double > lhs3 ( CSymList->at(csIt)[2] ), rhs3 ( ret.at(retIt)[2] );
                    const FloatingPoint< proshade_double > lhs4 ( CSymList->at(csIt)[3] ), rhs4 ( ret.at(retIt)[3] );
                    const FloatingPoint< proshade_double > lhs5 ( CSymList->at(csIt)[4] ), rhs5 ( ret.at(retIt)[4] );
                    const FloatingPoint< proshade_double > lhs6 ( CSymList->at(csIt)[5] ), rhs6 ( ret.at(retIt)[5] );
                    if ( lhs1.AlmostEquals ( rhs1 ) &&
                         lhs2.AlmostEquals ( rhs2 ) &&
                         lhs3.AlmostEquals ( rhs3 ) &&
                         lhs4.AlmostEquals ( rhs4 ) &&
                         lhs5.AlmostEquals ( rhs5 ) &&
                         lhs6.AlmostEquals ( rhs6 ) )
                    {
                        ProSHADE_internal_misc::addToUnsignVector ( &settings->allDetectedIAxes, csIt );
                    }
                }
            }
        }
    }

    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "I symmetry detection complete.", settings->messageShift );

    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function predicts a list of all I symmetry axes from the already computed C symmetries list.
 
    This function starts by checking if there is a pair of C3 and C5 symmetries with the icosahedron dihedral angle ( acos( std::sqrt ( ( 1.0 + 2.0 / std::sqrt ( 5.0 ) ) / 3.0 ) ) ). If found,
    it calls the  predictIcosAxes() function, which uses the knowledge of the two axes (C5 and C3) which are closest to the dihedral angle to find the best rotation matrix matching a
    pre-computed icosahedron model to the detected axes. After rotating the model, the model axes become the predicted axes for the structure and their peak heights are then
    computed. Once complete, all the predicted axes are in the ret variable.
 
    \warning This function does not check if the correct number of C axes was found, it is assumed this will be checked when the determination of
    which symmetry was found.
 
    \param[in] settings A pointer to settings class containing all the information required for symmetry detection.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[out] ret A vector of all the detected axes in the standard ProSHADE format with height either the detected value (for the detected ones) or 0 for the predicted ones.
 */
std::vector < std::vector< proshade_double* > > ProSHADE_internal_data::ProSHADE_data::getPredictedIcosahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList )
{
    //================================================ Initialise variables
    std::vector< std::vector< proshade_double* > > ret;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting I symmetry prediction.", settings->messageShift );
    
    //================================================ Are the basic requirements for icosahedral symmetry met?
    if ( ProSHADE_internal_symmetry::detectIcosahedralSymmetry ( CSymList, settings->axisErrTolerance, settings->minSymPeak ) )
    {
        //============================================ Generate the rest of the axes
        ProSHADE_internal_symmetry::predictIcosAxes   ( CSymList, &ret, settings->axisErrTolerance, settings->minSymPeak );

        //============================================ For each possible axes pair
        for ( size_t pIt = 0; pIt < ret.size(); pIt++ )
        {
            //======================================== Get heights for the predicted axes
            ProSHADE_internal_symmetry::findPredictedAxesHeights ( &(ret.at(pIt)), this, settings );
        }
    }
    
    //================================================ Sort by best peak height sum
    std::sort                                         ( ret.begin(), ret.end(), ProSHADE_internal_misc::sortISymByPeak );

    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "I symmetry prediction complete.", settings->messageShift );

    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function predicts a list of all O symmetry axes from the already computed C symmetries list.
 
    This function starts by checking if there is a pair of C3 and C4 symmetries with the octahedron dihedral angle ( acos( ( 1.0 / sqrt ( 3.0 ) ) ) ). If found,
    it calls the  predictOctaAxes() function, which uses the knowledge of the two axes (C4 and C3) which are closest to the dihedral angle to find the best rotation matrix matching a
    pre-computed octahedron model to the detected axes. After rotating the model, the model axes become the predicted axes for the structure and their peak heights are then
    computed. Once complete, all the predicted axes are in the ret variable.
 
    \warning This function does not check if the correct number of C axes was found, it is assumed this will be checked when the determination of
    which symmetry was found.
 
    \param[in] settings A pointer to settings class containing all the information required for symmetry detection.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[out] ret A vector of all the detected axes in the standard ProSHADE format with height either the detected value (for the detected ones) or 0 for the predicted ones.
 */
std::vector< proshade_double* > ProSHADE_internal_data::ProSHADE_data::getPredictedOctahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList )
{
    //================================================ Initialise variables
    std::vector< proshade_double* > ret;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting O symmetry prediction.", settings->messageShift );
    
    //================================================ Are the basic requirements for icosahedral symmetry met?
    if ( ProSHADE_internal_symmetry::detectOctahedralSymmetry ( CSymList, settings->axisErrTolerance, settings->minSymPeak ) )
    {
        //============================================ Generate the rest of the axes
        ProSHADE_internal_symmetry::predictOctaAxes   ( CSymList, &ret, settings->axisErrTolerance, settings->minSymPeak );

        //============================================ Get heights for the predicted axes
        ProSHADE_internal_symmetry::findPredictedAxesHeights ( &ret, this, settings );

        //============================================ Add predicted axes to detected C axes list and also to the settings Icosahedral symmetry list
        for ( proshade_unsign retIt = 0; retIt < static_cast < proshade_unsign > ( ret.size() ); retIt++ )
        {
            proshade_signed matchedPos                = ProSHADE_internal_symmetry::addAxisUnlessSame ( static_cast< proshade_unsign > ( ret.at(retIt)[0] ), ret.at(retIt)[1], ret.at(retIt)[2], ret.at(retIt)[3], ret.at(retIt)[5], CSymList, settings->axisErrTolerance );
            ProSHADE_internal_misc::addToUnsignVector ( &settings->allDetectedOAxes, static_cast < proshade_unsign > ( matchedPos ) );
        }
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "O symmetry prediction complete.", settings->messageShift );

    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function takes the list of C symmetries and decides whether basic requirements for isosahedral symmetry are there.
 
    This function first finds all the C5 symmetries in the C symmetries list and then it checks each present C5 against all C3 symmetries for having
    the angle between the pair equal to the dihedral angle of an icosahedron ( acos( sqrt(5)/3 ) ). If a single such pair is detected, this is
    likely an icosahedral symmetry and all other axes need to be located. Otherwise, false is returned.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] minPeakHeight The minimum average peak height required for symmetry axis to be considered.
    \param[out] X Boolean value telling whether there are C5 and C3 symmetries with icosahedral dihhedral angle.
 */
bool ProSHADE_internal_symmetry::detectIcosahedralSymmetry ( std::vector< proshade_double* >* CSymList, proshade_double axErr, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C5List;
    proshade_double dotProduct;
    
    //================================================ Find all C5 symmetries
    for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ )
    {
        const FloatingPoint< proshade_double > lhs1 ( CSymList->at(cSym)[0] ), rhs1 ( 5.0 );
        if ( lhs1.AlmostEquals ( rhs1 ) && CSymList->at(cSym)[5] >= minPeakHeight ) { ProSHADE_internal_misc::addToUnsignVector ( &C5List, cSym ); }
    }
    
    //================================================ For each unique pair of C5 and C3
    for ( proshade_unsign c5 = 0; c5 < static_cast<proshade_unsign> ( C5List.size() ); c5++ )
    {
        for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ )
        {
            //======================================== Compare only C3s to the C5List
            const FloatingPoint< proshade_double > lhs1 ( CSymList->at(cSym)[0] ), rhs1 ( 3.0 );
            if ( !lhs1.AlmostEquals ( rhs1 ) ) { continue; }
            
            //========================================  Check the angle between the C5 and C3 axes
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &CSymList->at(C5List.at(c5))[1],
                                                                                                     &CSymList->at(C5List.at(c5))[2],
                                                                                                     &CSymList->at(C5List.at(c5))[3],
                                                                                                     &CSymList->at(cSym)[1],
                                                                                                     &CSymList->at(cSym)[2],
                                                                                                     &CSymList->at(cSym)[3] );
            
            //======================================== Is the angle approximately the dihedral angle
            if ( std::abs ( std::abs( std::sqrt ( ( 1.0 + 2.0 / std::sqrt ( 5.0 ) ) / 3.0 ) ) - std::abs( dotProduct ) ) < axErr )
            {
                return                                ( true );
            }
        }
    }
    
    //================================================ Done
    return                                            ( false );
    
}

/*! \brief This function takes the list of C symmetries and finds the six C5 symmetries with given angles required for full icosahedral symmetry.
 
    This function searches the list of all detected C symmetries for the presence of six C5 symmetries, which have the angle of acos (0.5) to each other; this
    ability is specifically required for detection of icosahedral symmetry. This function allows for multiple groups of C5 symmetries, doing the missing symmetry
    axis checks and returning the group with highest average peak height. If successfull, the ret vector will have 6 entries, otherwise it will be empty.
 
    This function is specific to detecting the octahedral symmetry. It should be called once octahedral symmetry is suspected (by detecting its dihedral
    angles) and it needs to be fully described. This function specifically searches for the three C4 symmetries which must all be detected in order to fully
    describe octahedral symmetry. If all three are found, the ret vector will contain these as its only four entries, while it will be empty if some of the
    C4 symmetries are not found. The missing symmetry axis detection is implemented as part of this function as well.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret The vector .
    \param[in] axErr The error tolerance on angle matching.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
    \param[in] verbose How loud the announcments should be?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
 */
void ProSHADE_internal_symmetry::findIcos6C5s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C5PossibilitiesHlp;
    std::vector< std::vector< proshade_unsign > > C5Possibilities;
    bool groupMatched;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of six C5 axes.", messageShift );

    //================================================ For all symmetries in the C symmetries list
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
    {
        //============================================ Search only using C5s and check peak height
        const FloatingPoint< proshade_double > lhs1 ( CSymList->at(cIt)[0] ), rhs1 ( 5.0 );
        if ( !lhs1.AlmostEquals ( rhs1 ) || CSymList->at(cIt)[5] < minPeakHeight ) { continue; }

        //============================================ If second or more C5, check if it has the correct angle to all other already found C5s for each group
        groupMatched                                  = false;
        for ( proshade_unsign gIt = 0; gIt < static_cast<proshade_unsign> ( C5Possibilities.size() ); gIt++ )
        {
            if ( ProSHADE_internal_symmetry::testGroupAgainstSymmetry ( CSymList, &C5Possibilities.at(gIt), CSymList->at(cIt), axErr, 1.0/2.0, true, cIt ) ) { ProSHADE_internal_misc::addToUnsignVector ( &C5Possibilities.at(gIt), cIt ); groupMatched = true; break; }
        }

        //============================================ If no group matched, create a new group
        if ( !groupMatched ) { C5PossibilitiesHlp.clear(); ProSHADE_internal_misc::addToUnsignVector ( &C5PossibilitiesHlp, cIt ); ProSHADE_internal_misc::addToUnsignVectorVector ( &C5Possibilities, C5PossibilitiesHlp ); continue; }
    }
    
    //================================================ Test for missing symmetry axes, if need be
    ProSHADE_internal_symmetry::findMissingAxes       ( &C5Possibilities, CSymList, 6, axErr, 1.0 / 2.0, 5, dataObj, minPeakHeight );

    //=================================================Any group has 6 entries? If more such groups, take the one with highest average height.
    proshade_double maxHeight = 0.0; proshade_unsign maxGrp = 0;
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( C5Possibilities.size() ); iter++ ) { if ( C5Possibilities.at(iter).size() == 6 ) { if ( ( ( CSymList->at(C5Possibilities.at(iter).at(0))[5] + CSymList->at(C5Possibilities.at(iter).at(1))[5] + CSymList->at(C5Possibilities.at(iter).at(2))[5] + CSymList->at(C5Possibilities.at(iter).at(3))[5] + CSymList->at(C5Possibilities.at(iter).at(4))[5] + CSymList->at(C5Possibilities.at(iter).at(5))[5] ) / 6.0 ) > maxHeight ) { maxHeight = ( ( CSymList->at(C5Possibilities.at(iter).at(0))[5] + CSymList->at(C5Possibilities.at(iter).at(1))[5] + CSymList->at(C5Possibilities.at(iter).at(2))[5] + CSymList->at(C5Possibilities.at(iter).at(3))[5] + CSymList->at(C5Possibilities.at(iter).at(4))[5] + CSymList->at(C5Possibilities.at(iter).at(5))[5] ) / 6.0 ); maxGrp = iter; } } }

    if ( C5Possibilities.at(maxGrp).size() == 6 )
    {
        //============================================ Success! Save and exit
        for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( C5Possibilities.at(maxGrp).size() ); it++ ) { ProSHADE_internal_misc::addToDblPtrVector ( ret, CSymList->at(C5Possibilities.at(maxGrp).at(it)) ); }

        //============================================ Report progress
        ProSHADE_internal_messages::printProgressMessage ( verbose, 3, "Detection of six C5 axes successfull.", messageShift );

        //============================================ Done
        return ;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds all the pairs of axes conforming to the icosahedron dihedral angle.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
    \param[in] axErr The error tolerance on angle matching.
    \param[out] ret Vector of pairs containing the indices of all axes conforming to the required icosahedron dihedral angle.
 */
std::vector < std::pair< proshade_unsign, proshade_unsign > > findBestIcosDihedralPair ( std::vector< proshade_double* >* CSymList, proshade_double minPeakHeight, proshade_double axErr )
{
    //================================================ Initialise variables
    std::vector < std::pair< proshade_unsign, proshade_unsign > > ret;
    std::vector< proshade_unsign  > C5List;
    proshade_double dotProduct;
    
    //================================================ Find all C5 symmetries
    for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ ) { const FloatingPoint< proshade_double > lhs1 ( CSymList->at(cSym)[0] ), rhs1 ( 5.0 ); if ( lhs1.AlmostEquals ( rhs1 ) && CSymList->at(cSym)[5] >= minPeakHeight ) { ProSHADE_internal_misc::addToUnsignVector ( &C5List, cSym ); } }
    
    //================================================ For each unique pair of C5 and C3
    for ( proshade_unsign c5 = 0; c5 < static_cast<proshade_unsign> ( C5List.size() ); c5++ )
    {
        for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ )
        {
            //======================================== Compare only C3s to the C5List and only with decent average peak height
            const FloatingPoint< proshade_double > lhs1 ( CSymList->at(cSym)[0] ), rhs1 ( 3.0 );
            if ( !lhs1.AlmostEquals ( rhs1 ) ) { continue; }
            if ( CSymList->at(cSym)[5] < minPeakHeight ) { continue; }
            
            //========================================  Check the angle between the C5 and C3 axes
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &CSymList->at(C5List.at(c5))[1],
                                                                                                     &CSymList->at(C5List.at(c5))[2],
                                                                                                     &CSymList->at(C5List.at(c5))[3],
                                                                                                     &CSymList->at(cSym)[1],
                                                                                                     &CSymList->at(cSym)[2],
                                                                                                     &CSymList->at(cSym)[3] );
            
            //======================================== Is the angle approximately the dihedral angle?
            if ( std::abs ( std::abs( std::sqrt ( ( 1.0 + 2.0 / std::sqrt ( 5.0 ) ) / 3.0 ) ) - std::abs( dotProduct ) ) < axErr )
            {
                std::pair< proshade_unsign, proshade_unsign > hlp;
                hlp.first                             = C5List.at(c5);
                hlp.second                            = cSym;
                ret.emplace_back                      ( hlp );
            }
        }
    }
    
    //================================================ Done
    return                                            ( ret );
}

/*! \brief This function predicts all possible icosahedral point groups symmetry axes from the cyclic point groups list.
 
    This function starts by finding the rotation matrix corresponding to the angle between the predicted C5 axis and the C5 axis found in the pre-computed
    Icosahedron model available in the ProSHADE_precomputedValues file. It then proceeds to use this rotation matrix to rotate the pre-computed model C3
    axis to now be in the correct orientation to the detected C5 axis.
 
    Next, the function computes the rotation matrix corresponding to rotation along the detected C5 axis about the angle between the rotated pre-computed model
    C3 axis and the detected C3 axis. Finally, when these two rotation matrices are combined, the resulting rotation matrix is the optimal match rotation between the
    pre-computed model and the detected axes positions. This final rotation matrix is then used to rotate the model axes and these rotated model axes are then the
    predicted axes in the structre.
 
    Please note that the peak heights are set to 0.0 for all predicted axes, as they were not detected in the structure, but were only predicted.
 
    \warning This function assumes that the detectIcosahedralSymmetry() function has successfully run (i.e. returned true).
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret The vector containing all the axes forming icosahedral group or empty vector.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
 */
void ProSHADE_internal_symmetry::predictIcosAxes ( std::vector< proshade_double* >* CSymList, std::vector< std::vector< proshade_double* > >* ret, proshade_double axErr, proshade_double minPeakHeight )
{
    //================================================ Find the best axis combination with dihedral angle and correct folds
    std::vector < std::pair< proshade_unsign, proshade_unsign > > initAxes = findBestIcosDihedralPair ( CSymList, minPeakHeight, axErr );
    
    //================================================ For each pair of possible axis combinations
    for ( size_t pIt = 0; pIt < initAxes.size(); pIt++ )
    {
        //============================================ Create the tetrahedronAxes object
        ProSHADE_internal_precomputedVals::icosahedronAxes *icoAx = new ProSHADE_internal_precomputedVals::icosahedronAxes ( );
        
        //============================================ Find rotation between the detected C5 and the model C5 axes.
        proshade_double* rotMat                       = ProSHADE_internal_maths::findRotMatMatchingVectors ( icoAx->getValue ( 0, 1 ),
                                                                                                             icoAx->getValue ( 0, 2 ),
                                                                                                             icoAx->getValue ( 0, 3 ),
                                                                                                             CSymList->at(initAxes.at(pIt).first)[1],
                                                                                                             CSymList->at(initAxes.at(pIt).first)[2],
                                                                                                             CSymList->at(initAxes.at(pIt).first)[3] );
        
        //============================================ Rotate the model C3 to the correct orientation relative to the detected C5 axis.
        proshade_double* rotModelC3                   = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMat,
                                                                                                                        icoAx->getValue ( 6, 1 ),
                                                                                                                        icoAx->getValue ( 6, 2 ),
                                                                                                                        icoAx->getValue ( 6, 3 ) );
        
        //============================================ Find the angle betwen the rotated model C3 and the detected C3 axes along the detected C5 axis.
        proshade_double bestAng = 0.0, curAngDist, bestAngDist = 999.9;
        proshade_double* rotMatHlp                    = new proshade_double[9];
        ProSHADE_internal_misc::checkMemoryAllocation ( rotMatHlp, __FILE__, __LINE__, __func__ );
        for ( proshade_double ang = 0.0; ang < ( M_PI * 2.0 ); ang += 0.002 )
        {
            //============================================ Compute rotation matrix for this angle value
            ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( rotMatHlp, CSymList->at(initAxes.at(pIt).first)[1], CSymList->at(initAxes.at(pIt).first)[2], CSymList->at(initAxes.at(pIt).first)[3], ang );

            //======================================== Rotate the rotated C2 by the matrix
            proshade_double* rotRotModelC3            = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMatHlp,
                                                                                                                        rotModelC3[0],
                                                                                                                        rotModelC3[1],
                                                                                                                        rotModelC3[2] );

            //======================================== Find distance
            curAngDist                                = std::sqrt ( std::pow ( rotRotModelC3[0] - CSymList->at(initAxes.at(pIt).second)[1], 2.0 ) +
                                                                    std::pow ( rotRotModelC3[1] - CSymList->at(initAxes.at(pIt).second)[2], 2.0 ) +
                                                                    std::pow ( rotRotModelC3[2] - CSymList->at(initAxes.at(pIt).second)[3], 2.0 ) );

            //======================================== Save best angle
            if ( curAngDist < bestAngDist ) { bestAngDist = curAngDist; bestAng = ang; }

            //======================================== Release memory
            delete[] rotRotModelC3;
        }
        
        //============================================ Release memory
        delete[] rotMatHlp;
        
        //============================================ For the rotation matrix along the detected C5 axis with the same anlge as is between the rotated model C3 and the detected C3 axes.
        proshade_double* rotMat2                      = new proshade_double[9];
        ProSHADE_internal_misc::checkMemoryAllocation ( rotMat2, __FILE__, __LINE__, __func__ );
        ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( rotMat2, CSymList->at(initAxes.at(pIt).first)[1], CSymList->at(initAxes.at(pIt).first)[2], CSymList->at(initAxes.at(pIt).first)[3], bestAng );
        
        //============================================ Combine the two rotation matrices into a single rotation matrix
        proshade_double* rotMatFin                    = ProSHADE_internal_maths::compute3x3MatrixMultiplication ( rotMat2, rotMat );

        //============================================ For each model axis
        std::vector< proshade_double* > hlpAxes;
        for ( proshade_unsign iter = 0; iter < icoAx->getNoAxes ( ); iter++ )
        {
            //======================================== Rotate the model axis to fit the detected orientation
            proshade_double* rotAxis                  = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMatFin,
                                                                                                                        icoAx->getValue ( iter, 1 ),
                                                                                                                        icoAx->getValue ( iter, 2 ),
                                                                                                                        icoAx->getValue ( iter, 3 ) );

            //======================================== Create ProSHADE symmetry axis representation
            proshade_double* axis                     = new proshade_double[7];
            ProSHADE_internal_misc::checkMemoryAllocation ( axis, __FILE__, __LINE__, __func__ );

            axis[0]                                   = icoAx->getValue ( iter, 0 );
            axis[1]                                   = rotAxis[0];
            axis[2]                                   = rotAxis[1];
            axis[3]                                   = rotAxis[2];
            axis[4]                                   = ( 2.0 * M_PI ) / axis[0];
            axis[5]                                   = 0.0;
            axis[6]                                   = -std::numeric_limits < proshade_double >::infinity();

            //======================================== Save axis to ret
            ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( &hlpAxes, axis );

            //======================================== Release memory
            delete[] rotAxis;
            delete[] axis;
        }
        
        //============================================ Save to ret
        ret->emplace_back                             ( hlpAxes );

        //============================================ Release memory
        delete[] rotMat;
        delete[] rotMat2;
        delete[] rotMatFin;
        delete[] rotModelC3;
        delete   icoAx;
    }

    //================================================ Done
    return ;
    
}

/*! \brief This function finds the best pair of axes conforming to the octahedron dihedral angle.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
    \param[in] axErr The error tolerance on angle matching.
    \param[out] ret The pair of axes with closest angle to the required icosahedron dihedral angle.
 */
std::pair< proshade_unsign, proshade_unsign > findBestOctaDihedralPair ( std::vector< proshade_double* >* CSymList, proshade_double minPeakHeight, proshade_double axErr )
{
    //================================================ Initialise variables
    std::pair< proshade_unsign, proshade_unsign > ret;
    std::vector< proshade_unsign  > C4List;
    proshade_double bestHeightSum                     = 0.0;
    proshade_double dotProduct;
    
    //================================================ Find all C5 symmetries
    for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ ) { const FloatingPoint< proshade_double > lhs1 ( CSymList->at(cSym)[0] ), rhs1 ( 4.0 ); if ( lhs1.AlmostEquals ( rhs1 ) && CSymList->at(cSym)[5] >= minPeakHeight ) { ProSHADE_internal_misc::addToUnsignVector ( &C4List, cSym ); } }
    
    //================================================ For each unique pair of C5 and C3
    for ( proshade_unsign c4 = 0; c4 < static_cast<proshade_unsign> ( C4List.size() ); c4++ )
    {
        for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ )
        {
            //======================================== Compare only C3s to the C5List and only with decent average peak height
            const FloatingPoint< proshade_double > lhs1 ( CSymList->at(cSym)[0] ), rhs1 ( 3.0 );
            if ( !lhs1.AlmostEquals ( rhs1 ) ) { continue; }
            if ( CSymList->at(cSym)[5] < minPeakHeight ) { continue; }
            
            //========================================  Check the angle between the C5 and C3 axes
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &CSymList->at(C4List.at(c4))[1],
                                                                                                     &CSymList->at(C4List.at(c4))[2],
                                                                                                     &CSymList->at(C4List.at(c4))[3],
                                                                                                     &CSymList->at(cSym)[1],
                                                                                                     &CSymList->at(cSym)[2],
                                                                                                     &CSymList->at(cSym)[3] );
            
            //======================================== Is the angle approximately the dihedral angle?
            if ( ( ( 1.0 / sqrt ( 3.0 ) ) > ( std::abs( dotProduct ) - axErr ) ) && ( ( 1.0 / sqrt ( 3.0 ) ) < ( std::abs( dotProduct ) + axErr ) ) )
            {
                if ( bestHeightSum < ( CSymList->at(C4List.at(c4))[5] + CSymList->at(cSym)[5] ) )
                {
                    bestHeightSum                     = ( CSymList->at(C4List.at(c4))[5] + CSymList->at(cSym)[5] );
                    ret.first                         = C4List.at(c4);
                    ret.second                        = cSym;
                }
            }
        }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function predicts all octahedral point group symmetry axes from the cyclic point groups list.
 
    This function starts by finding the rotation matrix corresponding to the angle between the predicted C4 axis and the C4 axis found in the pre-computed
    octahedron model available in the ProSHADE_precomputedValues file. It then proceeds to use this rotation matrix to rotate the pre-computed model C3
    axis to now be in the correct orientation to the detected C4 axis.
 
    Next, the function computes the rotation matrix corresponding to rotation along the detected C4 axis about the angle between the rotated pre-computed model
    C3 axis and the detected C3 axis. Finally, when these two rotation matrices are combined, the resulting rotation matrix is the optimal match rotation between the
    pre-computed model and the detected axes positions. This final rotation matrix is then used to rotate the model axes and these rotated model axes are then the
    predicted axes in the structre.
 
    Please note that the peak heights are set to 0.0 for all predicted axes, as they were not detected in the structure, but were only predicted.
 
    \warning This function assumes that the detectIcosahedralSymmetry() function has successfully run (i.e. returned true).
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret The vector containing all the axes forming icosahedral group or empty vector.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
 */
void ProSHADE_internal_symmetry::predictOctaAxes ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, proshade_double minPeakHeight )
{
    //================================================ Create the tetrahedronAxes object
    ProSHADE_internal_precomputedVals::octahedronAxes *octAx = new ProSHADE_internal_precomputedVals::octahedronAxes ( );
    
    //================================================ Find the best axis combination with dihedral angle and correct folds
    std::pair< proshade_unsign, proshade_unsign > initAxes = findBestOctaDihedralPair ( CSymList, minPeakHeight, axErr );
    
    //================================================ Find rotation between the detected C4 and the model C4 axes.
    proshade_double* rotMat                           = ProSHADE_internal_maths::findRotMatMatchingVectors ( octAx->getValue ( 0, 1 ),
                                                                                                             octAx->getValue ( 0, 2 ),
                                                                                                             octAx->getValue ( 0, 3 ),
                                                                                                             CSymList->at(initAxes.first)[1],
                                                                                                             CSymList->at(initAxes.first)[2],
                                                                                                             CSymList->at(initAxes.first)[3] );

    //================================================ Rotate the model C3 to the correct orientation relative to the detected C4 axis.
    proshade_double* rotModelC3                       = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMat,
                                                                                                                        octAx->getValue ( 3, 1 ),
                                                                                                                        octAx->getValue ( 3, 2 ),
                                                                                                                        octAx->getValue ( 3, 3 ) );
    
    //================================================ Find the angle betwen the rotated model C3 and the detected C3 axes along the detected C4 axis.
    proshade_double bestAng = 0.0, curAngDist, bestAngDist = 999.9;
    proshade_double* rotMatHlp                    = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation ( rotMatHlp, __FILE__, __LINE__, __func__ );
    for ( proshade_double ang = 0.0; ang < ( M_PI * 2.0 ); ang += 0.002 )
    {
        //============================================ Compute rotation matrix for this angle value
        ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( rotMatHlp, CSymList->at(initAxes.first)[1], CSymList->at(initAxes.first)[2], CSymList->at(initAxes.first)[3], ang );
        
        //============================================ Rotate the rotated C2 by the matrix
        proshade_double* rotRotModelC3                = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMatHlp,
                                                                                                                        rotModelC3[0],
                                                                                                                        rotModelC3[1],
                                                                                                                        rotModelC3[2] );
        
        //============================================ Find distance
        curAngDist                                    = std::sqrt ( std::pow ( rotRotModelC3[0] - CSymList->at(initAxes.second)[1], 2.0 ) +
                                                                    std::pow ( rotRotModelC3[1] - CSymList->at(initAxes.second)[2], 2.0 ) +
                                                                    std::pow ( rotRotModelC3[2] - CSymList->at(initAxes.second)[3], 2.0 ) );
        
        //============================================ Save best angle
        if ( curAngDist < bestAngDist ) { bestAngDist = curAngDist; bestAng = ang; }
        
        //============================================ Release memory
        delete[] rotRotModelC3;
    }
    
    //============================================ Release memory
    delete[] rotMatHlp;
    
    //================================================ For the rotation matrix along the detected C5 axis with the same anlge as is between the rotated model C3 and the detected C3 axes.
    proshade_double* rotMat2                          = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat2, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( rotMat2, CSymList->at(initAxes.first)[1], CSymList->at(initAxes.first)[2], CSymList->at(initAxes.first)[3], bestAng );
    
    //================================================ Combine the two rotation matrices into a single rotation matrix
    proshade_double* rotMatFin                        = ProSHADE_internal_maths::compute3x3MatrixMultiplication ( rotMat2, rotMat );
    
    //================================================ For each model axis
    for ( proshade_unsign iter = 0; iter < octAx->getNoAxes ( ); iter++ )
    {
        //============================================ Rotate the model axis to fit the detected orientation
        proshade_double* rotAxis                      = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMatFin,
                                                                                                                        octAx->getValue ( iter, 1 ),
                                                                                                                        octAx->getValue ( iter, 2 ),
                                                                                                                        octAx->getValue ( iter, 3 ) );
        
        //============================================ Create ProSHADE symmetry axis representation
        proshade_double* axis                         = new proshade_double[7];
        ProSHADE_internal_misc::checkMemoryAllocation ( axis, __FILE__, __LINE__, __func__ );
        
        axis[0]                                       = octAx->getValue ( iter, 0 );
        axis[1]                                       = rotAxis[0];
        axis[2]                                       = rotAxis[1];
        axis[3]                                       = rotAxis[2];
        axis[4]                                       = ( 2.0 * M_PI ) / axis[0];
        axis[5]                                       = 0.0;
        axis[6]                                       = -std::numeric_limits < proshade_double >::infinity();
        
        //============================================ Save axis to ret
        ProSHADE_internal_misc::addToDblPtrVector     ( ret, axis );
        
        //============================================ Release memory
        delete[] rotAxis;
    }
    
    //================================================ Release memory
    delete[] rotMat;
    delete[] rotMat2;
    delete[] rotMatFin;
    delete[] rotModelC3;
    delete   octAx;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes the list of C symmetries and finds the ten C3 symmetries with correct angles required for full icosahedral symmetry.
 
    This function is specific to detecting the icosahedral symmetry. It should be called once icosahedral symmetry is suspected (by detecting its dihedral
    angles) and it needs to be fully described. This function specifically searches for the ten C3 symmetries which must all be detected in order to fully
    describe icosahedral symmetry. If all ten are found, the ret vector will have these ten axes added to the already present six C5 axes; alternatively,
    the ret array size will not change.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret The vector containing the already detected axes to which newly detected axes (if any) will be added.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
    \param[in] verbose How loud the announcments should be?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
 */
void ProSHADE_internal_symmetry::findIcos10C3s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > prospectiveC3s, retGrp;
    proshade_double dotProd;
    proshade_unsign noClose, noAway;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of ten C3 axes.", messageShift );

    //================================================ For each C3
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
    {
        //============================================ Use only C3s with hight enough average
        if ( CSymList->at(cIt)[0] != 3.0 || CSymList->at(cIt)[0] < minPeakHeight ) { continue; }

        //============================================ Check the C3 has acos ( std::sqrt ( ( 1.0 + 2.0 / std::sqrt ( 5.0 ) ) / 3.0 ) ) to 3 C5s and acos ( 1.0 - ( std::sqrt ( ( 1.0 + 2.0 / std::sqrt ( 5.0 ) ) / 3.0 ) ) ) to the other three C5s
        noClose = 0; noAway = 0;
        for ( proshade_unsign rIt = 0; rIt < 6; rIt++ )
        {
            dotProd                                   = ProSHADE_internal_maths::computeDotProduct ( &ret->at(rIt)[1],
                                                                                                     &ret->at(rIt)[2],
                                                                                                     &ret->at(rIt)[3],
                                                                                                     &CSymList->at(cIt)[1],
                                                                                                     &CSymList->at(cIt)[2],
                                                                                                     &CSymList->at(cIt)[3] );

            if ( ( std::abs ( dotProd ) > (       ( std::sqrt ( ( 1.0 + 2.0 / std::sqrt ( 5.0 ) ) / 3.0 ) ) - axErr ) ) && ( std::abs ( dotProd ) < (       ( std::sqrt ( ( 1.0 + 2.0 / std::sqrt ( 5.0 ) ) / 3.0 ) ) + axErr ) ) ) { noClose += 1; continue; }
            if ( ( std::abs ( dotProd ) > ( 1.0 - ( std::sqrt ( ( 1.0 + 2.0 / std::sqrt ( 5.0 ) ) / 3.0 ) ) - axErr ) ) && ( std::abs ( dotProd ) < ( 1.0 - ( std::sqrt ( ( 1.0 + 2.0 / std::sqrt ( 5.0 ) ) / 3.0 ) ) + axErr ) ) ) { noAway  += 1; continue; }
        }

        //============================================ If correct angles distribution is found, save the axis
        if ( ( noClose == 3 ) && ( noAway == 3 ) )
        {
            ProSHADE_internal_misc::addToUnsignVector ( &prospectiveC3s, cIt );
        }
    }

    //================================================ Search for missing axes
    for ( proshade_unsign iter = 0; iter < 6; iter++ ) { ProSHADE_internal_misc::addToUnsignVector ( &retGrp, iter ); }
    if ( !ProSHADE_internal_symmetry::findMissingAxesDual ( &prospectiveC3s, CSymList, ret, &retGrp, 10, axErr, 3, std::sqrt ( ( 1.0 + 2.0 / std::sqrt ( 5.0 ) ) / 3.0 ), 3, 1.0 - ( std::sqrt ( ( 1.0 + 2.0 / std::sqrt ( 5.0 ) ) / 3.0 ) ), 3, dataObj ) )
    {
        return ;
    }

    //================================================ Found correct number of axes! Now save the
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( prospectiveC3s.size() ); iter++ )
    {
        ProSHADE_internal_misc::addToDblPtrVector     ( ret, CSymList->at(prospectiveC3s.at(iter)) );
    }

    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "Detection of ten C3 axes successfull.", messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes the list of C symmetries and finds the fifteen C3 symmetries with correct angles required for full icosahedral symmetry.
 
    This function is specific to detecting the icosahedral symmetry. It should be called once icosahedral symmetry is suspected (by detecting its dihedral
    angles) and it needs to be fully described. This function specifically searches for the ten C3 symmetries which must all be detected in order to fully
    describe icosahedral symmetry. If all ten are found, the ret vector will have these ten axes added to the already present six C5 axes; alternatively,
    the ret array size will not change.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret The vector containing the already detected axes to which newly detected axes (if any) will be added.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
    \param[in] verbose How loud the announcments should be?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
 */
void ProSHADE_internal_symmetry::findIcos15C2s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > prospectiveC2s, retGrp;
    proshade_double dotProd;
    proshade_unsign noClose, noMidway, noAway;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of fifteen C2 axes.", messageShift );

    //================================================ For each C2
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
    {
        //============================================ Use only C2s
        const FloatingPoint< proshade_double > lhs999 ( CSymList->at(cIt)[5] ), rhs999 ( static_cast< proshade_double > ( -999.9 ) );
        if ( CSymList->at(cIt)[0] != 2.0 || ( ( CSymList->at(cIt)[5] < minPeakHeight ) && !( lhs999.AlmostEquals( rhs999 ) ) ) ) { continue; }

        //============================================ Check the C2 has acos ( 0.0 ) to 2 C5s, acos ( 0.5 ) to another 2 C5s and acos ( sqrt ( 3.0 ) / 2.0 ) to the last two C5s
        noClose = 0; noMidway = 0; noAway = 0;
        for ( proshade_unsign rIt = 0; rIt < 6; rIt++ )
        {
            dotProd                                   = ProSHADE_internal_maths::computeDotProduct ( &ret->at(rIt)[1],
                                                                                                     &ret->at(rIt)[2],
                                                                                                     &ret->at(rIt)[3],
                                                                                                     &CSymList->at(cIt)[1],
                                                                                                     &CSymList->at(cIt)[2],
                                                                                                     &CSymList->at(cIt)[3] );

            if ( ( std::abs ( dotProd ) > ( ( sqrt ( 3.0 ) / 2.0 ) - axErr ) ) && ( std::abs ( dotProd ) < ( ( sqrt ( 3.0 ) / 2.0 ) + axErr ) ) ) { noAway    += 1; continue; }
            if ( ( std::abs ( dotProd ) > ( ( 1.0 / 2.0 )          - axErr ) ) && ( std::abs ( dotProd ) < ( ( 1.0 / 2.0 )          + axErr ) ) ) { noMidway  += 1; continue; }
            if ( ( std::abs ( dotProd ) > ( ( 0.0 )                - axErr ) ) && ( std::abs ( dotProd ) < ( ( 0.0 )                + axErr ) ) ) { noClose   += 1; continue; }
        }

        //============================================ If correct angles distribution is found, save the axis
        if ( ( noClose == 2 ) && ( noMidway == 2 ) && ( noAway == 2 ) )
        {
            ProSHADE_internal_misc::addToUnsignVector ( &prospectiveC2s, cIt );
        }
    }

    //================================================ Search for missing axes
    for ( proshade_unsign iter = 0; iter < 6; iter++ ) { ProSHADE_internal_misc::addToUnsignVector ( &retGrp, iter ); }
    if ( !ProSHADE_internal_symmetry::findMissingAxesTriple ( &prospectiveC2s, CSymList, ret, &retGrp, 15, axErr, 2, 0.0, 2, 1.0/2.0, 2, sqrt ( 3.0 ) / 2.0, 2, dataObj ) )
    {
        return ;
    }
    
    //================================================ Found correct number of axes! Now save the
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( prospectiveC2s.size() ); iter++ )
    {
        ProSHADE_internal_misc::addToDblPtrVector     ( ret, CSymList->at(prospectiveC2s.at(iter)) );
    }

    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "Detection of fifteen C2 axes successfull.", messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function tries to find a particular symmetry axis which would complete a group of symmetries with three different angle requirement to another group.
 
    Assuming there is a group of symmetry axis, which have particular number of particular angles to each other, but some are missing, this function tries to find any such
    missing axes. This is a solution for the group of axes having three different angles to the other group members. For all newly detected group members, the average peak
    height and the uniqueness are both tested for.
 
    \param[in] possibilities A vector of already detected axis indices which should be extended.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret A list of already detected octahedral axes.
    \param[in] retGroup A vector of indices in the ret list which form the group to which new axes are compared to.
    \param[in] requiredNoAxes Number of axes required for positive result.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] noMatchesG1 The number of axes from ret that need to be matched with angle1.
    \param[in] angle1 The angle with which noMatchesG1 axes need to be matched with the retGroup axes.
    \param[in] noMatchesG2 The number of axes from ret that need to be matched with angle2.
    \param[in] angle2 The angle with which noMatchesG2 axes need to be matched with the retGroup axes.
    \param[in] noMatchesG3 The number of axes from ret that need to be matched with angle3.
    \param[in] angle3 The angle with which noMatchesG3 axes need to be matched with the retGroup axes.
    \param[in] fold The fold of the searched for axis.
    \param[in] dataObj The full data holding object pointer - this is to get access to self-rotation function values.
    \param[out] atLeastOne Boolean value speciying whether at least the minimum required number of axes was found.
 */
bool ProSHADE_internal_symmetry::findMissingAxesTriple ( std::vector< proshade_unsign >* possibilities, std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, std::vector< proshade_unsign >* retGroup, proshade_unsign requiredNoAxes, proshade_double axErr, proshade_unsign noMatchesG1, proshade_double angle1, proshade_unsign noMatchesG2, proshade_double angle2, proshade_unsign noMatchesG3, proshade_double angle3, proshade_unsign fold, ProSHADE_internal_data::ProSHADE_data* dataObj )
{
    //================================================ Initialise variables
    bool atLeastOne                                   = false;
    std::vector< proshade_double* > prosp;
    std::vector< proshade_double > sol;
    
    //================================================ Proceed only if need be
    if ( static_cast<proshade_unsign> ( possibilities->size() ) == requiredNoAxes ) { atLeastOne = true; return ( atLeastOne ); }
    
    //================================================ Copy already found to prospective
    for ( proshade_unsign prIt = 0; prIt < static_cast<proshade_unsign> ( possibilities->size() ); prIt++ )
    {
        ProSHADE_internal_symmetry::addAxisUnlessSame ( static_cast< proshade_unsign > ( CSymList->at(possibilities->at(prIt))[0] ),
                                                        CSymList->at(possibilities->at(prIt))[1],
                                                        CSymList->at(possibilities->at(prIt))[2],
                                                        CSymList->at(possibilities->at(prIt))[3],
                                                        CSymList->at(possibilities->at(prIt))[5], &prosp, axErr );
    }
    
    //================================================ Start generating possible solutions
    for ( proshade_unsign rgIt1 = 0; rgIt1 < static_cast<proshade_unsign> ( retGroup->size() ); rgIt1++ )
    {
        for ( proshade_unsign rgIt2 = 0; rgIt2 < static_cast<proshade_unsign> ( retGroup->size() ); rgIt2++ )
        {
            //======================================== Use unique combinations (order matters here!)
            if ( rgIt1 == rgIt2 ) { continue; }
            
            for ( proshade_unsign rgIt3 = 0; rgIt3 < static_cast<proshade_unsign> ( retGroup->size() ); rgIt3++ )
            {
                //==================================== Use unique combinations (order matters here!)
                if ( ( rgIt1 == rgIt3 ) || ( rgIt2 == rgIt3 ) ) { continue; }
                
                //==================================== Generate possible solution (1)
                sol                                   = ProSHADE_internal_maths::findVectorFromThreeVAndThreeD ( ret->at(rgIt1)[1], ret->at(rgIt1)[2], ret->at(rgIt1)[3],
                                                                                                                 ret->at(rgIt2)[1], ret->at(rgIt2)[2], ret->at(rgIt2)[3],
                                                                                                                 ret->at(rgIt3)[1], ret->at(rgIt3)[2], ret->at(rgIt3)[3], angle1, angle2, angle3 );
                
                //==================================== Check if solution fits the group completely
                ProSHADE_internal_symmetry::checkFittingAxisTripleAndSave ( retGroup, ret, fold, sol.at(0), sol.at(1), sol.at(2), &prosp, axErr, noMatchesG1, angle1, noMatchesG2, angle2, noMatchesG3, angle3, dataObj );
                if ( prosp.size() == requiredNoAxes ) { break; }
                
                //==================================== Generate possible solution (2)
                sol                                   = ProSHADE_internal_maths::findVectorFromThreeVAndThreeD ( ret->at(rgIt1)[1], ret->at(rgIt1)[2], ret->at(rgIt1)[3],
                                                                                                                 ret->at(rgIt2)[1], ret->at(rgIt2)[2], ret->at(rgIt2)[3],
                                                                                                                 ret->at(rgIt3)[1], ret->at(rgIt3)[2], ret->at(rgIt3)[3], -angle1, -angle2, -angle3 );
                
                //==================================== Check if solution fits the group completely
                ProSHADE_internal_symmetry::checkFittingAxisTripleAndSave ( retGroup, ret, fold, sol.at(0), sol.at(1), sol.at(2), &prosp, axErr, noMatchesG1, angle1, noMatchesG2, angle2, noMatchesG3, angle3, dataObj );
                if ( prosp.size() == requiredNoAxes ) { break; }
            }
            
            if ( prosp.size() == requiredNoAxes ) { break; }
        }
        
        if ( prosp.size() == requiredNoAxes ) { break; }
    }
    
    //================================================ Found all required axes
    if ( prosp.size() == requiredNoAxes )
    {
        //============================================ For each found missing axis
        for ( proshade_unsign axIt = static_cast<proshade_unsign> ( possibilities->size() ); axIt < static_cast<proshade_unsign> ( prosp.size() ); axIt++ )
        {
            if ( ProSHADE_internal_maths::isAxisUnique ( CSymList, prosp.at(axIt), axErr ) )
            {
                //======================================== Add
                ProSHADE_internal_misc::addToDblPtrVector ( CSymList, prosp.at(axIt) );
                ProSHADE_internal_misc::addToUnsignVector ( possibilities, static_cast<proshade_unsign> ( CSymList->size()-1 ) );
            }
        }
        
        atLeastOne                                    = true;
        return                                        ( atLeastOne );
    }
    else
    {
        //============================================ Delete all found, but unnecessary axes
        for ( proshade_unsign axIt = static_cast<proshade_unsign> ( possibilities->size() ); axIt < static_cast<proshade_unsign> ( prosp.size() ); axIt++ )
        {
            delete[] prosp.at(axIt);
        }
    }
    
    //================================================ Done
    return                                            ( atLeastOne );
    
}

/*! \brief This function takes a newly detected "missing" axis and tests it for belonging to the group, checking the height and replacing lower height members with better members.
 
    This function takes the list of already detected axes, information about the tested new axis and the conditions for belonging. It then proceeds to check if the
    new axis conforms to the conditions of belonging. If so, it then checks if the axis height is high enough to be considered as part of the group. Again, if so,
    it will save this new axis to the old set, replacing any old axis with this new one, if it is the same and has better height.
 
    \param[in] retGroup A vector of indices in the ret list which form the group to which new axes are compared to.
    \param[in] ret A list of already detected axes.
    \param[in] fold The fold of the searched for axis.
    \param[in] axX The x-axis element of the new axis.
    \param[in] axY The y-axis element of the new axis.
    \param[in] axZ The z-axis element of the new axis.
    \param[in] prosp The vector to which the axis is to be saved.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] noMatchesG1 The number of axes from ret that need to be matched with angle1.
    \param[in] angle1 The angle with which noMatchesG1 axes need to be matched with the retGroup axes.
    \param[in] noMatchesG2 The number of axes from ret that need to be matched with angle2.
    \param[in] angle2 The angle with which noMatchesG2 axes need to be matched with the retGroup axes.
    \param[in] noMatchesG3 The number of axes from ret that need to be matched with angle3.
    \param[in] angle3 The angle with which noMatchesG3 axes need to be matched with the retGroup axes.
    \param[in] dataObj The full data holding object pointer - this is to get access to self-rotation function values.
 */
void ProSHADE_internal_symmetry::checkFittingAxisTripleAndSave ( std::vector< proshade_unsign >* retGroup, std::vector< proshade_double* >* ret, proshade_unsign fold, proshade_double axX, proshade_double axY, proshade_double axZ, std::vector< proshade_double* >* prosp, proshade_double axErr, proshade_unsign noMatchesG1, proshade_double angle1, proshade_unsign noMatchesG2, proshade_double angle2, proshade_unsign noMatchesG3, proshade_double angle3, ProSHADE_internal_data::ProSHADE_data* dataObj )
{
    //================================================ Initialise variables
    proshade_unsign noG1                              = 0;
    proshade_unsign noG2                              = 0;
    proshade_unsign noG3                              = 0;
    proshade_double dotProd                           = 0.0;
    proshade_double axHeight                          = 0.0;
    
    //================================================ Find the angle and count dual matching frequencies
    for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( retGroup->size() ); rIt++ )
    {
        dotProd                                       = ProSHADE_internal_maths::computeDotProduct ( &ret->at(retGroup->at(rIt))[1],
                                                                                                     &ret->at(retGroup->at(rIt))[2],
                                                                                                     &ret->at(retGroup->at(rIt))[3],
                                                                                                     &axX, &axY, &axZ );
        
        if ( ( std::abs ( dotProd ) > ( angle1 - axErr ) ) && ( std::abs ( dotProd ) < ( angle1 + axErr ) ) ) { noG1 += 1; continue; }
        if ( ( std::abs ( dotProd ) > ( angle2 - axErr ) ) && ( std::abs ( dotProd ) < ( angle2 + axErr ) ) ) { noG2 += 1; continue; }
        if ( ( std::abs ( dotProd ) > ( angle3 - axErr ) ) && ( std::abs ( dotProd ) < ( angle3 + axErr ) ) ) { noG3 += 1; continue; }
    }
    
    //================================================ If correct frequencies are matched, check height.
    if ( ( noG1 == noMatchesG1 ) && ( noG2 == noMatchesG2 ) && ( noG3 == noMatchesG3 ) )
    {
        //============================================ Is the height good enough?
        axHeight                                      = ProSHADE_internal_symmetry::missingAxisHeight ( axX, axY, axZ, dataObj, fold, axErr );
        
        //============================================ If so, save
        if ( axHeight > 0.1 )
        {
            ProSHADE_internal_symmetry::addAxisUnlessSame ( fold, axX, axY, axZ, axHeight, prosp, axErr );
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function obtains a list of all C symmetries from the angle-axis space mapped rotation function values.

    This function oversees the full search for cyclic point groups in the self-rotation function. It starts with finding all prime numbers up to the
    user specified limit. It then checks for each of the prime numbers whether there is a cyclic point group with fold equal to the prime number.
 
    If any such point groups are found, the function searches for nultiples of these folds, making use of the fact that any structure with cyclic
    point group of fold n must also contain a point group of fold n/2 if n/2 is an integer. In this manner, cyclic point group with any fold can be
    found using a small number of specific fold searches, thus eliminating the need to determine which folds should be considered.

    \param[in] settings A pointer to settings class containing all the information required for symmetry detection.
    \param[out] ret A vector of arrays containing all detected cyclic point groups in the standard ProSHADE format, i.e. [0] = fold, [1] = x-axis, [2] = y-axis, [3] = z-axis, [4] = angle, [5] = average peak height.
 */
std::vector< proshade_double* > ProSHADE_internal_data::ProSHADE_data::getCyclicSymmetriesListFromAngleAxis ( ProSHADE_settings* settings )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign  > primes            = ProSHADE_internal_maths::findAllPrimes ( settings->maxSymmetryFold );
    std::vector< proshade_double* > ret, tmpHolder;
    std::vector< proshade_unsign  > testedFolds;
    proshade_double symThres;
    proshade_unsign foldToTest;
    bool foldDone, anyNewSyms = true;
    
    //================================================ For each found prime number in the limit
    for ( proshade_unsign prIt = 0; prIt < static_cast< proshade_unsign > ( primes.size() ); prIt++ )
    {
        //============================================ Report progress
        std::stringstream hlpSS;
        hlpSS << "Searching for prime fold symmetry C" << primes.at(prIt) << ".";
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 3, hlpSS.str(), settings->messageShift );
        
        //============================================ Get all symmetries for this prime fold
        std::vector< proshade_double* > prSyms        = this->findRequestedCSymmetryFromAngleAxis ( settings, primes.at(prIt), &symThres );
        
        //============================================ Save the detected C symmetries
        for ( size_t axIt = 0; axIt < prSyms.size(); axIt++ )
        {
            //======================================== Is this symmetry passing the threshold?
            if ( prSyms.at(axIt)[5] >= symThres )
            {
                //==================================== Add this symmetry to final list
                if ( ProSHADE_internal_maths::isAxisUnique ( &ret, prSyms.at(axIt), settings->axisErrTolerance, true ) )
                {
                    ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( &ret, prSyms.at(axIt) );
                }
            }
        
            //======================================== Release memory
            delete[] prSyms.at(axIt);
        }
    }
    
    //================================================ Was anything found?
    if ( ret.size() < 1 ) { return ( ret ); }
    
    //================================================ Check for prime symmetry fold multiples
    while ( anyNewSyms )
    {
        //============================================ Initialise new iteration
        anyNewSyms                                    = false;
        
        //============================================ For each passing symmetry, look if there are any combinations of symmetries that would contain it
        for ( proshade_unsign axIt1 = 0; axIt1 < static_cast< proshade_unsign > ( ret.size() ); axIt1++ )
        {
            for ( proshade_unsign axIt2 = 0; axIt2 < static_cast< proshade_unsign > ( ret.size() ); axIt2++ )
            {
                //==================================== Initialise iteration
                foldToTest                            = static_cast< proshade_unsign > ( ret.at(axIt1)[0] * ret.at(axIt2)[0] );
                if ( foldToTest > settings->maxSymmetryFold ) { continue; }
                
                //==================================== Was this fold tested already?
                foldDone                              = false;
                for ( proshade_unsign fIt = 0; fIt < static_cast< proshade_unsign > ( testedFolds.size() ); fIt++ ) { if ( testedFolds.at(fIt) == foldToTest ) { foldDone = true; break; } }
                if ( foldDone ) { continue; }
                else            { ProSHADE_internal_misc::addToUnsignVector ( &testedFolds, foldToTest ); }
                
                //==================================== Report progress
                std::stringstream hlpSS2;
                hlpSS2 << "Searching for fold combination of detected folds " << ret.at(axIt1)[0] << " and " << ret.at(axIt2)[0] << ", i.e. C" << foldToTest << ".";
                ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 3, hlpSS2.str(), settings->messageShift );
                
                //==================================== Get all symmetries for this fold
                std::vector< proshade_double* > prSyms = this->findRequestedCSymmetryFromAngleAxis ( settings, foldToTest, &symThres );
                
                //==================================== For each detected group with the required fold
                for ( size_t newAxIt = 0; newAxIt < prSyms.size(); newAxIt++ )
                {
                    if ( prSyms.at(newAxIt)[5] >= symThres )
                    {
                        //================================ Add to detected axes
                        if ( ProSHADE_internal_maths::isAxisUnique ( &ret, prSyms.at(newAxIt), settings->axisErrTolerance, true ) )
                        {
                            ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( &tmpHolder, prSyms.at(newAxIt) );
                        }
                    }

                    //==================================== Release memory
                    delete[] prSyms.at(newAxIt);
                }
            }
        }
        
        //============================================ Add newly found groups and repeat if need be
        if ( tmpHolder.size() > 0 )
        {
            for ( proshade_unsign tmpIt = 0; tmpIt < static_cast< proshade_unsign > ( tmpHolder.size() ); tmpIt++ )
            {
                ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( &ret, tmpHolder.at(tmpIt) );
                delete[] tmpHolder.at(tmpIt);
            }
            
            anyNewSyms                                = true;
            tmpHolder.clear                           ( );
        }
    }
    
    //================================================ Sort the vector
    std::sort                                         ( ret.begin(), ret.end(), ProSHADE_internal_misc::sortSymHlpInv );
    
    //================================================ Done
    return                                            ( ret );
}

/*! \brief This function allows using std::sort to sort vectors of ProSHADE symmetry format..
 
    \param[in] a Pointer to a ProSHADE symmetry formatted array.
    \param[in] b Pointer to a ProSHADE symmetry formatted array.
    \param[out] X  Boolean whether a is larger than b.
 */
bool sortProSHADESymmetryByPeak ( proshade_double* a, proshade_double* b)
{
    //================================================ Done
    return ( a[5] > b[5] );

}

/*! \brief This function searches the angle-axis representation of the rotation function for a cyclic point group with given fold.
 
    This function is a simplification of the getCyclicSymmetriesListFromAngleAxis() function, where this function does not map the whole
    rotation function, but rothar only to the spheres it knows it will required. Moreover, it does not search for all cyclic point groups, but only
    those which span all the spheres (angles) and therefore have the required fold.
 
    In terms of operations, this function interpolates the rotation function values onto the spheres it requires, then it finds peaks and removes
    the small peaks, so that these can then be grouped. For each group which spans all the angles it then finds the index with highest sum of
    peak height over all spheres. It can then do the bi-cubic interpolation if requested. Finally, all the detected peaks are sorted by the peak
    height and returned.
 
    \param[in] settings ProSHADE_settings object containing all the settings for this run.
    \param[in] fold The fold which should be sought for by the function.
    \param[in] peakThres The threshold used to cut peaks.
    \param[out] ret  Vector of double pointers to arrays having the standard ProSHADE symmetry group structure.
 */
std::vector < proshade_double* > ProSHADE_internal_data::ProSHADE_data::findRequestedCSymmetryFromAngleAxis ( ProSHADE_settings* settings, proshade_unsign fold, proshade_double* peakThres )
{    
    //================================================ Initialise variables
    proshade_double soughtAngle;
    std::vector< proshade_double  > allPeakHeights;
    std::vector< ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup* > peakGroups;
    std::vector< proshade_double* > ret;
    bool newPeak;
    
    //================================================ Make sure we have a clean start
    this->sphereMappedRotFun.clear();
    
    //================================================ Convert rotation function to only the required angle-axis space spheres and find all peaks
    for ( proshade_double angIt = 1.0; angIt < static_cast < proshade_double > ( fold ); angIt += 1.0 )
    {
        //============================================ Figure the angles to form the symmetry
        soughtAngle                                   = angIt * ( 2.0 * M_PI / static_cast<proshade_double> ( fold ) );
        
        //============================================ Create the angle-axis sphere with correct radius (angle)
        this->sphereMappedRotFun.emplace_back         ( new ProSHADE_internal_spheres::ProSHADE_rotFun_sphere ( soughtAngle,
                                                                                                                M_PI / static_cast < proshade_double > ( this->maxShellBand ),
                                                                                                                this->maxShellBand * 2,
                                                                                                                soughtAngle,
                                                                                                                static_cast < proshade_unsign > ( angIt - 1.0 ) ) );
        
        //=========================================== Interpolate rotation function onto the sphere
        this->sphereMappedRotFun.at(static_cast < size_t > ( angIt - 1.0 ))->interpolateSphereValues ( this->getInvSO3Coeffs ( ) );
        
        //============================================ Find all peaks for this sphere
        this->sphereMappedRotFun.at(static_cast < size_t > ( angIt - 1.0 ))->findAllPeaks ( static_cast< proshade_signed > ( settings->peakNeighbours ), &allPeakHeights );
    }
    
    //============================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "Found a total of " << std::pow ( static_cast< proshade_double > ( this->maxShellBand ) * 2.0 * ( static_cast< proshade_double > ( fold ) - 1.0 ), 2.0 ) - static_cast< proshade_double > ( allPeakHeights.size() ) << " non-peaks for thresholding.";
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 4, hlpSS.str(), settings->messageShift );
    
    //================================================ Determine the threshold for significant peaks
   *peakThres                                         = std::min ( settings->minSymPeak, determinePeakThreshold ( allPeakHeights, settings->noIQRsFromMedianNaivePeak ) );
    
    //============================================ Report progress
    std::stringstream hlpSS2;
    hlpSS2 << "Determined peak threshold " << *peakThres << ".";
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 4, hlpSS2.str(), settings->messageShift );
    
    //================================================ Remove small peaks
    for ( proshade_unsign shIt = 0; shIt < static_cast<proshade_unsign> ( this->sphereMappedRotFun.size() ); shIt++ )
    {
        this->sphereMappedRotFun.at(shIt)->removeSmallPeaks ( *peakThres );
    }
    
    //================================================ Group peaks
    for ( proshade_unsign sphIt = 0; sphIt < static_cast<proshade_unsign> ( this->sphereMappedRotFun.size() ); sphIt++ )
    {
        //============================================ For each peak
        for ( proshade_unsign pkIt = 0; pkIt < static_cast<proshade_unsign> ( this->sphereMappedRotFun.at(sphIt)->getPeaks().size() ); pkIt++ )
        {
            //======================================== Check if peak belongs to an already detected peak group
            newPeak                                   = true;
            for ( proshade_unsign pkGrpIt = 0; pkGrpIt < static_cast<proshade_unsign> ( peakGroups.size() ); pkGrpIt++ )
            {
                if ( peakGroups.at(pkGrpIt)->checkIfPeakBelongs ( static_cast< proshade_double > ( this->sphereMappedRotFun.at(sphIt)->getPeaks().at(pkIt).first  ),
                                                                  static_cast< proshade_double > ( this->sphereMappedRotFun.at(sphIt)->getPeaks().at(pkIt).second ),
                                                                  sphIt, settings->axisErrTolerance, settings->verbose, settings->messageShift ) ) { newPeak = false; break; }
            }
            
            //======================================== If already added, go to next one
            if ( !newPeak ) { continue; }
            
            //======================================== If not, create a new group with this peak
            peakGroups.emplace_back                   ( new ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup ( static_cast< proshade_double > ( this->sphereMappedRotFun.at(sphIt)->getPeaks().at(pkIt).first ),
                                                                                                                         static_cast< proshade_double > ( this->sphereMappedRotFun.at(sphIt)->getPeaks().at(pkIt).second ),
                                                                                                                         sphIt,
                                                                                                                         this->sphereMappedRotFun.at(sphIt)->getAngularDim() ) );
        }
    }
    
    //================================================ For each peak group, look for the requested fold
    for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( peakGroups.size() ); grIt++ )
    {
        //============================================ Report progress
        std::stringstream hlpSS3;
        hlpSS3 << "Now considering group with LAT " << peakGroups.at(grIt)->getLatFromIndices() << " - " << peakGroups.at(grIt)->getLatToIndices() << " and LON " << peakGroups.at(grIt)->getLonFromIndices() << " - " << peakGroups.at(grIt)->getLonToIndices() << " spanning spheres ";
        for ( proshade_unsign sphIt = 0; sphIt < static_cast<proshade_unsign> ( peakGroups.at(grIt)->getSpherePositions().size() ); sphIt++ ) { hlpSS3 << peakGroups.at(grIt)->getSpherePositions().at(sphIt) << " ; "; }
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 5, hlpSS3.str(), settings->messageShift );
        
        //============================================ Find point groups in the peak group
        peakGroups.at(grIt)->findCyclicPointGroupsGivenFold ( this->sphereMappedRotFun, &ret, settings->useBiCubicInterpolationOnPeaks, fold, settings->verbose, settings->messageShift );

        //============================================ Release the memory
        delete peakGroups.at(grIt);
    }
    
    //================================================ Sort ret by peak height
    std::sort                                         ( ret.begin(), ret.end(), sortProSHADESymmetryByPeak );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function finds the rotation function value for all axes supplied in the ret parameter.
 
    This function supplements the polyhedral symmetry prediction functions, as these functions predict the symmetry axes, but do not
    find their peak heights. This function, then, firstly finds all the individual folds in the symmetry axes set and for each fold computes
    the appropriate angles. Next. it computes the sphere mappings of the rotation function for all detected angles and for each symmetry
    axes, it finds the rotation function average as well as the average for the whole symmetry (i.e. over all axes). Finally, the function attempts
    to locally optimise the detected symmetry group by searching for slightly rotated axes having higher rotation function average.
 
    \param[in] ret The list of axes for which the heights are to be found.
    \param[in] dataObj The structure object with computed rotation function in which the peaks are to be found.
    \param[in] settings ProSHADE_settings object containing all the settings for this run.
    \param[in] improve Should the axes position be optimised, or is prediction sufficient?
 */
void ProSHADE_internal_symmetry::findPredictedAxesHeights ( std::vector< proshade_double* >* ret, ProSHADE_internal_data::ProSHADE_data* dataObj, ProSHADE_settings* settings, bool improve )
{
    //================================================ Initialise variables
    std::vector < proshade_unsign > folds;
    std::vector < proshade_double > angs, applicableAngs;
    bool alreadyFound = false;
    size_t corAngIt = 0;
    proshade_double lat = 0.0, lon = 0.0, radRange = 0.0, searchRangeInDeg = 0.0, axSum = 0.0, curSum = 0.0, maxSum = 0.0, bestXRot = 0.0, bestYRot = 0.0, bestZRot = 0.0, finXRotChan = 0.0, finYRotChan = 0.0, finZRotChan = 0.0;
    proshade_double latSamlUnit                       = ( 2.0 * M_PI ) / ( static_cast< proshade_double > ( dataObj->maxShellBand ) * 2.0 );
    proshade_double lonSamlUnit                       = ( 1.0 * M_PI ) / ( static_cast< proshade_double > ( dataObj->maxShellBand ) * 2.0 );
    
    //================================================ Determine all the folds for which rotation function mapping will be required
    for ( proshade_unsign iter = 0; iter < static_cast < proshade_unsign > ( ret->size() ); iter++ )
    {
        alreadyFound                                  = false;
        for ( proshade_unsign it = 0; it < static_cast < proshade_unsign > ( folds.size() ); it++ ) { const FloatingPoint< proshade_double > lhs1 ( static_cast< proshade_double > ( folds.at(it) ) ), rhs1 ( ret->at(iter)[0] ); if ( lhs1.AlmostEquals ( rhs1 ) ) { alreadyFound = true; break; } }
        
        if ( !alreadyFound ) { ProSHADE_internal_misc::addToUnsignVector ( &folds, static_cast< proshade_unsign > ( ret->at(iter)[0] ) ); }
    }
    
    //================================================ Generate vector of all angles
    for ( proshade_unsign foldIt = 0; foldIt < static_cast < proshade_unsign > ( folds.size() ); foldIt++ ) { for ( proshade_double angIt = 1.0; angIt < static_cast<proshade_double> ( folds.at(foldIt) ); angIt += 1.0 ) { ProSHADE_internal_misc::addToDoubleVector ( &angs, angIt * ( 2.0 * M_PI / static_cast<proshade_double> ( folds.at(foldIt) ) ) ); } }
    std::sort                                         ( angs.begin(), angs.end() );
    
    //================================================ Remove redundant angles from the list
    for ( int angIt = static_cast< int > ( angs.size() - 2 ); angIt >= 0; angIt-- ) { const FloatingPoint< proshade_double > lhs1 ( angs.at(static_cast< size_t > ( angIt ) ) ), rhs1 ( angs.at(static_cast< size_t > ( angIt + 1 ) ) ); if ( lhs1.AlmostEquals ( rhs1 ) ) { angs.erase ( angs.begin() + (angIt+1) ); } }
    
    //================================================ Generate all sphere mapped rotation function
    dataObj->sphereMappedRotFun.clear();
    for ( size_t angIt = 0; angIt < angs.size(); angIt++ )
    {
        //============================================ Decide the range in which the sphere operates
        if ( ( angIt == 0 ) && ( angs.size() > 1 ) )  { radRange = ( angs.at(1) - angs.at(0) ) / 2; }
        else { if ( ( angIt == ( angs.size() - 1 ) ) && ( angs.size() > 1 ) ) { radRange = ( angs.at(angIt) - angs.at(angIt-1) ) / 2; }
               else { if ( angs.size() > 2 ) { radRange = std::min ( ( angs.at(angIt) - angs.at(angIt-1) ) / 2, ( angs.at(angIt+1) - angs.at(angIt) ) / 2 ); }
                      else { radRange = 0.5; } } }
        
        //============================================ Create the sphere
        dataObj->sphereMappedRotFun.emplace_back      ( new ProSHADE_internal_spheres::ProSHADE_rotFun_sphere ( angs.at(angIt),
                                                                                                                radRange,
                                                                                                                dataObj->maxShellBand * 2,
                                                                                                                angs.at(angIt),
                                                                                                                static_cast<proshade_unsign> ( angIt ) ) );

        //=========================================== Interpolate rotation function onto the sphere
        dataObj->sphereMappedRotFun.at( static_cast < proshade_unsign > ( angIt ))->interpolateSphereValues ( dataObj->getInvSO3Coeffs ( ) );
    }
    
    //================================================ Check for improved sum
    searchRangeInDeg                                  = 360.0 / ( dataObj->getMaxBand() * 2.0 );
    proshade_double* rotMat, *newAxis;
    while ( searchRangeInDeg > 0.09 )
    {
        //============================================ For change along each dimension
        for ( proshade_double xChan = -searchRangeInDeg; xChan < ( 1.5 * searchRangeInDeg ); xChan += searchRangeInDeg )
        {
            for ( proshade_double yChan = -searchRangeInDeg; yChan < ( 1.5 * searchRangeInDeg ); yChan += searchRangeInDeg )
            {
                for ( proshade_double zChan = -searchRangeInDeg; zChan < ( 1.5 * searchRangeInDeg ); zChan += searchRangeInDeg )
                {
                    //================================ Initialise local variables
                    curSum                            = 0.0;
                    
                    //================================ Find the rotation matrix of the appropriate rotation
                    rotMat                            = ProSHADE_internal_maths::build3x3MatrixFromXYZRotations ( xChan + finXRotChan, yChan + finYRotChan, zChan + finZRotChan );
                    
                    //================================ For each axis, find new position and its RF value
                    for ( proshade_unsign axIt = 0; axIt < static_cast< proshade_unsign > ( ret->size() ); axIt++ )
                    {
                        //============================ Find rotated axis
                        newAxis                       = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMat, ret->at(axIt)[1], ret->at(axIt)[2], ret->at(axIt)[3] );
                        
                        //============================ Convert XYZ to lat and lon INDICES
                        lat                           = ( std::atan2( newAxis[1], newAxis[0] ) / latSamlUnit );
                        lon                           = ( std::acos ( newAxis[2] ) / lonSamlUnit );

                        if ( lat < 0.0 )                  { lat += ( static_cast< proshade_double > ( dataObj->maxShellBand ) * 2.0 ); }
                        if ( lon < 0.0 )                  { lon += ( static_cast< proshade_double > ( dataObj->maxShellBand ) * 2.0 ); }

                        //============================ Generate all angles for this fold
                        applicableAngs.clear          ( );
                        for ( proshade_double angIt = 1.0; angIt < ret->at(axIt)[0]; angIt += 1.0 ) { ProSHADE_internal_misc::addToDoubleVector ( &applicableAngs, angIt * ( 2.0 * M_PI / ret->at(axIt)[0] ) ); }

                        //============================ For each shpere with the correct angle, average the peak heights
                        axSum                         = 1.0;
                        for ( size_t angIt = 0; angIt < angs.size(); angIt++ )
                        {
                            //======================== Find the correct sphere
                            alreadyFound              = false;
                            for ( size_t aIt = 0; aIt < applicableAngs.size(); aIt++ ) { if ( alreadyFound ) { break; } const FloatingPoint< proshade_double > lhs1 ( angs.at(angIt) ), rhs1 ( applicableAngs.at(aIt) ); if ( lhs1.AlmostEquals ( rhs1 ) ) { alreadyFound = true; corAngIt = angIt; } }

                            if ( !alreadyFound ) { continue; }

                            //======================== Get its peak height for the longitude and latitude
                            axSum                    += dataObj->sphereMappedRotFun.at(corAngIt)->getSphereLatLonLinearInterpolationPos ( lat, lon );
                        }

                        //============================ And average the peak heights over the axis
                        axSum                        /= ret->at(axIt)[0];
                        curSum                       += axSum;
                        
                        //============================ Release memory
                        delete[] newAxis;
                    }
                    
                    //================================ And average the peak heights over all axes
                    curSum                           /= static_cast< proshade_double > ( ret->size() );
            
                    //================================ If improved, save
                    if ( curSum >  maxSum )
                    {
                        maxSum                        = curSum;
                        bestXRot                      = xChan;
                        bestYRot                      = yChan;
                        bestZRot                      = zChan;
                    }
                    
                    //================================ Release memory
                    delete[] rotMat;
                    
                    
                }
            }
        }
    
        //============================================ Prepare for next iteration
        searchRangeInDeg                             /= 2.0;
        finXRotChan                                  += bestXRot;
        finYRotChan                                  += bestYRot;
        finZRotChan                                  += bestZRot;
        bestXRot                                      = 0.0;
        bestYRot                                      = 0.0;
        bestZRot                                      = 0.0;
    }
    
    //================================================ Apply the optimisation
    rotMat                                            = ProSHADE_internal_maths::build3x3MatrixFromXYZRotations ( finXRotChan, finYRotChan, finZRotChan );
    for ( proshade_unsign axIt = 0; axIt < static_cast< proshade_unsign > ( ret->size() ); axIt++ )
    {
        //============================================ Find the rotated axis
        newAxis                                       = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMat, ret->at(axIt)[1], ret->at(axIt)[2], ret->at(axIt)[3] );
        
        //============================================ Change axes
        ret->at(axIt)[1]                              = newAxis[0];
        ret->at(axIt)[2]                              = newAxis[1];
        ret->at(axIt)[3]                              = newAxis[2];
    }
    
    //================================================ Release memory
    delete[] rotMat;
    
    //================================================ For each ret axis, compute predicted position
    for ( proshade_unsign axIt = 0; axIt < static_cast< proshade_unsign > ( ret->size() ); axIt++ )
    {
        //============================================ Convert XYZ to lat and lon INDICES
        lat                                           = std::atan2( ret->at(axIt)[2], ret->at(axIt)[1] ) / latSamlUnit;
        lon                                           = std::acos ( ret->at(axIt)[3] ) / lonSamlUnit;

        if ( lat < 0.0 )                              { lat += ( static_cast< proshade_double > ( dataObj->maxShellBand ) * 2.0 ); }
        if ( lon < 0.0 )                              { lon += ( static_cast< proshade_double > ( dataObj->maxShellBand ) * 2.0 ); }

        //============================================ Generate all angles for this fold
        applicableAngs.clear                          ( );
        for ( proshade_double angIt = 1.0; angIt < ret->at(axIt)[0]; angIt += 1.0 ) { ProSHADE_internal_misc::addToDoubleVector ( &applicableAngs, angIt * ( 2.0 * M_PI / ret->at(axIt)[0] ) ); }
        
        //============================================ For each shpere with the correct angle, average the peak heights
        ret->at(axIt)[5]                              = 0.0;
        for ( size_t angIt = 0; angIt < angs.size(); angIt++ )
        {
            //======================================== Find the correct sphere
            alreadyFound                              = false;
            for ( size_t aIt = 0; aIt < applicableAngs.size(); aIt++ ) { if ( alreadyFound ) { break; } const FloatingPoint< proshade_double > lhs1 ( angs.at(angIt) ), rhs1 ( applicableAngs.at(aIt) ); if ( lhs1.AlmostEquals ( rhs1 ) ) { alreadyFound = true; corAngIt = angIt; } }
            
            if ( !alreadyFound ) { continue; }
            
            //======================================== Get its peak height for the longitude and latitude
            ret->at(axIt)[5]                         += dataObj->sphereMappedRotFun.at(corAngIt)->getSphereLatLonLinearInterpolationPos ( lat, lon );
        }
        
        //============================================ And average the peak heights over the axis
        ret->at(axIt)[5]                             /= ( ret->at(axIt)[0] - 1.0 );
        maxSum                                       += ret->at(axIt)[5];
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, "Peak height detection and rotation optimisation complete.", settings->messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes two axes with almost dihedral angle and optimises their relative positions as well as orientation with respect to the optimal angle and the rotation function.
 
    This function has two parts. Firstly, it makes the assumption that the two axes supplied are almost perpendicular, but that there may be an error. To correct for this, the function will find
    the vector perpendicular to the plane formed by the two supplied axes and then it will find a vector perpendicular to the plane between this new vector and the axis with higher rotation
    function value. The resulting new vector is perpendicular to the axis with higher rotation function value and replaces the axis vector with the lower rotation function values (after normalisation).
 
    Secondly, this function supplies the two axes to the findPredictedAxesHeights() function, which does the maximisation of the rotation function values for both supplied axes at the same time.
    The resulting optimised axes are then returned in place of the inputted axes.
 
    \param[in] ret The list of axes for which the optimisation is to be done.
    \param[in] dataObj The structure object with computed rotation function in which the peaks are to be found.
    \param[in] settings ProSHADE_settings object containing all the settings for this run.
 */
void ProSHADE_internal_symmetry::optimiseDGroupAngleFromAxesHeights ( std::vector < std::vector< proshade_double > >* ret, ProSHADE_internal_data::ProSHADE_data* dataObj, ProSHADE_settings* settings )
{
    //================================================ Sanity check
    if ( ret->size() != 2 )
    {
        throw ProSHADE_exception ( "Attempted to optimise less than two axes for dihedral\n                    : group.", "ES00070", __FILE__, __LINE__, __func__, "The function for optimisation of dihedral angle of D\n                    : group was called on group with less than two axes. This\n                    : seems like a programming bug and should not happen - \n                    : contact author if you ever see this." );
    }
    
    //================================================ Set the angle to the correct dihedral group position - i.e. 90 deg
    proshade_double *crossProd, *perpVec, normFactor;
    size_t higherRFIndex                              = 0;
    
    // ... Find vector perperndicular to the plane given by the two axes
    crossProd                                         = ProSHADE_internal_maths::computeCrossProduct ( &ret->at(0).at(1), &ret->at(0).at(2), &ret->at(0).at(3), &ret->at(1).at(1), &ret->at(1).at(2), &ret->at(1).at(3) );
    
    // ... Find a vector perpendicular to the plane between the new vector and the vector with higher rotation function value
    if ( ret->at(1).at(5) > ret->at(0).at(5) ) { higherRFIndex = 1; }
    perpVec                                           = ProSHADE_internal_maths::computeCrossProduct ( &ret->at(higherRFIndex).at(1), &ret->at(higherRFIndex).at(2), &ret->at(higherRFIndex).at(3), &crossProd[0], &crossProd[1], &crossProd[2] );
    
    // ... Normalise the new vector
    normFactor                                        = std::sqrt ( pow ( perpVec[0], 2.0 ) + pow ( perpVec[1], 2.0 ) + pow ( perpVec[2], 2.0 ) );
    perpVec[0] /= normFactor; perpVec[1] /= normFactor; perpVec[2] /= normFactor;
    
    // ... Which vector are we to over-write?
    if ( higherRFIndex == 0 ) { higherRFIndex = 1; }
    else                      { higherRFIndex = 0; }
    
    // ... Set largest axis element to positive
    const FloatingPoint< proshade_double > lhs1 ( std::max ( std::abs ( perpVec[0] ), std::max( std::abs ( perpVec[1] ), std::abs ( perpVec[2] ) ) ) );
    const FloatingPoint< proshade_double > rhs1 ( std::abs ( perpVec[0] ));
    const FloatingPoint< proshade_double > rhs2 ( std::abs ( perpVec[1] ) );
    const FloatingPoint< proshade_double > rhs3 ( std::abs ( perpVec[2] ) );
    if ( ( lhs1.AlmostEquals ( rhs1 ) && ( perpVec[0] < 0.0 ) ) ||
         ( lhs1.AlmostEquals ( rhs2 ) && ( perpVec[1] < 0.0 ) ) ||
         ( lhs1.AlmostEquals ( rhs3 ) && ( perpVec[2] < 0.0 ) ) )
    {
        perpVec[0]                                   *= -1.0;
        perpVec[1]                                   *= -1.0;
        perpVec[2]                                   *= -1.0;
    }
    
    // ... Over-write the old vector with the better one
    ret->at(higherRFIndex).at(1) = perpVec[0]; ret->at(higherRFIndex).at(2) = perpVec[1]; ret->at(higherRFIndex).at(3) = perpVec[2];
    
    // ... Release memory
    delete[] perpVec;
    delete[] crossProd;
    
    //================================================ Convert input to pointers
    std::vector< proshade_double* > convVec;
    for ( size_t axIt = 0; axIt < 2; axIt++ )
    {
        //============================================ Allocate memory
        proshade_double* axVals                       = new proshade_double[7];
        ProSHADE_internal_misc::checkMemoryAllocation ( axVals, __FILE__, __LINE__, __func__ );
        
        //============================================ Copy values
        for ( size_t elIt = 0; elIt < 7; elIt++ )
        {
            axVals[elIt]                              = ret->at(axIt).at(elIt);
        }
        
        //============================================ Save
        convVec.push_back                             ( axVals );
    }
    
    //================================================ Run normal optimisation
    ProSHADE_internal_symmetry::findPredictedAxesHeights ( &convVec, dataObj, settings );
    
    //================================================ Convert back and release memory
    for ( size_t axIt = 0; axIt < 2; axIt++ )
    {
        //============================================ Copy values
        for ( size_t elIt = 0; elIt < 7; elIt++ )
        {
            ret->at(axIt).at(elIt)                     = convVec.at(axIt)[elIt];
        }
        
        //============================================ Release memory
        delete[] convVec.at(axIt);
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes two axes with almost dihedral angle and optimises their relative positions as well as orientation with respect to the optimal angle and the rotation function.
 
    This function is basically a wrapper to the overloaded function with the same name. It starts with converting the axes list and index vector to a single vactor containing the axes in the correct format
    and then it proceeds to call the overloaded version of this function to do the actual work. Once it completes, this function saves the results into the original input vector and terminates itself.
 
    \param[in] allCs The list of axes for a subset of which the optimisation is to be done.
    \param[in] selection A vector of two indices (of the allCs vector) specifying which axes to optimise.
    \param[in] dataObj The structure object with computed rotation function in which the peaks are to be found.
    \param[in] settings ProSHADE_settings object containing all the settings for this run.
 */
void ProSHADE_internal_symmetry::optimiseDGroupAngleFromAxesHeights ( std::vector < std::vector< proshade_double > >* allCs, std::vector< proshade_unsign > selection, ProSHADE_internal_data::ProSHADE_data* dataObj, ProSHADE_settings* settings )
{
    //================================================ Initialise local variables
    std::vector < std::vector< proshade_double > > ortPair;
    std::vector< proshade_double > hlpVec;
    
    //================================================ Convert the indices and the list into a single vector containing only the axes to be optimised.
    hlpVec.push_back ( allCs->at(selection.at(0))[0] ); hlpVec.push_back ( allCs->at(selection.at(0))[1] ); hlpVec.push_back ( allCs->at(selection.at(0))[2] );
    hlpVec.push_back ( allCs->at(selection.at(0))[3] ); hlpVec.push_back ( allCs->at(selection.at(0))[4] ); hlpVec.push_back ( allCs->at(selection.at(0))[5] );
    hlpVec.push_back ( allCs->at(selection.at(0))[6] );
    ortPair.push_back ( hlpVec ); hlpVec.clear ( );
    hlpVec.push_back ( allCs->at(selection.at(1))[0] ); hlpVec.push_back ( allCs->at(selection.at(1))[1] ); hlpVec.push_back ( allCs->at(selection.at(1))[2] );
    hlpVec.push_back ( allCs->at(selection.at(1))[3] ); hlpVec.push_back ( allCs->at(selection.at(1))[4] ); hlpVec.push_back ( allCs->at(selection.at(1))[5] );
    hlpVec.push_back ( allCs->at(selection.at(1))[6] );
    ortPair.push_back ( hlpVec );
    
    //================================================ Run the optimisation proper
    optimiseDGroupAngleFromAxesHeights                ( &ortPair, dataObj, settings );
    
    //================================================ Save the results back to the vector
    allCs->at(selection.at(0))[1] = ortPair.at(0).at(1); allCs->at(selection.at(0))[2] = ortPair.at(0).at(2); allCs->at(selection.at(0))[3] = ortPair.at(0).at(3);
    allCs->at(selection.at(0))[5] = ortPair.at(0).at(5); allCs->at(selection.at(0))[6] = ortPair.at(0).at(6);
    
    allCs->at(selection.at(1))[1] = ortPair.at(1).at(1); allCs->at(selection.at(1))[2] = ortPair.at(1).at(2); allCs->at(selection.at(1))[3] = ortPair.at(1).at(3);
    allCs->at(selection.at(1))[5] = ortPair.at(1).at(5); allCs->at(selection.at(1))[6] = ortPair.at(1).at(6);
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds the rotation function value for a single axis.
 
    This function is a simplified version of the findPredictedAxesHeights, except this one searches for the density map
    peak height for a single supplied axis (with the format of the array being x = [0], y = [1] and z = [2] and the fold being
    supplied separately).
 
    \param[in] axis A single axis for which the height is to be found.
    \param[in] fold The fold the axis should have.
    \param[in] dataObj The structure object with computed rotation function in which the peaks are to be found.
    \param[in] settings ProSHADE_settings object containing all the settings for this run.
    \param[out] height The height for this axis.
 */
proshade_double ProSHADE_internal_symmetry::findPredictedSingleAxisHeight ( proshade_double* axis, proshade_double fold, ProSHADE_internal_data::ProSHADE_data* dataObj, ProSHADE_settings* settings )
{
    //================================================ Initialise variables
    proshade_double height                            = 0.0;
    proshade_double lat, lon;
    proshade_double latSamlUnit                       = ( 2.0 * M_PI ) / ( static_cast< proshade_double > ( dataObj->maxShellBand ) * 2.0 );
    proshade_double lonSamlUnit                       = ( 1.0 * M_PI ) / ( static_cast< proshade_double > ( dataObj->maxShellBand ) * 2.0 );
    
    //================================================ Make sure we have a clean start
    dataObj->sphereMappedRotFun.clear                 ( );
    
    //================================================ Convert rotation function to only the required angle-axis space spheres and find all peaks
    for ( proshade_double angIt = 1.0; angIt < fold; angIt += 1.0 )
    {
        //============================================ Create the angle-axis sphere with correct radius (angle)
        dataObj->sphereMappedRotFun.emplace_back      ( new ProSHADE_internal_spheres::ProSHADE_rotFun_sphere ( angIt * ( 2.0 * M_PI / fold ),
                                                                                                                M_PI / fold,
                                                                                                                dataObj->maxShellBand * 2,
                                                                                                                angIt * ( 2.0 * M_PI / fold ),
                                                                                                                static_cast<proshade_unsign> ( angIt - 1.0 ) ) );
        
        //============================================ Interpolate rotation function onto the sphere
        dataObj->sphereMappedRotFun.at( static_cast < proshade_unsign > ( angIt - 1.0 ))->interpolateSphereValues ( dataObj->getInvSO3Coeffs ( ) );
    }
    
    //================================================ Convert XYZ to lat and lon INDICES
    lat                                               = std::atan2( axis[1], axis[0] ) / latSamlUnit;
    lon                                               = std::acos ( axis[2] ) / lonSamlUnit;
    
    if ( lat < 0.0 ) { lat += ( static_cast< proshade_double > ( dataObj->maxShellBand ) * 2.0 ); }
    if ( lon < 0.0 ) { lon += ( static_cast< proshade_double > ( dataObj->maxShellBand ) * 2.0 ); }
    
    lat                                               = std::round ( lat );
    lon                                               = std::round ( lon );
    
    //================================================ Initialise the peak group
    ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup* grp = nullptr;
    
    //================================================ Construct a peak group with entry from each sphere with the axis as the peak
    for ( proshade_unsign sphIt = 0; sphIt < static_cast<proshade_unsign> ( dataObj->sphereMappedRotFun.size() ); sphIt++ )
    {
        if ( sphIt == 0 )
        {
            //======================================== If first sphere, create the peak group
            grp                                       = new ProSHADE_internal_spheres::ProSHADE_rotFun_spherePeakGroup ( lat, lon, sphIt, dataObj->sphereMappedRotFun.at(sphIt)->getAngularDim() );
        }
        else
        {
            //======================================== Add to the existing object
            grp->checkIfPeakBelongs                   ( lat, lon, sphIt, settings->axisErrTolerance, settings->verbose, settings->messageShift );
        }
    }
    
    //================================================ Find the peak height
    std::vector < proshade_double* > detectedAxis;
    grp->findCyclicPointGroupsGivenFold               ( dataObj->sphereMappedRotFun, &detectedAxis, settings->useBiCubicInterpolationOnPeaks, static_cast< proshade_unsign > ( fold ), settings->verbose, settings->messageShift );
    
    //================================================ Save it!
    if ( detectedAxis.size() > 0 )                    { height = detectedAxis.at(0)[5]; }
    else                                              { height = 0.0; }
    
    //================================================ Release memory
    for ( proshade_unsign i = 0; i < static_cast < proshade_unsign > ( detectedAxis.size() ); i++ ) { delete detectedAxis.at(i); }
    delete grp;
    
    //================================================ Done
    return                                            ( height );
    
}

/*! \brief This function finds the best pair of axes conforming to the tetrahedron dihedral angle.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
    \param[in] axErr The error tolerance on angle matching.
    \param[out] ret The pair of axes with closest angle to the required icosahedron dihedral angle.
 */
std::pair< proshade_unsign, proshade_unsign > findBestTetraDihedralPair ( std::vector< proshade_double* >* CSymList, proshade_double minPeakHeight, proshade_double axErr )
{
    //================================================ Initialise variables
    std::pair< proshade_unsign, proshade_unsign > ret;
    std::vector< proshade_unsign  > C3List;
    proshade_double bestHeightSum                     = 0.0;
    proshade_double dotProduct;
    
    //================================================ Find all C3 symmetries
    for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ ) { const FloatingPoint< proshade_double > lhs1 ( CSymList->at(cSym)[0] ), rhs1 ( 3.0 ); if ( lhs1.AlmostEquals ( rhs1 ) && CSymList->at(cSym)[5] >= minPeakHeight ) { ProSHADE_internal_misc::addToUnsignVector ( &C3List, cSym ); } }
    
    //================================================ For each unique pair of C3 and C2
    for ( proshade_unsign c3 = 0; c3 < static_cast<proshade_unsign> ( C3List.size() ); c3++ )
    {
        for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ )
        {
            //======================================== Compare only C2s to the C3List and only with decent average peak height
            const FloatingPoint< proshade_double > lhs1 ( CSymList->at(cSym)[0] ), rhs1 ( 2.0 );
            if ( !lhs1.AlmostEquals ( rhs1 ) ) { continue; }
            if ( CSymList->at(cSym)[5] < minPeakHeight ) { continue; }
            
            //========================================  Check the angle between the C5 and C3 axes
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &CSymList->at(C3List.at(c3))[1],
                                                                                                     &CSymList->at(C3List.at(c3))[2],
                                                                                                     &CSymList->at(C3List.at(c3))[3],
                                                                                                     &CSymList->at(cSym)[1],
                                                                                                     &CSymList->at(cSym)[2],
                                                                                                     &CSymList->at(cSym)[3] );
            
            //======================================== Is the angle approximately the dihedral angle?
            if ( ( ( 1.0 / sqrt ( 3.0 ) ) > ( std::abs( dotProduct ) - axErr ) ) && ( ( 1.0 / sqrt ( 3.0 ) ) < ( std::abs( dotProduct ) + axErr ) ) )
            {
                if ( bestHeightSum < ( CSymList->at(C3List.at(c3))[5] + CSymList->at(cSym)[5] ) )
                {
                    bestHeightSum                     = CSymList->at(C3List.at(c3))[5] + CSymList->at(cSym)[5];
                    ret.first                         = C3List.at(c3);
                    ret.second                        = cSym;
                }
            }
        }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function predicts a list of all T symmetry axes from the already computed C symmetries list.
 
    This function starts by checking if there is a pair of C3 and C2 symmetries with the tetrahedron dihedral angle ( acos( ( 1.0 / sqrt ( 3.0 ) ) ). If found,
    it calls the  predictTetraAxes() function, which uses the knowledge of the two axes (C3 and C2) which are closest to the dihedral angle to find the best rotation matrix matching a
    pre-computed tetrhedron model to the detected axes. After rotating the model, the model axes become the predicted axes for the structure and their peak heights are then
    computed. Once complete, all the predicted axes are in the ret variable.
 
    \warning This function does not check if the correct number of C axes was found, it is assumed this will be checked when the determination of
    which symmetry was found.
 
    \param[in] settings A pointer to settings class containing all the information required for symmetry detection.
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[out] ret A vector of all the detected axes in the standard ProSHADE format with height either the detected value (for the detected ones) or 0 for the predicted ones.
 */
std::vector< proshade_double* > ProSHADE_internal_data::ProSHADE_data::getPredictedTetrahedralSymmetriesList ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSymList )
{
    //================================================ Initialise variables
    std::vector< proshade_double* > ret;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting T symmetry prediction.", settings->messageShift );
    
    //================================================ Are the basic requirements for icosahedral symmetry met?
    if ( ProSHADE_internal_symmetry::detectTetrahedralSymmetry ( CSymList, settings->axisErrTolerance, settings->minSymPeak ) )
    {
        //============================================ Generate the rest of the axes
        ProSHADE_internal_symmetry::predictTetraAxes  ( CSymList, &ret, settings->axisErrTolerance, settings->minSymPeak );

        //============================================ Get heights for the predicted axes
        ProSHADE_internal_symmetry::findPredictedAxesHeights ( &ret, this, settings );

        //============================================ Add predicted axes to detected C axes list and also to the settings Icosahedral symmetry list
        for ( proshade_unsign retIt = 0; retIt < static_cast < proshade_unsign > ( ret.size() ); retIt++ )
        {
            proshade_signed matchedPos                = ProSHADE_internal_symmetry::addAxisUnlessSame ( static_cast< proshade_unsign > ( ret.at(retIt)[0] ), ret.at(retIt)[1], ret.at(retIt)[2], ret.at(retIt)[3], ret.at(retIt)[5], CSymList, settings->axisErrTolerance );
            ProSHADE_internal_misc::addToUnsignVector ( &settings->allDetectedTAxes, static_cast < proshade_unsign > ( matchedPos ) );
        }
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "T symmetry prediction complete.", settings->messageShift );

    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function predicts all tetrahedral point group symmetry axes from the cyclic point groups list.
 
    This function starts by finding the rotation matrix corresponding to the angle between the predicted C3 axis and the C3 axis found in the pre-computed
    tetrahedron model available in the ProSHADE_precomputedValues file. It then proceeds to use this rotation matrix to rotate the pre-computed model C2
    axis to now be in the correct orientation to the detected C3 axis.
 
    Next, the function computes the rotation matrix corresponding to rotation along the detected C3 axis about the angle between the rotated pre-computed model
    C2axis and the detected C2 axis. Finally, when these two rotation matrices are combined, the resulting rotation matrix is the optimal match rotation between the
    pre-computed model and the detected axes positions. This final rotation matrix is then used to rotate the model axes and these rotated model axes are then the
    predicted axes in the structre.
 
    Please note that the peak heights are set to 0.0 for all predicted axes, as they were not detected in the structure, but were only predicted.
 
    \warning This function assumes that the detectIcosahedralSymmetry() function has successfully run (i.e. returned true).
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] ret The vector containing all the axes forming icosahedral group or empty vector.
    \param[in] axErr The error tolerance on angle matching.
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
 */
void ProSHADE_internal_symmetry::predictTetraAxes ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, proshade_double minPeakHeight )
{
    //================================================ Create the tetrahedronAxes object
    ProSHADE_internal_precomputedVals::tetrahedronAxes *tetAx = new ProSHADE_internal_precomputedVals::tetrahedronAxes ( );
    
    //================================================ Find the best axis combination with dihedral angle and correct folds
    std::pair< proshade_unsign, proshade_unsign > initAxes = findBestTetraDihedralPair ( CSymList, minPeakHeight, axErr );

    //================================================ Find rotation between the detected C3 and the model C3 axes.
    proshade_double* rotMat                           = ProSHADE_internal_maths::findRotMatMatchingVectors ( tetAx->getValue ( 0, 1 ),
                                                                                                             tetAx->getValue ( 0, 2 ),
                                                                                                             tetAx->getValue ( 0, 3 ),
                                                                                                             CSymList->at(initAxes.first)[1],
                                                                                                             CSymList->at(initAxes.first)[2],
                                                                                                             CSymList->at(initAxes.first)[3] );

    //================================================ Rotate the model C2 to the correct orientation relative to the detected C3 axis.
    proshade_double* rotModelC2                       = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMat,
                                                                                                                        tetAx->getValue ( 4, 1 ),
                                                                                                                        tetAx->getValue ( 4, 2 ),
                                                                                                                        tetAx->getValue ( 4, 3 ) );

    //================================================ Find the angle betwen the rotated model C2 and the detected C2 axes along the detected C3 axis.
    proshade_double bestAng = 0.0, curAngDist, bestAngDist = 999.9;
    proshade_double* rotMatHlp                        = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMatHlp, __FILE__, __LINE__, __func__ );
    for ( proshade_double ang = 0.0; ang < ( M_PI * 2.0 ); ang += 0.002 )
    {
        //============================================ Compute rotation matrix for this angle value
        ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( rotMatHlp, CSymList->at(initAxes.first)[1], CSymList->at(initAxes.first)[2], CSymList->at(initAxes.first)[3], ang );
        
        //============================================ Rotate the rotated C2 by the matrix
        proshade_double* rotRotModelC2                = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMatHlp,
                                                                                                                        rotModelC2[0],
                                                                                                                        rotModelC2[1],
                                                                                                                        rotModelC2[2] );
        
        //============================================ Find distance
        curAngDist                                    = std::sqrt ( std::pow ( rotRotModelC2[0] - CSymList->at(initAxes.second)[1], 2.0 ) +
                                                                    std::pow ( rotRotModelC2[1] - CSymList->at(initAxes.second)[2], 2.0 ) +
                                                                    std::pow ( rotRotModelC2[2] - CSymList->at(initAxes.second)[3], 2.0 ) );
        
        //============================================ Save best angle
        if ( curAngDist < bestAngDist ) { bestAngDist = curAngDist; bestAng = ang; }
        
        //============================================ Release memory
        delete[] rotRotModelC2;
    }
    
    //================================================ Release memory
    delete[] rotMatHlp;
    
    //================================================ For the rotation matrix along the detected C5 axis with the same anlge as is between the rotated model C3 and the detected C3 axes.
    proshade_double* rotMat2                          = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat2, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( rotMat2, CSymList->at(initAxes.first)[1], CSymList->at(initAxes.first)[2], CSymList->at(initAxes.first)[3], bestAng );
    
    //================================================ Combine the two rotation matrices into a single rotation matrix
    proshade_double* rotMatFin                        = ProSHADE_internal_maths::compute3x3MatrixMultiplication ( rotMat2, rotMat );
    
    //================================================ For each model axis
    for ( proshade_unsign iter = 0; iter < tetAx->getNoAxes( ); iter++ )
    {
        //============================================ Rotate the model axis to fit the detected orientation
        proshade_double* rotAxis                      = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMatFin,
                                                                                                                        tetAx->getValue ( iter, 1 ),
                                                                                                                        tetAx->getValue ( iter, 2 ),
                                                                                                                        tetAx->getValue ( iter, 3 ) );

        //============================================ Create ProSHADE symmetry axis representation
        proshade_double* axis                         = new proshade_double[7];
        ProSHADE_internal_misc::checkMemoryAllocation ( axis, __FILE__, __LINE__, __func__ );

        axis[0]                                       = tetAx->getValue ( iter, 0 );
        axis[1]                                       = rotAxis[0];
        axis[2]                                       = rotAxis[1];
        axis[3]                                       = rotAxis[2];
        axis[4]                                       = ( 2.0 * M_PI ) / axis[0];
        axis[5]                                       = 0.0;
        axis[6]                                       = -std::numeric_limits < proshade_double >::infinity();

        //============================================ Save axis to ret
        ProSHADE_internal_misc::addToDblPtrVector     ( ret, axis );

        //============================================ Release memory
        delete[] rotAxis;
    }

    //================================================ Release memory
    delete[] rotMat;
    delete[] rotMat2;
    delete[] rotMatFin;
    delete[] rotModelC2;
    delete   tetAx;

    //================================================ Done
    return ;
    
}

/*! \brief This function checks the list of detected axes (presumably from phaseless symmetry detection) and returns the best dihedral (or cyclic, if no dihedral is found) point
    group, or empty vector if nothing is found.
 
    This function starts by computing a stringent rotation function height threshold and proceeds to search for orthogonal pair of axes with at least this threshold RF value.
    If more than one pair is found, the pair with the highest sum of RF value and FSC value is chosen and returned, while if no pair is found, a single axis passing the threshold
    (and with the highest RF + FSC value sum) is returned. If no axis passes, empty vector is returned.
 
    \warning This function is intended to be used internally and specifically with phase-less symmetry detection, but it does not check for these conditions to be met!
 
    \param[in] allCs A list of all detected symmetries in the phase-less map.
    \param[in] verbose How loud the function should be in the standard output?
    \param[in] messageShift Are we in a subprocess, so that the log should be shifted for this function call? If so, by how much?
    \param[in] tolerance What is the tolerance on the perpendicularity of two axes in terms of the dot product?
    \param[out] ret A vector containing a) zero entries if no reliable symmetry axis was found, b) a single symmetry axis index if only a reliable cyclic symmetry axis was found or c)
                    two axes indices in the case where two perpendicular reliable axes were detected.
 */
std::vector< proshade_unsign > ProSHADE_internal_symmetry::findReliableUnphasedSymmetries ( std::vector < std::vector< proshade_double > >* allCs, proshade_signed verbose, proshade_signed messageShift, proshade_double tolerance )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Deciding whether any axis(es) is/are reliable.", messageShift );
    
    for ( int i = 0; i < allCs->size(); i++ ) { std::cout << allCs->at(i).at(0) << " | " << allCs->at(i).at(1) << " x " << allCs->at(i).at(2) << " x " << allCs->at(i).at(3) << " | " << allCs->at(i).at(5) << " | " << allCs->at(i).at(6) << std::endl; }
    
    //================================================ Initialise local variables
    std::vector< proshade_unsign > ret;
    
    //================================================ Find the threshold
    proshade_double bestHistPeakStart                 = ProSHADE_internal_maths::findTopGroupSmooth ( allCs, 5, 0.01, 0.03, 5 );
    proshade_double bestFSCPeakStart                  = ProSHADE_internal_maths::findTopGroupSmooth ( allCs, 6, 0.01, 0.005, 5, 0.94 );
    std::cout << "Peak threshold: " << bestHistPeakStart << std::endl << "FSC Threshold: " << bestFSCPeakStart << std::endl;
    if ( bestHistPeakStart > 0.9 ) { bestHistPeakStart = 0.9; }
    
    //================================================ Are any axes orthogonal
    proshade_double dotProduct, maxOrtSum = 0.0, curOrtSum = 0.0;
    proshade_unsign maxOrtAx1 = 0, maxOrtAx2 = 0;
    for ( size_t relAx1 = 0; relAx1 < allCs->size(); relAx1++ )
    {
        //============================================ Consider only reliable axes
        if ( allCs->at(relAx1)[5] < bestHistPeakStart ) { continue; }
        if ( allCs->at(relAx1)[6] < bestFSCPeakStart  ) { continue; }
        
        for ( size_t relAx2 = 1; relAx2 < allCs->size(); relAx2++ )
        {
            //======================================== Ignore same axes
            if ( relAx1 >= relAx2 ) { continue; }
            
            //======================================== Consider only reliable axes
            if ( allCs->at(relAx2)[5] < bestHistPeakStart ) { continue; }

            //======================================== Are the two axes orthogonal?
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &allCs->at(relAx1)[1], &allCs->at(relAx1)[2],
                                                                                                     &allCs->at(relAx1)[3], &allCs->at(relAx2)[1],
                                                                                                     &allCs->at(relAx2)[2], &allCs->at(relAx2)[3] );
            
            //======================================== If close to zero, these two axes are perpendicular
            if ( std::abs( dotProduct ) < tolerance )
            {
                //==================================== Find sum
                curOrtSum                             = allCs->at(relAx1)[5] + allCs->at(relAx1)[6] + allCs->at(relAx2)[5] + allCs->at(relAx2)[6];

                //==================================== If best, save it
                if ( curOrtSum > maxOrtSum )
                {
                    maxOrtSum                         = curOrtSum;
                    maxOrtAx1                         = static_cast< proshade_unsign > ( relAx1 );
                    maxOrtAx2                         = static_cast< proshade_unsign > ( relAx2 );
                }
            }
        }
    }
    
    //================================================ If any orthogonal pair was found, return it
    if ( maxOrtAx2 != 0 )
    {
        //================================================ Report progress
        ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "Found two, orthogonal and reliable axes.", messageShift );
        std::cout << allCs->at(maxOrtAx1)[0] << " | " << allCs->at(maxOrtAx1)[1] << " x " << allCs->at(maxOrtAx1)[2] << " x " << allCs->at(maxOrtAx1)[3] << " || " << allCs->at(maxOrtAx1)[5] << " || " << allCs->at(maxOrtAx1)[6] << std::endl;
        std::cout << allCs->at(maxOrtAx2)[0] << " | " << allCs->at(maxOrtAx2)[1] << " x " << allCs->at(maxOrtAx2)[2] << " x " << allCs->at(maxOrtAx2)[3] << " || " << allCs->at(maxOrtAx2)[5] << " || " << allCs->at(maxOrtAx2)[6] << std::endl;
        
        
        exit(0); // REMOVE ALSO THE ONE AT THE END OF THE FUNCTION!
        
        ProSHADE_internal_misc::addToUnsignVector     ( &ret, maxOrtAx1 );
        ProSHADE_internal_misc::addToUnsignVector     ( &ret, maxOrtAx2 );
        return                                        ( ret );
    }
    
    //================================================ Well, no orthogonal axes. Is there at least one good axis?
    curOrtSum = 0.0; maxOrtSum = 0.0;
    for ( size_t relAx1 = 0; relAx1 < allCs->size(); relAx1++ )
    {
        //============================================ Consider only reliable axes in terms of RF
        if ( allCs->at(relAx1)[5] < bestHistPeakStart ) { continue; }
        
        //============================================ Consider only reasonable axes in terms of FSC
        if ( allCs->at(relAx1)[6] < 0.3 ) { continue; }
        
        //============================================ Get the sum
        curOrtSum                                     = allCs->at(relAx1)[5] + allCs->at(relAx1)[6];
        
        //============================================ If highest sum, save
        if ( curOrtSum > maxOrtSum )
        {
            maxOrtSum                                 = curOrtSum;
            maxOrtAx1                                 = static_cast< proshade_unsign > ( relAx1 );
        }
    }
    
    //================================================ If anything was found, save it
    if ( maxOrtSum > 0.1 )
    {
        //============================================ Report progress
        ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "Found single reliable axis.", messageShift );
        
        ProSHADE_internal_misc::addToUnsignVector     ( &ret, maxOrtAx1 );
    }
    else
    {
        //================================================ Report progress
        ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "Found no reliable axis.", messageShift );
    }
    exit(0);
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function allocates the required memory for the Fourier transforms required to find the centre of the map.
 
    \param[in] xDim The size of the x-axis in indices.
    \param[in] yDim The size of the y-axis in indices.
    \param[in] zDim The size of the z-axis in indices.
    \param[in] origMap Array to which the original map will be saved before Fourier transform computation.
    \param[in] origCoeffs Array to which the result of the Fourier transform of the original map will be saved into.
    \param[in] rotMapComplex Array to which the rotated map will be saved before Fourier transform computation.
    \param[in] rotCoeffs Array to which the result of the Fourier transform of the rotated map will be saved into.
    \param[in] trFunc Array to which the results of inverse Fourier transform of the conbined coefficients will be saved.
    \param[in] trFuncCoeffs Array to which the two maps coefficients will be combined into before inverse Fourier transform computation.
    \param[in] planForwardFourier FFTW3 plan for the original map Fourier transform.
    \param[in] planForwardFourierRot FFTW3 plat for the rotated map Fourier transform.
    \param[in] planReverseFourierComb FFTW3 plan for the inverse Fourier transform from the combined coefficients to translation function.
 */
void ProSHADE_internal_symmetry::allocateCentreOfMapFourierTransforms ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, fftw_complex *&origMap, fftw_complex *&origCoeffs, fftw_complex *&rotMapComplex, fftw_complex *&rotCoeffs, fftw_complex *&trFunc, fftw_complex *&trFuncCoeffs, fftw_plan *planForwardFourier, fftw_plan *planForwardFourierRot, fftw_plan *planReverseFourierComb )
{
    //================================================ Allocate the memory
    origMap                                           = new fftw_complex [xDim * yDim * zDim];
    origCoeffs                                        = new fftw_complex [xDim * yDim * zDim];
    rotMapComplex                                     = new fftw_complex [xDim * yDim * zDim];
    rotCoeffs                                         = new fftw_complex [xDim * yDim * zDim];
    trFunc                                            = new fftw_complex [xDim * yDim * zDim];
    trFuncCoeffs                                      = new fftw_complex [xDim * yDim * zDim];
    
    //================================================ Check the memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( origMap,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( origCoeffs,    __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMapComplex, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotCoeffs,     __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( trFunc,        __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( trFuncCoeffs,  __FILE__, __LINE__, __func__ );
    
    //================================================
   *planForwardFourier                                = fftw_plan_dft_3d ( static_cast< int > ( xDim ), static_cast< int > ( yDim ), static_cast< int > ( zDim ), origMap,       origCoeffs, FFTW_FORWARD,   FFTW_ESTIMATE );
   *planForwardFourierRot                             = fftw_plan_dft_3d ( static_cast< int > ( xDim ), static_cast< int > ( yDim ), static_cast< int > ( zDim ), rotMapComplex, rotCoeffs,  FFTW_FORWARD,   FFTW_ESTIMATE );
   *planReverseFourierComb                            = fftw_plan_dft_3d ( static_cast< int > ( xDim ), static_cast< int > ( yDim ), static_cast< int > ( zDim ), trFuncCoeffs,  trFunc,     FFTW_BACKWARD,  FFTW_ESTIMATE );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function releases the allocated memory for the Fourier transforms used to find the centre of the map.
 
    \param[in] origMap Array to which the original map will be saved before Fourier transform computation.
    \param[in] origCoeffs Array to which the result of the Fourier transform of the original map will be saved into.
    \param[in] rotMapComplex Array to which the rotated map will be saved before Fourier transform computation.
    \param[in] rotCoeffs Array to which the result of the Fourier transform of the rotated map will be saved into.
    \param[in] trFunc Array to which the results of inverse Fourier transform of the conbined coefficients will be saved.
    \param[in] trFuncCoeffs Array to which the two maps coefficients will be combined into before inverse Fourier transform computation.
    \param[in] planForwardFourier FFTW3 plan for the original map Fourier transform.
    \param[in] planForwardFourierRot FFTW3 plat for the rotated map Fourier transform.
    \param[in] planReverseFourierComb FFTW3 plan for the inverse Fourier transform from the combined coefficients to translation function.
 */
void ProSHADE_internal_symmetry::releaseCentreOfMapFourierTransforms ( fftw_complex *origMap, fftw_complex *origCoeffs, fftw_complex *rotMapComplex, fftw_complex *rotCoeffs, fftw_complex *trFunc, fftw_complex *trFuncCoeffs, fftw_plan planForwardFourier, fftw_plan planForwardFourierRot, fftw_plan planReverseFourierComb )
{
    //================================================ Destroy the FFTW3 plans
    fftw_destroy_plan                                 ( planReverseFourierComb );
    fftw_destroy_plan                                 ( planForwardFourier );
    fftw_destroy_plan                                 ( planForwardFourierRot );
    
    //================================================ Set pointers to NULL
    planReverseFourierComb                            = nullptr;
    planForwardFourier                                = nullptr;
    planForwardFourierRot                             = nullptr;
    
    //================================================ Release the memory
    delete[] origMap;
    delete[] origCoeffs;
    delete[] rotMapComplex;
    delete[] rotCoeffs;
    delete[] trFunc;
    delete[] trFuncCoeffs;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes a single rotation matrix and procceds to compute the optimal translation between the original map and a map rotated by the supplied rotation matrix.
 
    \param[in] symStr A ProSHADE_data structure containing the structure for which the line on which the centre of rotation lies is to be found.
    \param[in] symElems Vector containing single symmetry element (rotaiton matrix) which is not identity.
    \param[in] origCoeffs The Fourier coefficients of the original (non-rotated) map.
    \param[in] rotMapComplex Array to which the rotated map will be saved and from which the Fourier transform plan (planForwardFourierRot) is prepared.
    \param[in] rotCoeffs Array to which the result of the Fourier transform of the rotated map will be saved into by the supplied plan (planForwardFourierRot).
    \param[in] planForwardFourierRot FFTW3 plan for forward Fourier transform from rotMapComplex to rotCoeffs.
    \param[in] trFuncCoeffs The array to which the combined Fourier coefficients for translation function will be saved into and also for which the inverse Fourier transform plan (planReverseFourierComb) is prepared for.
    \param[in] trFunc The array to which the translation function will be saved into by the reverse Fourier transform planned by the plan (planReverseFourierComb).
    \param[in] planReverseFourierComb FFTW3 plan for reverse Fourier transform from the combined coefficients (trFuncCoeffs) to the translation function array (trFunc).
    \param[in] verbose How loud the function should be?
    \param[out] trsVec A vector containing the optimal translation between the original and the rotated maps in Angstroms.
*/
std::vector< proshade_double > ProSHADE_internal_symmetry::findTranslationBetweenRotatedAndOriginalMap ( ProSHADE_internal_data::ProSHADE_data* symStr, std::vector < proshade_double > symElem, fftw_complex *origCoeffs, fftw_complex* rotMapComplex, fftw_complex* rotCoeffs, fftw_plan planForwardFourierRot, fftw_complex* trFuncCoeffs, fftw_complex* trFunc, fftw_plan planReverseFourierComb, proshade_signed verbose )
{
    //================================================ Initialise local variables
    proshade_double axX, axY, axZ, axAng, mapPeak, trsX, trsY, trsZ;
    std::vector< proshade_double > trsVec;
    
    //=============================================== Rotate the map by the rotation matrix
    proshade_double *rotMap;
    ProSHADE_internal_maths::getAxisAngleFromRotationMatrix ( &symElem, &axX, &axY, &axZ, &axAng );
    symStr->rotateMapRealSpace                        ( axX, axY, axZ, axAng, rotMap );
    
    //================================================ Convert rotated map to Fourier space
    for ( size_t it = 0; it < static_cast< size_t > ( symStr->getXDim() * symStr->getYDim() * symStr->getZDim() ); it++ ) { rotMapComplex[it][0] = rotMap[it]; rotMapComplex[it][1] = 0.0; }
    fftw_execute                                      ( planForwardFourierRot );

    //================================================ Combine coeffs for translation function
    ProSHADE_internal_maths::combineFourierForTranslation ( origCoeffs, rotCoeffs, trFuncCoeffs, symStr->getXDim(), symStr->getYDim(), symStr->getZDim() );

    //================================================ Compute translation function
    fftw_execute                                      ( planReverseFourierComb );

    //================================================ Find peak
    mapPeak                                           = 0.0;
    ProSHADE_internal_maths::findHighestValueInMap    ( trFunc, symStr->getXDim(), symStr->getYDim(), symStr->getZDim(), &trsX, &trsY, &trsZ, &mapPeak );
    
    //================================================ Convert to Angstroms
    trsX                                             *= ( symStr->getXDimSize() / symStr->getXDim() );
    trsY                                             *= ( symStr->getYDimSize() / symStr->getYDim() );
    trsZ                                             *= ( symStr->getZDimSize() / symStr->getZDim() );
    
    //================================================ Do not translate over half
    if ( trsX > ( static_cast< proshade_double > ( symStr->getXDimSize() ) / 2.0 ) ) { trsX = trsX - static_cast< proshade_double > ( symStr->getXDimSize() ); }
    if ( trsY > ( static_cast< proshade_double > ( symStr->getYDimSize() ) / 2.0 ) ) { trsY = trsY - static_cast< proshade_double > ( symStr->getYDimSize() ); }
    if ( trsZ > ( static_cast< proshade_double > ( symStr->getZDimSize() ) / 2.0 ) ) { trsZ = trsZ - static_cast< proshade_double > ( symStr->getZDimSize() ); }

    //================================================ Save line point
    ProSHADE_internal_misc::addToDoubleVector         ( &trsVec, trsX );
    ProSHADE_internal_misc::addToDoubleVector         ( &trsVec, trsY );
    ProSHADE_internal_misc::addToDoubleVector         ( &trsVec, trsZ );
    
    //================================================ Release memory
    delete[] rotMap;

    //================================================ Done
    return                                            ( trsVec );
    
}

/*! \brief This function computes the average of optimal translations for a cyclic point group.
 
    This function takes a single cyclic point group elements and proceeds to compute all optimal translations between the original map and map
    rotated by each point group element. The sum of these translations divided by the number of the point group elements (including the identity
    element) then gives a point that must lie on the symmetry axis.
 
    \param[in] symStr A ProSHADE_data structure containing the structure for which the line on which the centre of rotation lies is to be found.
    \param[in] symElems Vector containing single symmetry element (rotaiton matrix) which is not identity.
    \param[in] origCoeffs The Fourier coefficients of the original (non-rotated) map.
    \param[in] rotMapComplex Array to which the rotated map will be saved and from which the Fourier transform plan (planForwardFourierRot) is prepared.
    \param[in] rotCoeffs Array to which the result of the Fourier transform of the rotated map will be saved into by the supplied plan (planForwardFourierRot).
    \param[in] planForwardFourierRot FFTW3 plan for forward Fourier transform from rotMapComplex to rotCoeffs.
    \param[in] trFuncCoeffs The array to which the combined Fourier coefficients for translation function will be saved into and also for which the inverse Fourier transform plan (planReverseFourierComb) is prepared for.
    \param[in] trFunc The array to which the translation function will be saved into by the reverse Fourier transform planned by the plan (planReverseFourierComb).
    \param[in] planReverseFourierComb FFTW3 plan for reverse Fourier transform from the combined coefficients (trFuncCoeffs) to the translation function array (trFunc).
    \param[in] verbose How loud the function should be?
    \param[out] pointOnLine A vector specifying a point that lies on the symmetry axis (given by the averaged sum of the translations of the rotated maps).
*/
std::vector< proshade_double > ProSHADE_internal_symmetry::findPointFromTranslations ( ProSHADE_internal_data::ProSHADE_data* symStr, std::vector < std::vector < proshade_double > > symElems, fftw_complex *origCoeffs, fftw_complex* rotMapComplex, fftw_complex* rotCoeffs, fftw_plan planForwardFourierRot, fftw_complex* trFuncCoeffs, fftw_complex* trFunc, fftw_plan planReverseFourierComb, proshade_signed verbose )
{
    //================================================ Initialise local variables
    std::vector< proshade_double > pointOnLine        ( 3, 0.0 );
    std::vector< proshade_double > identityMat        ( 9, 0.0 ); identityMat.at(0) = 1.0; identityMat.at(4) = 1.0; identityMat.at(8) = 1.0;
    
    //================================================ For each symmetry element in this cyclic group
    for ( size_t gEl = 0; gEl < symElems.size(); gEl++ )
    {
        //============================================ Ignore identity element
        if ( ProSHADE_internal_maths::rotationMatrixSimilarity ( &symElems.at(gEl), &identityMat, 0.01 ) ) { continue; }
        
        //============================================ Find translation difference between rotated and original map
        std::vector< proshade_double > trsCenHlp      = ProSHADE_internal_symmetry::findTranslationBetweenRotatedAndOriginalMap ( symStr,
                                                                                                                                  symElems.at(gEl),
                                                                                                                                  origCoeffs, rotMapComplex,
                                                                                                                                  rotCoeffs, planForwardFourierRot,
                                                                                                                                  trFuncCoeffs, trFunc,
                                                                                                                                  planReverseFourierComb,
                                                                                                                                  verbose );
        
        //============================================ Sum translations over the whole axis
        pointOnLine.at(0)                            += trsCenHlp.at(0);
        pointOnLine.at(1)                            += trsCenHlp.at(1);
        pointOnLine.at(2)                            += trsCenHlp.at(2);
    }
    
    //================================================ Average over all symmetry elements (including the identity one)
    pointOnLine.at(0)                                /= static_cast< proshade_double > ( symElems.size() );
    pointOnLine.at(1)                                /= static_cast< proshade_double > ( symElems.size() );
    pointOnLine.at(2)                                /= static_cast< proshade_double > ( symElems.size() );
    
    //================================================ Done
    return                                            ( pointOnLine );
    
}
