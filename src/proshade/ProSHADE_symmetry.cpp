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
    \version   0.7.4.4
    \date      OCT 2020
 */

//==================================================== ProSHADE
#include "ProSHADE_symmetry.hpp"

/*! \brief This function computes the self-rotation function for this structure.
 
    This function assumes that the spherical harmonics have been computed for a data object. It can be then called on this
    object and it will proceed to compute the E matrices for this object against itself. From these "self E matrices", the
    function will generate the SO(3) transform coefficients and finally it will invert transform these coefficients back,
    thus getting the self-rotation function.
 
    \param[in] settings A pointer to settings class containing all the information required for map self-rotation function computation.
 */
void ProSHADE_internal_data::ProSHADE_data::getRotationFunction ( ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting self-rotation function computation." );
    
    //================================================ Compute un-weighted E matrices and their weights
    ProSHADE_internal_distances::computeEMatrices     ( this, this, settings );
    
    //================================================ Normalise E matrices by the magnitudes
    ProSHADE_internal_distances::normaliseEMatrices   ( this, this, settings );
    
    //================================================ Generate SO(3) coefficients
    ProSHADE_internal_distances::generateSO3CoeffsFromEMatrices ( this, this, settings );
    
    //================================================ Compute the inverse SO(3) Fourier Transform (SOFT) on the newly computed coefficients
    ProSHADE_internal_distances::computeInverseSOFTTransform ( this, this, settings );
    
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Self-rotation function obtained." );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function obtains a list of all C symmetries from already computed self-rotation map.
 
    This function starts by finding all peaks in the self-rotation map, which are outliers in terms of height. It then proceeds to
    group these by the height, searching for C symmetries in each peak height group (thus making sure symmetries with higher peak heights
    are found first). The symmetry detection proceeds by detecting possible C symmetry folds and searching whether the all peaks are present
    to support the prospective C symmetry. If only few are missing, it will even search for the missing peaks. Finally, the function returns
    all detected symmetries in the order of decreasing average peak height.
 
    \param[in] settings A pointer to settings class containing all the information required for symmetry detection.
 */
std::vector< proshade_double* > ProSHADE_internal_data::ProSHADE_data::getCyclicSymmetriesList ( ProSHADE_settings* settings )
{
    //================================================ Initialise variables
    std::vector< proshade_double* > ret;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting C symmetry detection." );
    
    //================================================ Get list of peaks in the self-rotation map
    std::vector< proshade_double* > allPeaks          = ProSHADE_internal_peakSearch::getAllPeaksNaive ( this->getInvSO3Coeffs (), this->getMaxBand() * 2,
                                                                                                         settings->peakNeighbours,
                                                                                                         settings->noIQRsFromMedianNaivePeak );
    
    //================================================ Convert peaks to angle-axis
    std::vector< proshade_double* > peaksAA           = ProSHADE_internal_symmetry::getPeaksAngleAxisPositions ( allPeaks, settings->verbose );
    
    //================================================ Sort peaks by height groups
    std::vector< proshade_double > peakGroupsBoundaries = ProSHADE_internal_symmetry::findPeaksByHeightBoundaries ( peaksAA, settings->smoothingFactor );
    
    //================================================ Get symmetry per group
    std::vector< std::vector< proshade_unsign > > detectedCSymmetries;
    for ( proshade_signed iter = static_cast<proshade_unsign> ( peakGroupsBoundaries.size() - 1 ); iter >= 0; iter-- )
    {
        //============================================ Get peaks group peaks only
        std::vector< proshade_double* > symPeaks;
        for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( peaksAA.size() ); it++ )
        {
            if ( peaksAA.at(it)[4] > peakGroupsBoundaries.at(iter) ) { ProSHADE_internal_misc::addToDblPtrVector ( &symPeaks, peaksAA.at(it) ); }
        }
        
        //============================================ Search for symmetry in these peaks
        detectedCSymmetries                           = ProSHADE_internal_symmetry::findPeaksCSymmetry ( &symPeaks, settings->verbose, this->getMaxBand(),
                                                                                                         settings->symMissPeakThres,
                                                                                                         settings->axisErrTolerance, settings->axisErrToleranceDefault, this );
        
        //============================================ Print detected symmetries
        for ( proshade_unsign detIt = 0; detIt < static_cast<proshade_unsign> ( detectedCSymmetries.size() ); detIt++ ) { ProSHADE_internal_symmetry::printSymmetryGroup ( detectedCSymmetries.at(detIt), symPeaks, settings->verbose ); }
        
        //============================================ Save detected
        ProSHADE_internal_symmetry::saveAllCSymmetries ( detectedCSymmetries, symPeaks, &ret, settings->axisErrTolerance );
    }

    //================================================ Release memory
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( peaksAA.size() ); iter++ ) { delete[] allPeaks.at(iter); delete[] peaksAA.at(iter); }
    
    //================================================ Report completion
    ProSHADE_internal_symmetry::printSymmetryCompletion ( static_cast<proshade_unsign>( ret.size() ), settings->verbose );

    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function converts peaks ZXZ Euler anles to angle-axis representation for further processing.
 
    The only functionality here is taking a vector of Euler ZXZ angles and converting these though the rotation matrices to a vector
    of angle-axis representation of the same angles.
 
    \param[in] allPeaks A vector of pointers where Euler ZXZ representations of the peaks are saved.
    \param[in] verbose How loud the standard output of this run should be?
    \param[out] X A vector of pointers where angle-axis representations of the peaks will be saved.
 */
std::vector< proshade_double* > ProSHADE_internal_symmetry::getPeaksAngleAxisPositions ( std::vector< proshade_double* > allPeaks, proshade_unsign verbose )
{
    //================================================ Initialise variables
    std::vector< proshade_double* > ret;
    proshade_double* hlpP                             = NULL;
    proshade_double* rotMat                           = new proshade_double [9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    
    //================================================ For each peak
    for ( proshade_unsign peakIter = 0; peakIter < static_cast<proshade_unsign> ( allPeaks.size() ); peakIter++ )
    {
        //============================================ Convert Euler ZXZ angles to rotation matrix
        ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( allPeaks.at(peakIter)[0], allPeaks.at(peakIter)[1], allPeaks.at(peakIter)[2], rotMat );
        
        //============================================ Allocate pointer to results
        hlpP                                          = new proshade_double [5];
        ProSHADE_internal_misc::checkMemoryAllocation ( hlpP, __FILE__, __LINE__, __func__ );
        
        //============================================ Convert rotation matrix to Angle-axis reporesentation
        ProSHADE_internal_maths::getAxisAngleFromRotationMatrix ( rotMat, &hlpP[0], &hlpP[1], &hlpP[2], &hlpP[3] );
        hlpP[4]                                       = allPeaks.at(peakIter)[3];
        
        //============================================ Save results
        ProSHADE_internal_misc::addToDblPtrVector     ( &ret, hlpP );
    }
    
    //================================================ Release memory
    delete[] rotMat;
    
    //================================================ Report progress
    std::stringstream hlpSSP;
    hlpSSP << "Found " << ret.size() << " possible peaks.";
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, hlpSSP.str() );
    
    //================================================ Warning if no peaks!
    if ( ret.size() < 1 )
    {
        ProSHADE_internal_messages::printWarningMessage ( verbose, "!!! ProSHADE WARNING !!! Failed to detect any symmetries. There are no reasonable peaks in the self-rotation map. If you believe there should be some symmetry, you can try decreasing the resolution or changing the peak IQR threshold.", "WS00029" );
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function groups the peaks by height and returns the boundaries between such groups.
 
    This function allows for a list of peaks to be divided into multiple groups based on the peak heights, so that only the most
    confident values would be used for symmetry detection first.
 
    \param[in] allPeaks A vector of pointers where angle-axis representations of the peaks is saved.
    \param[in] smoothing Value determining how smooth the distribution of peaks should be made. Larger number means more groups.
    \param[out] X The boundaries for peak groups by height as determined by 1D grouping.
 */
std::vector< proshade_double > ProSHADE_internal_symmetry::findPeaksByHeightBoundaries ( std::vector< proshade_double* > allPeaks, proshade_double smoothing )
{
    //================================================ Initialise variables
    std::vector< proshade_double > boundaries;
    ProSHADE_internal_misc::addToDoubleVector         ( &boundaries, 0.0 );
    proshade_double peakContribution                  = 0.0;
    
    //================================================ Generate Probability Density function (PDF)
    std::vector< proshade_double > pdf;
    for ( proshade_double iter = 0.0; iter <= 1.0; iter += 0.01 )
    {
        //============================================ Initialise point
        peakContribution                              = 0.0;
        
        //============================================ Sum peak contributions
        for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( allPeaks.size() ); it++ )
        {
            peakContribution                         += ProSHADE_internal_maths::normalDistributionValue ( allPeaks.at(it)[4], smoothing, iter );
        }
        
        //============================================ Save result
        ProSHADE_internal_misc::addToDoubleVector     ( &pdf, peakContribution );
    }
    
    //================================================ Find boundaries
    proshade_double prev                              = pdf.at(0);
    for ( proshade_unsign iter = 1; iter < static_cast<proshade_unsign> ( pdf.size() - 1 ); iter ++ )
    {
        //============================================ Check for local minima
        if ( ( prev > pdf.at(iter) ) && ( pdf.at(iter+1) > pdf.at(iter) ) )
        {
            ProSHADE_internal_misc::addToDoubleVector ( &boundaries, iter * 0.01 );
        }
        
        //============================================ Prepare next iteration
        prev                                          = pdf.at(iter);
    }

    //================================================ Done
    return                                            ( boundaries );
    
}

/*! \brief This function searches the list of peaks for presence of cyclic symmetry.
 
    This function takes a set of peaks and a bunch of settings parameters and proceeds to search these peaks for containing any Cyclic
    (C) symmetries. It contains all the functionality including missing peaks searching and automatic possible fold detection including
    allowing for errors. It will finally save all the results in the vector of vectors it returns.
 
    \param[in] peaks A vector of pointers where angle-axis representations of the peaks is saved.
    \param[in] verbose How loud the standard output of this run should be?
    \param[in] band The bandwidth of these computations.
    \param[in] missPeakThres Threshold for the percentage of missing peaks there can be to warrant a full search for missing peaks.
    \param[in] axisErrTolerance Tolerance for symmetry axis identity.
    \param[in] axisErrToleranceDef Should the automatic axis tolerance decrease be applied?
    \param[in] dataObj The data object for which symmetry is being searched. This is only needed for missing peaks search, but needed nonetheless.
    \param[out] X Vector of vectors with first number being the detected fold and all remaining numbers being the indices of peaks forming the symmetry.
 */
std::vector< std::vector< proshade_unsign > > ProSHADE_internal_symmetry::findPeaksCSymmetry ( std::vector< proshade_double* >* peaks, proshade_signed verbose, proshade_unsign band, proshade_double missPeakThres, proshade_double axisErrTolerance, bool axisErrToleranceDef, ProSHADE_internal_data::ProSHADE_data* dataObj )
{
    //======================================== Initialise variables
    std::vector< std::vector< proshade_unsign > > ret;
    std::vector< proshade_double > triedAlready;
    std::vector< proshade_unsign > angsToTry, testedAlready;
    proshade_double angDist, angDivisionRemainder, angDivisionBasis, nextSymmetryError, nextPeakError = ( M_PI * 2.0 ) / ( static_cast<proshade_double> ( band ) * 2.0 );
    
    //================================================ Sanity check
    if ( peaks->size() < 1 ) { return ( ret ); }
    
    //================================================ Group peaks by axes
    std::vector< std::vector< proshade_unsign > > sameAxesGroups = ProSHADE_internal_symmetry::groupSameAxes ( *peaks, axisErrTolerance );
    
    //================================================ For each axis group
    for ( proshade_unsign grpIt = 0; grpIt < static_cast<proshade_unsign> ( sameAxesGroups.size() ); grpIt++ )
    {
        //============================================ Print axis group if need be
        ProSHADE_internal_symmetry::printSymmetryPeaks ( sameAxesGroups.at(grpIt), *peaks, verbose, grpIt );
 
        //============================================ While there are distances between rotation angles in the group
        triedAlready.clear                            ( );
        testedAlready.clear                           ( );
        while ( ProSHADE_internal_symmetry::smallestDistanceBetweenAngles ( sameAxesGroups.at(grpIt), *peaks, &triedAlready, &angDist ) )
        {
            //======================================== Check if testable fold value exists, otherwise test other distances
            angsToTry.clear                           ( );
            if ( !ProSHADE_internal_symmetry::determineFoldToTry ( angDist, &angDivisionBasis, &angDivisionRemainder, nextPeakError, &nextSymmetryError, &angsToTry ) ) { continue; }
            
            //======================================== If reasonable folds are found, test these for being complete symmetries
            ProSHADE_internal_symmetry::findSymmetryUsingFold ( dataObj, &angsToTry, &sameAxesGroups.at(grpIt), peaks, &ret, &testedAlready, axisErrTolerance, axisErrToleranceDef, missPeakThres, verbose );
        }
        
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function groups the peaks by their axes of rotation.
 
    This function takes the list of peaks so far detected and groups these by their axis or rotation, ignoring peaks
    with zero rotation angle. The return value is a vector of vectors of the groups and their members, but not re-organised
    list of peaks. This function also adds a zero angle peak to all peak groups (so that the zero angle peak has the same axis
    as all other group members for all groups).
 
    \param[in] peaks A vector of pointers where angle-axis representations of the peaks is saved.
    \param[in] errTolerance A value within which two axes are considered equal.
    \param[out] X A vector of peak groups with each group entry being a vector of groups member indices in the peaks vector.
 */
std::vector< std::vector< proshade_unsign > > ProSHADE_internal_symmetry::groupSameAxes ( std::vector< proshade_double* >& peaks, proshade_double errTolerance )
{
    //================================================ Initialise variables
    std::vector< std::vector< proshade_unsign > > ret;
    bool sameAxisFound                                = false;
    
    //================================================ Set all largest axis value to positive (this will make the 0,0,1 and 0,0,-1 axes the same)
    ProSHADE_internal_symmetry::giveOppositeAxesSameDirection ( peaks );
    
    //================================================ For each axis
    for ( proshade_unsign peakIter = 0; peakIter < static_cast<proshade_unsign> ( peaks.size() ); peakIter++ )
    {
        //============================================ Initialise variables for next peak
        sameAxisFound                                 = false;
        
        //============================================ Ignore zero angle peaks
        if ( ( peaks.at(peakIter)[3]-errTolerance <= 0.0 ) && ( peaks.at(peakIter)[3]+errTolerance > 0.0 ) ) { continue; }
        
        //============================================ Ignore zero axis peaks
        if ( ( ( peaks.at(peakIter)[0]-errTolerance <= 0.0 ) && ( peaks.at(peakIter)[0]+errTolerance > 0.0 ) ) && ( ( peaks.at(peakIter)[1]-errTolerance <= 0.0 ) && ( peaks.at(peakIter)[1]+errTolerance > 0.0 ) ) && ( ( peaks.at(peakIter)[2]-errTolerance <= 0.0 ) && ( peaks.at(peakIter)[2]+errTolerance > 0.0 ) ) ) { continue; }
        
        //============================================ Compare to all already detected axes groups
        for ( proshade_unsign sameAxisGrp = 0; sameAxisGrp < static_cast<proshade_unsign> ( ret.size() ); sameAxisGrp++ )
        {
            //======================================== and all their members
            for ( proshade_unsign sameAxis = 0; sameAxis < static_cast<proshade_unsign> ( ret.at(sameAxisGrp).size() ); sameAxis++ )
            {
                //==================================== Is this identical axis to the tested one?
                if ( ( ( (peaks.at(ret.at(sameAxisGrp).at(sameAxis))[0]-errTolerance) <= peaks.at(peakIter)[0] ) &&
                       ( (peaks.at(ret.at(sameAxisGrp).at(sameAxis))[0]+errTolerance) >  peaks.at(peakIter)[0] ) ) &&
                     ( ( (peaks.at(ret.at(sameAxisGrp).at(sameAxis))[1]-errTolerance) <= peaks.at(peakIter)[1] ) &&
                       ( (peaks.at(ret.at(sameAxisGrp).at(sameAxis))[1]+errTolerance) >  peaks.at(peakIter)[1] ) ) &&
                     ( ( (peaks.at(ret.at(sameAxisGrp).at(sameAxis))[2]-errTolerance) <= peaks.at(peakIter)[2] ) &&
                       ( (peaks.at(ret.at(sameAxisGrp).at(sameAxis))[2]+errTolerance) >  peaks.at(peakIter)[2] ) ) )
                {
                    sameAxisFound                     = true;
                    ProSHADE_internal_misc::addToUnsignVector ( &ret.at(sameAxisGrp), peakIter );
                    break;
                }
            }
        }
        
        //============================================ If same axis was found, do nothing
        if ( sameAxisFound ) { continue; }
        
        //============================================ No similar axis was found
        std::vector<proshade_unsign> hlpVec;
        ProSHADE_internal_misc::addToUnsignVector     ( &hlpVec, peakIter );
        ProSHADE_internal_misc::addToUnsignVectorVector ( &ret, hlpVec );
    }
    
    //================================================ Add zero peak to all axes
    ProSHADE_internal_symmetry::addZeroPeakToGroups   ( ret, peaks );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function modifiest the axes so that the highest vector element is always positive.
 
    This function modifies the angle-axis representation of the peak positions so that the leargest dimmension of the rotation axis would
    be positive. This is important in order to make sure that the AA representations [0,0,1;3.14] and [0,0,-1;-3.14] are equal and not considered
    as completely different.
 
    \param[in] peaks A vector of pointers where angle-axis representations of the peaks is saved.
 */
void ProSHADE_internal_symmetry::giveOppositeAxesSameDirection ( std::vector< proshade_double* > peaks )
{
    //================================================ Apply to all peaks
    for ( proshade_unsign i = 0; i < static_cast<proshade_unsign> ( peaks.size() ); i++ )
    {
        if ( ( ( std::max ( std::abs ( peaks.at(i)[0] ), std::max( std::abs ( peaks.at(i)[1] ), std::abs ( peaks.at(i)[2] ) ) ) == std::abs ( peaks.at(i)[0] ) ) && ( peaks.at(i)[0] < 0.0 ) ) ||
             ( ( std::max ( std::abs ( peaks.at(i)[0] ), std::max( std::abs ( peaks.at(i)[1] ), std::abs ( peaks.at(i)[2] ) ) ) == std::abs ( peaks.at(i)[1] ) ) && ( peaks.at(i)[1] < 0.0 ) ) ||
             ( ( std::max ( std::abs ( peaks.at(i)[0] ), std::max( std::abs ( peaks.at(i)[1] ), std::abs ( peaks.at(i)[2] ) ) ) == std::abs ( peaks.at(i)[2] ) ) && ( peaks.at(i)[2] < 0.0 ) ) )
        {
            peaks.at(i)[0]                           *= -1.0;
            peaks.at(i)[1]                           *= -1.0;
            peaks.at(i)[2]                           *= -1.0;
            peaks.at(i)[3]                           *= -1.0;
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function simply prints the symmetry axis group supplied in the first parameter from the second parameter values.
 
    \param[in] grp A single symmetry axis group indices to be printed.
    \param[in] peaks The vector of all peaks from which the indices are drawn.
    \param[in] verbose How loud the run should be and therefore if anything should be printed at all.
 */
void ProSHADE_internal_symmetry::printSymmetryPeaks ( std::vector< proshade_unsign > grp, std::vector< proshade_double* > peaks, proshade_signed verbose, proshade_unsign groupNo )
{
    //================================================ Symmetry group output header
    std::stringstream hlpSS;
    hlpSS << "Symmetry axis group " << groupNo;
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 6, hlpSS.str() );
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 6, "Peak index\t\tx\t y\t z\tAngle\tPeak heiht" );
    
    //================================================ Print the symmetry group
    for ( proshade_unsign axIt = 0; axIt < static_cast<proshade_unsign> ( grp.size() ); axIt++ )
    {
        std::stringstream SS;
        SS << "    " << axIt << "\t      " << static_cast<int>( peaks.at(grp.at(axIt))[0] * 100.0 ) / 100.0 << "\t" << static_cast<int>( peaks.at(grp.at(axIt))[1] * 100.0 ) / 100.0 << "\t" << static_cast<int>( peaks.at(grp.at(axIt))[2] * 100.0 ) / 100.0 << "\t" << static_cast<int>( peaks.at(grp.at(axIt))[3] * 100.0 ) / 100.0 << "\t" << static_cast<int>( peaks.at(grp.at(axIt))[4] * 100.0 ) / 100.0;
        ProSHADE_internal_messages::printProgressMessage ( verbose, 6, SS.str() );
    }
    
    //================================================ Done
    return ;

}

/*! \brief This function finds the smallest distance between the rotation angles within a group.
 
    This function is used to control the while loop in the findPeaksCSymmetry() function. It has two outputs, the standard returned value
    is a boolean stating whether a new distance between group rotation angles was found; the second output is the distance itself, which is
    saved in the dist variable.
 
    \param[in] grp A single symmetry axis group indices to be printed.
    \param[in] peaks The vector of all peaks from which the indices are drawn.
    \param[in] tried A vector of doubles holding the already tried distances and group combinations, so that they would not be tried again.
    \param[in] dist A pointer to the variable where the smallest distance (if found) will be saved.
    \param[out] X Bool whether a new distance was found.
 */
bool ProSHADE_internal_symmetry::smallestDistanceBetweenAngles ( std::vector< proshade_unsign > grp, std::vector< proshade_double* > peaks, std::vector< proshade_double >* tried, proshade_double* dist )
{
    //================================================ Initialise variables
    bool ret                                          = false;
    bool skip                                         = false;
    proshade_unsign g1 = 0, g2 = 0;
   *dist                                              = 999.9;
    
    //================================================ For each group pair
    for ( proshade_unsign gr1It = 0; gr1It < static_cast<proshade_unsign> ( grp.size() ); gr1It++ )
    {
        for ( proshade_unsign gr2It = 1; gr2It < static_cast<proshade_unsign> ( grp.size() ); gr2It++ )
        {
            //======================================== Unique pairs only
            if ( gr1It >= gr2It ) { continue; }
            
            //======================================== Have we tried this already?
            skip                                      = false;
            for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( tried->size() ); iter += 3 )
            {
                //==================================== Avoid already tested combinations
                if ( ( gr2It == tried->at( iter + 1 ) ) && ( gr1It == tried->at( iter ) ) ) { skip = true; }
                
                //==================================== Also avoid distances very close to already tested  distances (no problem until approx C700)
                if ( !skip &&
                     ( ( std::abs( std::abs ( peaks.at(grp.at(gr1It))[3] ) - std::abs ( peaks.at(grp.at(gr2It))[3] ) ) - 0.01 ) < tried->at( iter + 2 ) ) &&
                     ( ( std::abs( std::abs ( peaks.at(grp.at(gr1It))[3] ) - std::abs ( peaks.at(grp.at(gr2It))[3] ) ) + 0.01 ) > tried->at( iter + 2 ) ) )
                {
                    skip                              = true;
                }
            }
            if ( skip ) { continue; }
            
            //======================================== Is this the smallest distance?
            if ( std::abs( std::abs ( peaks.at(grp.at(gr1It))[3] ) - std::abs ( peaks.at(grp.at(gr2It))[3] ) ) < (*dist) )
            {
                //==================================== Avoid very small angle distances as they would just take time (the hardcoded value would only be a problem for C700 and larger symmetries...
                if ( std::abs( std::abs ( peaks.at(grp.at(gr1It))[3] ) - std::abs ( peaks.at(grp.at(gr2It))[3] ) ) > 0.01 )
                {
                    g1                                = gr1It;
                    g2                                = gr2It;
                   *dist                              = std::abs( std::abs ( peaks.at(grp.at(gr1It))[3] ) - std::abs ( peaks.at(grp.at(gr2It))[3] ) );
                }
            }
        }
    }
    
    //================================================ If new dist found, save to tried and return success
    if ( *dist != 999.9 )
    {
        ret                                           = true;
        ProSHADE_internal_misc::addToDoubleVector     ( tried, g1 );
        ProSHADE_internal_misc::addToDoubleVector     ( tried, g2 );
        ProSHADE_internal_misc::addToDoubleVector     ( tried, *dist );
    }
    
    //================================================ Done
    return                                            ( ret );

}

/*! \brief This function takes the peak groups and adds zero peak to each of them.
 
    This function takes all of the detected peak axis groups and the list of peaks. It then proceeds to add a single peak per a group to
    the peaks list; this newly added peak has the same axis as the group, but zero angle. The function also adds the index of this new
    peak to the peak group, so that the group now has a new member, a peak with zero angle and the same axis.
 
    \param[in] grpsVec A list of all symmetry axis groups.
    \param[in] peaks The vector of all peaks from which the group indices are to be drawn.
 */
void ProSHADE_internal_symmetry::addZeroPeakToGroups ( std::vector< std::vector< proshade_unsign > >& grpsVec, std::vector< proshade_double* >& peaks )
{
    //================================================ Do your job
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( grpsVec.size() ); iter++ )
    {
        proshade_double* hlpP                         = new proshade_double [5];
        ProSHADE_internal_misc::checkMemoryAllocation ( hlpP, __FILE__, __LINE__, __func__ );
        hlpP[0]                                       = peaks.at(grpsVec.at(iter).at(0))[0];
        hlpP[1]                                       = peaks.at(grpsVec.at(iter).at(0))[1];
        hlpP[2]                                       = peaks.at(grpsVec.at(iter).at(0))[2];
        hlpP[3]                                       = 0.0;
        hlpP[4]                                       = peaks.at(grpsVec.at(iter).at(0))[4];
        ProSHADE_internal_misc::addToUnsignVector     ( &grpsVec.at(iter), static_cast<proshade_unsign> ( peaks.size() ) );
        ProSHADE_internal_misc::addToDblPtrVector     ( &peaks, hlpP );
    }
    
    //================================================ Done
    return ;

}

/*! \brief This function determines the symmetry fold to be searched for.
 
    This function detects which fold would belong to the rotation angle distance supplied. This is done by finding the division basis
    for the simple 2pi/dist equation and minimising the remainder. The function then checks whether the remainder is smaller than a
    threshold and whether the error on fold detection is not close to fold+1 value in terms of peak misplacement in the map - if it is,
    then surrounding fold values are also added to be tested. Finally, the function returns boolean value stating whether at least one
    testable fold value passed the checks.
 
    \param[in] dist The distance between rotation angles that should form the symmetry group.
    \param[in] divBasis Pointer to where to save the basis of the division 2pi/dist.
    \param[in] divRem Pointer to where to save the remainder of the division 2pi/dist.
    \param[in] peakErr The error in radians which would be the result of misplacing a peak by single map index.
    \param[in] symmErr Pointer to where to save the error which would be caused by mis-predicting the fold by 1.
    \param[in] angsToTry A vector where all the suggested fold values to be tested are saved.
    \param[out] X Boolean value whether at least single testable fold value was found
 */
bool ProSHADE_internal_symmetry::determineFoldToTry ( proshade_double dist, proshade_double* divBasis, proshade_double* divRem, proshade_double peakErr, proshade_double* symmErr, std::vector< proshade_unsign >* angsToTry )
{
    //================================================ Initialise variables
    bool ret                                          = false;
    
    //================================================ Find the basis and remainder of the 2pi/dist equation
   *divRem                                            = std::modf ( static_cast<double> ( ( 2.0 * M_PI ) / std::abs ( dist ) ), divBasis );
    
    //================================================ If the remainder would be smaller for larger basis, so this basis
    if ( *divRem > 0.5 )
    {
       *divRem                                       -= 1.0;
       *divBasis                                     += 1.0;
    }
    
    //================================================ Determine errors on peaks and on folds
   *symmErr                                           = ( M_PI * 2.0 / *divBasis ) - ( M_PI * 2.0 / ( *divBasis + 1.0 ) );
    proshade_double angTolerance                      = ( peakErr / *symmErr  );
    
    //================================================ Is remainder small enough?
    if ( ( *divRem < ( 0.0 + angTolerance ) ) && ( *divRem > ( 0.0 - angTolerance ) ) )
    {
        //============================================ Are we sure about the fold determination accuracy
        proshade_signed angTolRound                   = std::min ( ProSHADE_internal_mapManip::myRound ( angTolerance ), static_cast<proshade_signed> ( 10 ) );
        for ( proshade_signed iter = -angTolRound; iter <= angTolRound; iter++ )
        {
            ProSHADE_internal_misc::addToUnsignVector ( angsToTry, static_cast<proshade_unsign> ( std::max ( *divBasis + iter, 2.0 ) ) );
        }
    }
    
    //================================================ Return indication of whether testable fold value(s) was found.
    if ( angsToTry->size() == 0 ) { ret = false; }
    else { ret = true; }
    
    //================================================ Done
    return                                    ( ret );

}

/*! \brief This function computes the expected peak rotations for given fold.
 
    This function computes the expected peak rotation angle values for the peak range between -180 to +180 degrees plus one
    distance on both sides for a good measure. The resulting values are then saved to the second parameter vector.
 
    \param[in] fold The fold for which peak rotation angles should be predicted.
    \param[in] expAngs A vector where the expected peak rotation values will be saved to.
 */
void ProSHADE_internal_symmetry::findExpectedPeakRotations ( proshade_unsign fold, std::vector< proshade_double >* expAngs )
{
    //================================================ Initialise variables
    proshade_double groupAngle                        = ( 2.0 * M_PI ) / static_cast<proshade_double> ( fold );
    
    //================================================ Generate expected angles
    for ( proshade_signed iter = static_cast<proshade_signed> ( -( static_cast<proshade_double> ( fold ) / 2.0 + 1.0) ); iter <= static_cast<proshade_signed> ( static_cast<proshade_double> ( fold )/2.0 + 1.0 ); iter++ )
    {
        ProSHADE_internal_misc::addToDoubleVector     ( expAngs, iter * groupAngle );
    }
    
    //================================================ Done
    return ;

}

/*! \brief This function computes the expected peak rotations for given fold.
 
    This function compares the expected and the detected peak rotation angle values to check if the complete C symmetry is found within this peak axis group.
    It also saves the indices of the matched and missing peaks and returns the number of consecutive mathes.
 
    \param[in] grp A single symmetry axis group indices to be processed.
    \param[in] peaks The vector of all peaks from which the indices are drawn.
    \param[in] expAngs A vector where the expected peak rotation values are saved.
    \param[in] matchedAngs A vector where the indices of matched peaks will be saved.
    \param[in] missingAngs A vector where the indices of missing peaks will be saved.
    \param[in] axisTol The tolerance for matching the expected and found peak rotation angles.
    \param[out] X An integer with the longest consecutive streak of matched values.
 */
proshade_unsign ProSHADE_internal_symmetry::checkExpectedAgainstFound ( std::vector< proshade_unsign > grp, std::vector< proshade_double* > peaks, std::vector< proshade_double >* expAngs, std::vector< proshade_unsign >* matchedAngs, std::vector< proshade_unsign >* missingAngs, proshade_double axisTol )
{
    //================================================ Initialise variables
    proshade_unsign ret                               = 0;
    proshade_unsign retHlp                            = 0;
    proshade_double groupAngle                        = expAngs->at(1) - expAngs->at(0);
    bool matchedThisPeak                              = false;
    bool noDoubleMatches                              = false;
    std::vector < proshade_unsign > matchedAlready;
    
    //================================================ For each expected peak rotation angle value
    for ( proshade_unsign expAngIt = 0; expAngIt < static_cast<proshade_unsign> ( expAngs->size() ); expAngIt++ )
    {
        //============================================ For each peak in the group
        matchedThisPeak                               = false;
        for ( proshade_unsign peakIt = 0; peakIt < static_cast<proshade_unsign> ( grp.size() ); peakIt++ )
        {
            if ( ( expAngs->at(expAngIt) < ( peaks.at(grp.at(peakIt))[3] + axisTol ) ) &&
                 ( expAngs->at(expAngIt) > ( peaks.at(grp.at(peakIt))[3] - axisTol ) ) )
            {
                noDoubleMatches                       = false;
                for ( proshade_unsign ndm = 0; ndm < static_cast<proshade_unsign> ( matchedAlready.size() ); ndm++ )
                {
                    if ( matchedAlready.at(ndm) == grp.at(peakIt) ) { noDoubleMatches = true; break; }
                }
                
                if ( !noDoubleMatches )
                {
                    ProSHADE_internal_misc::addToUnsignVector ( matchedAngs, grp.at(peakIt) );
                    ProSHADE_internal_misc::addToUnsignVector ( &matchedAlready, grp.at(peakIt) );
                    matchedThisPeak                   = true;
                    break;
                }
            }
        }
        
        //============================================ If not matched, add to missing
        if ( !matchedThisPeak )
        {
            ProSHADE_internal_misc::addToUnsignVector ( missingAngs, expAngIt );
        }
    }
    
    //================================================ Find the number of consecutive matches
    if ( matchedAngs->size () > 1 )
    {
        for ( proshade_unsign iter = 1; iter < static_cast<unsigned int> ( matchedAngs->size () ); iter++ )
        {
            if ( ( ( peaks.at(matchedAngs->at(iter-1))[3] + groupAngle ) < ( peaks.at(matchedAngs->at(iter))[3] + axisTol ) ) &&
                 ( ( peaks.at(matchedAngs->at(iter-1))[3] + groupAngle ) > ( peaks.at(matchedAngs->at(iter))[3] - axisTol ) ) )
            {
                retHlp                       += 1;
            }
            else
            {
                retHlp                        = 0;
            }
            if ( retHlp > ret ) { ret = retHlp; }
        }
    }

    //================================================ Done
    return                                            ( ret + 1 ); // This is because the count of matches is the count of intervals between numbers, so +1 to get the count of matched numbers.

}

/*! \brief This function checks for the high of the correlation for particular rotation angle and axis.
 
    This is the core of missing peaks procedure. This function takes the angle-axis representation of the sought after peak/rotation
    and searches the data objects (respectivelly its inverse SO(3) Fourier Transform map) for the highest point conforming to these
    specifications. It then returns the highest value found, so that it could be decided whether the symmetry search has been successfully
    completed or whether the symmetry was not found.
 
    \param[in] dataObj The data object for which symmetry is being searched.
    \param[in] x The x-axis element of the searched for rotation angle-axis representation.
    \param[in] y The y-axis element of the searched for rotation angle-axis representation.
    \param[in] z The z-axis element of the searched for rotation angle-axis representation.
    \param[in] angle The angle element of the searched for rotation angle-axis representation.
    \param[in] heightThres The required self-rotation map height for this rotation.
    \param[in] axTol The tolerance on axis matching when searching for the rotation.
    \param[out] X The height of highest matching map point.
 */
proshade_double ProSHADE_internal_symmetry::checkForMissingPeak ( ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_double x, proshade_double y, proshade_double z, proshade_double angle, proshade_double heightThres, proshade_double axTol )
{
    //================================================ Initialise variables
    proshade_double ret                               = 0.0;
    proshade_unsign arrIndex                          = 0;
    proshade_double* rotMat                           = new proshade_double [9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    proshade_double pointHeight, euA, euB, euG, xPk, yPk, zPk, anglPk;
    
    //================================================ Search the self-rotation map
    for ( proshade_unsign xIt = 0; xIt < ( dataObj->getMaxBand() * 2 ); xIt++ )
    {
        for ( proshade_unsign yIt = 0; yIt < ( dataObj->getMaxBand() * 2 ); yIt++ )
        {
            for ( proshade_unsign zIt = 0; zIt < ( dataObj->getMaxBand() * 2 ); zIt++ )
            {
                //==================================== Get height and check against threshold
                arrIndex                              = zIt  + ( dataObj->getMaxBand() * 2 ) * ( yIt  + ( dataObj->getMaxBand() * 2 ) * xIt );
                pointHeight                           = pow( dataObj->getInvSO3Coeffs()[arrIndex][0], 2.0 ) + pow( dataObj->getInvSO3Coeffs()[arrIndex][1], 2.0 );
                if ( pointHeight < heightThres ) { continue; }
                
                //==================================== Get angle-axis values
                ProSHADE_internal_maths::getEulerZXZFromSOFTPosition ( dataObj->getMaxBand(), static_cast<proshade_signed> ( xIt ),
                                                                       static_cast<proshade_signed> ( yIt ), static_cast<proshade_signed> ( zIt ),
                                                                       &euA, &euB, &euG );
                ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( euA, euB, euG, rotMat );
                ProSHADE_internal_maths::getAxisAngleFromRotationMatrix ( rotMat, &xPk, &yPk, &zPk, &anglPk );
                
                //==================================== Check for matching angle
                if ( ( ( std::abs( anglPk ) - axTol ) < std::abs ( angle ) ) && ( ( std::abs( anglPk ) + axTol ) > std::abs ( angle ) ) )
                {
                    //================================ Make sure vector direction is the same
                    if ( ( ( std::max( std::abs( xPk ), std::max( std::abs( yPk ), std::abs( zPk ) ) ) == std::abs( xPk ) ) && ( xPk < 0.0 ) ) ||
                         ( ( std::max( std::abs( xPk ), std::max( std::abs( yPk ), std::abs( zPk ) ) ) == std::abs( yPk ) ) && ( yPk < 0.0 ) ) ||
                         ( ( std::max( std::abs( xPk ), std::max( std::abs( yPk ), std::abs( zPk ) ) ) == std::abs( zPk ) ) && ( zPk < 0.0 ) ) )
                    {
                        xPk                          *= -1.0;
                        yPk                          *= -1.0;
                        zPk                          *= -1.0;
                        anglPk                       *= -1.0;
                    }
                    
                    //================================ Compare axis elements
                    if ( ( ( ( xPk - axTol ) < x ) && ( ( xPk + axTol ) > x ) ) &&
                         ( ( ( yPk - axTol ) < y ) && ( ( yPk + axTol ) > y ) ) &&
                         ( ( ( zPk - axTol ) < z ) && ( ( zPk + axTol ) > z ) ) )
                    {
                        if ( ret < pointHeight ) { ret = pointHeight; }
                    }
                }
            }
        }
    }
    
    //================================================ Done
    return                                            ( ret );

}

/*! \brief This function saves a detected symmetry for reporting to the user.
 
    This function simply saves the supplied group members and fold value to the main output vector of vectors (also supplied). It makes sure the saving
    format (fold first, then all symmetry peak indices) is upheld.
 
    \param[in] fold This is the fold value of the detected C symmetry.
    \param[in] matchedPeaks A vector containing the indices of all peaks forming this symmetry.
    \param[in] ret The vector of vectors to be returned by findPeaksCSymmetry() and containing all detected symmetries (to which we are saving here).
    \param[in] verbose How loud the standard output of this run should be?
 */
void ProSHADE_internal_symmetry::saveDetectedCSymmetry ( proshade_unsign fold, std::vector< proshade_unsign >* matchedPeaks, std::vector< std::vector< proshade_unsign > >* ret, proshade_signed verbose )
{
    //================================================ Save fold as first vector value
    std::vector< proshade_unsign > hlpVec;
    ProSHADE_internal_misc::addToUnsignVector         ( &hlpVec, fold );
    
    //================================================ and follow it with indices of all symmetry forming peaks
    for ( proshade_unsign pIt = 0; pIt < static_cast<proshade_unsign> ( matchedPeaks->size() ); pIt++ )
    {
        ProSHADE_internal_misc::addToUnsignVector     ( &hlpVec, matchedPeaks->at(pIt) );
    }
    ProSHADE_internal_misc::addToUnsignVectorVector   ( ret, hlpVec );
    
    //================================================ Report finding symmetry
    std::stringstream hlpS;
    hlpS << "Found symmetry C" << fold;
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 5, hlpS.str() );
    
    //================================================ Done
    return ;

}

/*! \brief This function does the complete missing peak searching and filling in the missing peaks.
 
    This function does the complete work on missing peaks and detection of symmetries affected by them. It firstly decides on the threshold
    for a missing peak height and it then proceeds to check all missing peaks for being in the inverse SO(3) FT map. Any detected peaks will
    be saved to all appropriate variables (as supplied) and finally, if all missing peaks were sucessfully found, it will return true, otherwise
    false.
 
    \param[in] dataObj The data object for which symmetry is being searched.
    \param[in] fold This is the fold value of the detected C symmetry.
    \param[in] grp Vector with the indices of members of this symmetry axis group.
    \param[in] peaks A vector of pointers where angle-axis representations of the peaks is saved.
    \param[in] missingPeaks Vector with the indices of missing rotation angles (indices are from the expected peaks vector, not peaks vector!).
    \param[in] expectedAngles Vector with the expected rotation angle values.
    \param[in] axErrTolerance The allowed error on matching axes.
    \param[in] verbose How loud the standard output of this run should be?
    \param[out] X Was the missing symmetry part completion successfull?
 */
bool ProSHADE_internal_symmetry::completeMissingCSymmetry ( ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_unsign fold, std::vector< proshade_unsign >* grp, std::vector< proshade_double* >* peaks, std::vector< proshade_unsign >* missingPeaks, std::vector< proshade_double >* expectedAngles, std::vector< proshade_unsign >* matchedPeaks, proshade_double axErrTolerance, proshade_unsign verbose )
{
    //================================================ Initialise variables
    bool ret                                          = true;
    
    //================================================ Report searching for missing peaks
    std::stringstream hlpSSP;
    hlpSSP << "Searching for missing peaks for symmetry C" << fold;
    ProSHADE_internal_messages::printProgressMessage ( verbose, 4, hlpSSP.str() );
    
    //================================================ Height threshold for missing peak
    proshade_double heightThreshold                   = 0.0;
    for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( grp->size() ); grIt++ ) { heightThreshold += peaks->at(grp->at(grIt))[4]; }
    heightThreshold                                  /= static_cast<proshade_double> ( grp->size() );
    heightThreshold                                  *= 0.5;
    
    //================================================ For each missing value
    for ( proshade_unsign misPkIt = 0; misPkIt < static_cast<proshade_unsign> ( missingPeaks->size() ); misPkIt++ )
    {
        //============================================ Ignore the extra values in the expected values
        if ( expectedAngles->at(missingPeaks->at(misPkIt)) >   M_PI ) { continue; }
        if ( expectedAngles->at(missingPeaks->at(misPkIt)) <  -M_PI ) { continue; }
        
        //============================================ Search for the missing peak
        proshade_double misHeight = ProSHADE_internal_symmetry::checkForMissingPeak ( dataObj, peaks->at(grp->at(0))[0], peaks->at(grp->at(0))[1], peaks->at(grp->at(0))[2], expectedAngles->at(missingPeaks->at(misPkIt)), heightThreshold, axErrTolerance );
        if ( misHeight != 0.0 )
        {
            //======================================== Missing peak detected - save it to the group
            proshade_double* hlpP                     = new proshade_double [5];
            ProSHADE_internal_misc::checkMemoryAllocation ( hlpP, __FILE__, __LINE__, __func__ );
            hlpP[0]                                   = peaks->at(grp->at(0))[0];
            hlpP[1]                                   = peaks->at(grp->at(0))[1];
            hlpP[2]                                   = peaks->at(grp->at(0))[2];
            hlpP[3]                                   = expectedAngles->at(missingPeaks->at(misPkIt));
            hlpP[4]                                   = misHeight;
            ProSHADE_internal_misc::addToUnsignVector ( matchedPeaks, static_cast<proshade_unsign> ( peaks->size() ) );
            ProSHADE_internal_misc::addToDblPtrVector ( peaks, hlpP );
        }
        else
        {
            ret                                       = false;
        }
    }
    
    //================================================ Done
    return                                            ( ret );

}

/*! \brief This function tests all supplied folds for being supported by the peaks (i.e. and being complete present symmetry).
 
    This function takes all the possible folds which could be in the set of peaks and checks if these are indeed full symmetries,
    or whether these were random.
 
    \param[in] dataObj The data object for which symmetry is being searched.
    \param[in] angsToTry This vector contains all the folds that should be attempted.
    \param[in] grp Vector with the indices of members of this symmetry axis group.
    \param[in] peaks A vector of pointers where angle-axis representations of the peaks is saved.
    \param[in] ret The final variable holding all results (i.e. detected symmetries).
    \param[in] testedAlready A vector in which the already tested folds for this symmetry axis are saved.
    \param[in] axErrTolerance The allowed error on matching axes.
    \param[in] axErrToleranceDefault Should the axErrTolerance be decreased with increasing fold?
    \param[in] missPeakThres Threshold for the percentage of missing peaks there can be to warrant a full search for missing peaks.
    \param[in] verbose How loud the standard output of this run should be?
 */
void ProSHADE_internal_symmetry::findSymmetryUsingFold ( ProSHADE_internal_data::ProSHADE_data* dataObj, std::vector< proshade_unsign >* angsToTry, std::vector< proshade_unsign >* grp, std::vector< proshade_double* >* peaks, std::vector< std::vector< proshade_unsign > >* ret, std::vector< proshade_unsign >* testedAlready, proshade_double axErrTolerance, bool axErrToleranceDefault, proshade_double missPeakThres, proshade_unsign verbose )
{
    //================================================ Initialise variables
    bool skipFold                                     = false;
    std::vector< proshade_unsign > matchedPeaks, missingPeaks;
    std::vector< proshade_double > expectedAngles;
    proshade_double axErrToleranceOrig                = axErrTolerance;
    
    //================================================ Testing folds for being supported by peaks
    for ( proshade_unsign fIt = 0; fIt < static_cast<proshade_unsign> ( angsToTry->size() ); fIt++ )
    {
        //============================================ Was this fold already found?
        skipFold                                      = false;
        for ( proshade_unsign ftIt = 0; ftIt < static_cast<proshade_unsign> ( testedAlready->size() ); ftIt++ ) { if ( testedAlready->at(ftIt) == angsToTry->at(fIt) ) { skipFold = true; } }
        if ( skipFold ) { continue; }
        else { ProSHADE_internal_misc::addToUnsignVector( testedAlready, angsToTry->at(fIt) ); }
        
        //============================================ Set axis tolerance based on fold (if required)
        if ( axErrToleranceDefault )
        {
            axErrTolerance                            = std::max ( std::min ( ( ( 360.0 / static_cast<double> ( angsToTry->at(fIt) ) ) -
                                                                                ( 360.0 / static_cast<double> ( angsToTry->at(fIt) + 1 ) ) ) /
                                                                             acos ( axErrTolerance ) * axErrTolerance, axErrTolerance ), 0.02 );
        }
        
        //============================================ Find expected peak rotation angles
        expectedAngles.clear                          ( );
        ProSHADE_internal_symmetry::findExpectedPeakRotations ( angsToTry->at(fIt), &expectedAngles );
        
        //============================================ Compare group to expected angles
        matchedPeaks.clear                            ( );
        missingPeaks.clear                            ( );
        proshade_unsign consecMatches                 = ProSHADE_internal_symmetry::checkExpectedAgainstFound ( *grp, *peaks, &expectedAngles,
                                                                                                                &matchedPeaks, &missingPeaks, axErrTolerance );
    
        //============================================ If enough consecutive matches, symmetry was found. Save it
        if ( consecMatches >= angsToTry->at(fIt) )
        {
            ProSHADE_internal_symmetry::saveDetectedCSymmetry ( angsToTry->at(fIt), &matchedPeaks, ret, verbose );
        }
        else
        {
            if ( ( static_cast<proshade_double> ( matchedPeaks.size() ) / static_cast<proshade_double> ( angsToTry->at(fIt) ) ) >= ( 1.0 - missPeakThres ) )
            {
                //==================================== Attempt completing the symmetry  using missing peaks
                if ( ProSHADE_internal_symmetry::completeMissingCSymmetry ( dataObj, angsToTry->at(fIt), grp, peaks, &missingPeaks,
                                                                           &expectedAngles, &matchedPeaks, axErrTolerance, verbose ) )
                {
                    ProSHADE_internal_symmetry::saveDetectedCSymmetry ( angsToTry->at(fIt), &matchedPeaks, ret, verbose );
                }
            }
            else
            {
                //=================================== Symmetry not detected
                continue;
            }
        }
    }
    
    //=============================================== Done
    return ;

}

/*! \brief This function simply prints the detected symmetry and all its supporting peaks.
 
    \param[in] grp A single symmetry axis group indices to be printed.
    \param[in] peaks The vector of all peaks from which the indices are drawn.
    \param[in] verbose How loud the run should be and therefore if anything should be printed at all.
 */
void ProSHADE_internal_symmetry::printSymmetryGroup ( std::vector< proshade_unsign > grp, std::vector< proshade_double* > peaks, proshade_signed verbose )
{
    //================================================ Detected symmetry table header
    std::stringstream ss;
    ss << "Detected C" << grp.at(0) << " symmetry with following peaks:";
    ProSHADE_internal_messages::printProgressMessage ( verbose, 5, ss.str() );
    ProSHADE_internal_messages::printProgressMessage ( verbose, 5, "\tx\t y\t z\tAngle\tPeak height" );
    
    //================================================ Now print all supporting peaks
    for ( proshade_unsign pkIt = 1; pkIt < static_cast<proshade_unsign> ( grp.size() ); pkIt++ )
    {
        std::stringstream SS;
        SS << "  " << static_cast<int>( peaks.at(grp.at(pkIt))[0] * 100.0 ) / 100.0 << "\t" << static_cast<int>( peaks.at(grp.at(pkIt))[1] * 100.0 ) / 100.0 << "\t" << static_cast<int>( peaks.at(grp.at(pkIt))[2] * 100.0 ) / 100.0 << "\t" << static_cast<int>( peaks.at(grp.at(pkIt))[3] * 100.0 ) / 100.0 << "\t" << static_cast<int>( peaks.at(grp.at(pkIt))[4] * 100.0 ) / 100.0;
        ProSHADE_internal_messages::printProgressMessage ( verbose, 5, SS.str() );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function simply prints the summary and warnings for cyclic symmetries detection completion.
 
    \param[in] noSyms The number of symmetries that were detected.
    \param[in] verbose How loud the run should be and therefore if anything should be printed at all.
 */
void ProSHADE_internal_symmetry::printSymmetryCompletion ( proshade_unsign noSyms, proshade_unsign verbose )
{
    //================================================ Report completion of symmetry detection
    std::stringstream ss;
    ss << "Detected " << noSyms << " Cyclic symmetries.";
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, ss.str() );
    
    //================================================ If no symmetries were found, print warning
    if ( noSyms < 1 )
    {
        ProSHADE_internal_messages::printWarningMessage ( verbose, "!!! ProSHADE WARNING !!! Failed to detect any symmetries. If you believe there should be one, you can try decreasing the resolution or checking that the map is centred on the centry of symmetry (or use map centering option in ProSHADE).", "WS00030" );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes the detected symmetries indices and peaks and saves these in the main cyclic symmetries detection output format.
 
    This function uses the indices of peaks forming a detected symmetry along with the peak values corresponding to these indices in order to
    compute the symmetry description - that is the fold, average x, y and z-axis elements, angle (2pi/fold) and the average peak height. With
    all this computed for each detected symmetry, it saves these as double arrays to the output vector of double arrays for further processing.
    The function also does not save redundant symmetries.
 
    \param[in] detected This is a vector of vectors with the indices of detected symmetry peaks.
    \param[in] peaks These are the peaks and their values which come together to form the detected symmetry.
    \param[in] ret This is the variable where the results will be saved. It is a vector of double[6] arrays with the following meaning: [0] = fold, [1] = x-axis, [2] = y-axis, [3] = z-axis, [4] = angle, [5] = average peak height.
    \param[in] axErr The tolerance on axis matching.
 */
void ProSHADE_internal_symmetry::saveAllCSymmetries ( std::vector< std::vector< proshade_unsign > > detected, std::vector< proshade_double* > peaks, std::vector< proshade_double* >* ret, proshade_double axErr )
{
    //================================================ Initialise variables
    proshade_double sumX, sumY, sumZ, sumH;
    
    //================================================ Start saving
    for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( detected.size() ); symIt++ )
    {
        //============================================ Allocate the memory
        proshade_double* hlpP                         = new proshade_double [6];
        ProSHADE_internal_misc::checkMemoryAllocation ( hlpP, __FILE__, __LINE__, __func__ );
        
        //============================================ Set obvious values
        hlpP[0]                                       = static_cast<proshade_double> ( detected.at(symIt).at(0) );
        hlpP[4]                                       = static_cast<proshade_double> ( ( 2.0 * M_PI ) / hlpP[0] );
        
        //============================================ Compute peak averages for rest
        sumX = 0.0; sumY = 0.0; sumZ = 0.0; sumH = 0.0;
        for ( proshade_unsign pkIt = 1; pkIt < static_cast<proshade_unsign> ( detected.at(symIt).size() ); pkIt++ )
        {
            sumX                                     += peaks.at(detected.at(symIt).at(pkIt))[0];
            sumY                                     += peaks.at(detected.at(symIt).at(pkIt))[1];
            sumZ                                     += peaks.at(detected.at(symIt).at(pkIt))[2];
            sumH                                     += peaks.at(detected.at(symIt).at(pkIt))[4];
        }
        sumX                                         /= static_cast<proshade_double> ( detected.at(symIt).size() - 1 );
        sumY                                         /= static_cast<proshade_double> ( detected.at(symIt).size() - 1 );
        sumZ                                         /= static_cast<proshade_double> ( detected.at(symIt).size() - 1 );
        sumH                                         /= static_cast<proshade_double> ( detected.at(symIt).size() - 1 );
        
        //============================================ And add these as well
        hlpP[1]                                       = sumX;
        hlpP[2]                                       = sumY;
        hlpP[3]                                       = sumZ;
        hlpP[5]                                       = sumH;
        
        //============================================ Save the complete symmetry description to the vector, unless already there
        if ( !ProSHADE_internal_symmetry::isSymmetrySame ( ret, hlpP, axErr ) )
        {
            ProSHADE_internal_misc::addToDblPtrVector ( ret, hlpP );
        }
        else
        {
            delete[] hlpP;
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function checks if a very similar symmetry is not already saved.
 
    This is a simple function comparing a single double array of 6 to a vector of these, returning whether the vector already
    contains a very similar entry to the rested one. If the new has better height, replacement will take place.
 
    \param[in] ret This is the variable where the tested array will be saved if passed. It is a vector of double[6] arrays with the following meaning: [0] = fold, [1] = x-axis, [2] = y-axis, [3] = z-axis, [4] = angle, [5] = average peak height.
    \param[in] sym This is a double array of 6 which is to be compared to all the vector entries.
    \param[in] simThres The threshold for axis-wise comparison similarity.
    \param[out] X Boolean value stating whether a similar entry has been found (true = it was, false = it was not).
 */
bool ProSHADE_internal_symmetry::isSymmetrySame ( std::vector< proshade_double* >* ret, proshade_double* sym, proshade_double simThres )
{
    //================================================ Initialise variables
    proshade_double dotProduct                        = 0.0;
    
    //================================================ Check
    for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( ret->size() ); symIt++ )
    {
        //============================================ Minor speed-up => only test for same folds
        if ( ret->at(symIt)[0] == sym[0] )
        {
            //======================================== Is axis the same?
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &ret->at(symIt)[1], &ret->at(symIt)[2],
                                                                                                     &ret->at(symIt)[3], &sym[1], &sym[2], &sym[3] );
            if ( ( ( 1.0 > ( dotProduct - simThres ) ) && ( 1.0 < ( dotProduct + simThres ) ) ) || ( ( -1.0 > ( dotProduct - simThres ) ) && ( -1.0 < ( dotProduct + simThres ) ) ) )
            {
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
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Starting D symmetry detection." );
    
    //================================================If not enough axes, just end here
    if ( CSymList->size() < 2 ) { return ( ret ); }
    
    //================================================ For each unique pair of axes
    for ( proshade_unsign ax1 = 0; ax1 < static_cast<proshade_unsign> ( CSymList->size() ); ax1++ )
    {
        //============================================ Ignore small axes
        if ( CSymList->at(ax1)[5] < settings->minSymPeak ) { continue; }
        
        for ( proshade_unsign ax2 = 1; ax2 < static_cast<proshade_unsign> ( CSymList->size() ); ax2++ )
        {
            //======================================= Use unique pairs only
            if ( ax1 >= ax2 ) { continue; }
            
            //======================================== Ignore small axes
            if ( CSymList->at(ax2)[5] < settings->minSymPeak ) { continue; }
            
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
                else { ProSHADE_internal_symmetry::saveDSymmetry ( &ret, CSymList, ax2, ax1 ); }
            }
        }
    }
    
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "Detected " << ret.size() << " D symmetries.";
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, hlpSS.str() );
    
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
    proshade_double* hlpP                             = new proshade_double [12];
    ProSHADE_internal_misc::checkMemoryAllocation     ( hlpP, __FILE__, __LINE__, __func__ );
    
    //================================================ Set obvious values
    hlpP[0]                                           = static_cast<proshade_double> ( CSymList->at(axisOne)[0] );
    hlpP[4]                                           = static_cast<proshade_double> ( ( 2.0 * M_PI ) / hlpP[0] );
    hlpP[6]                                           = static_cast<proshade_double> ( CSymList->at(axisTwo)[0] );
    hlpP[10]                                          = static_cast<proshade_double> ( ( 2.0 * M_PI ) / hlpP[6] );
    
    //================================================ Set the axis and heights
    hlpP[1]                                           = CSymList->at(axisOne)[1];
    hlpP[2]                                           = CSymList->at(axisOne)[2];
    hlpP[3]                                           = CSymList->at(axisOne)[3];
    hlpP[5]                                           = CSymList->at(axisOne)[5];
    hlpP[7]                                           = CSymList->at(axisTwo)[1];
    hlpP[8]                                           = CSymList->at(axisTwo)[2];
    hlpP[9]                                           = CSymList->at(axisTwo)[3];
    hlpP[11]                                          = CSymList->at(axisTwo)[5];
    
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
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting T symmetry detection." );
    
    //================================================ Are the basic requirements for tetrahedral symmetry met?
    if ( ProSHADE_internal_symmetry::detectTetrahedralSymmetry ( CSymList, settings->axisErrTolerance, settings->minSymPeak ) )
    {
        //============================================ Search for all the symmetry axes
        ProSHADE_internal_symmetry::findTetra4C3s     ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->minSymPeak );
        if ( ret.size() != 4 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        
        ProSHADE_internal_symmetry::findTetra3C2s     ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->minSymPeak );
        if ( ret.size() != 7 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        else
        {
            for ( proshade_unsign csIt = 0; csIt < static_cast<proshade_unsign> ( CSymList->size() ); csIt++ )
            {
                for ( proshade_unsign retIt = 0; retIt < static_cast<proshade_unsign> ( ret.size() ); retIt++ )
                {
                    if ( ( CSymList->at(csIt)[0] == ret.at(retIt)[0] ) &&
                         ( CSymList->at(csIt)[1] == ret.at(retIt)[1] ) &&
                         ( CSymList->at(csIt)[2] == ret.at(retIt)[2] ) &&
                         ( CSymList->at(csIt)[3] == ret.at(retIt)[3] ) &&
                         ( CSymList->at(csIt)[4] == ret.at(retIt)[4] ) &&
                         ( CSymList->at(csIt)[5] == ret.at(retIt)[5] ) )
                    {
                        ProSHADE_internal_misc::addToUnsignVector ( &settings->allDetectedTAxes, csIt );
                    }
                }
            }
        }
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "T symmetry detection complete." );
    
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
        if ( CSymList->at(cSym)[0] == 3 && CSymList->at(cSym)[5] >= minPeakHeight ) { ProSHADE_internal_misc::addToUnsignVector ( &C3List, cSym ); }
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
    \param[in] verobse How loud the announcments should be?
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
 */
void ProSHADE_internal_symmetry::findTetra4C3s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_unsign verbose, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C3PossibilitiesHlp;
    std::vector< std::vector< proshade_unsign > > C3Possibilities;
    bool groupMatched;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage ( verbose, 2, "Starting detection of four C3 axes." );
    
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
        ProSHADE_internal_messages::printProgressMessage ( verbose, 3, "Detection of four C3 axes successfull." );
        
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
                    //================================ Add
                    ProSHADE_internal_misc::addToDblPtrVector ( CSymList, hlpVec.at(axIt) );
                    ProSHADE_internal_misc::addToUnsignVector ( &possibilities->at(gIt), static_cast<proshade_unsign> ( CSymList->size()-1 ) );
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
    proshade_double angStep                           = axErr;
    std::vector< proshade_double* > angVec;
    
    //================================================ Find map points conforming to the axis
    angVec                                            = ProSHADE_internal_symmetry::findMissingAxisPoints ( xVal, yVal, zVal, dataObj, axErr );
    
    //================================================ Sort points by angle
    std::sort                                         ( angVec.begin(), angVec.end(), ProSHADE_internal_symmetry::sortArrVecHlp );
    
    //================================================ Find the best X peaks with correct distances
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( std::floor ( ( 2.0 * M_PI / angStep ) / fold ) ); iter++ )
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
                if ( angVec.at(angIt)[0] < ( ( iter*angStep )     + ( ( 2.0 * M_PI / fold ) * angCmb ) ) ) { continue; }
                if ( angVec.at(angIt)[0] > ( ( (iter+1)*angStep ) + ( ( 2.0 * M_PI / fold ) * angCmb ) ) ) { break; }

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
                ProSHADE_internal_maths::getEulerZXZFromSOFTPosition ( dataObj->getMaxBand(), static_cast<proshade_signed> ( xIt ),
                                                                       static_cast<proshade_signed> ( yIt ), static_cast<proshade_signed> ( zIt ),
                                                                       &euA, &euB, &euG );
                ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( euA, euB, euG, rotMat );
                ProSHADE_internal_maths::getAxisAngleFromRotationMatrix ( rotMat, &xPk, &yPk, &zPk, &anglPk );
                
                //==================================== Set largest axis element to positive
                if ( ( ( std::max ( std::abs ( xPk ), std::max( std::abs ( yPk ), std::abs ( zPk ) ) ) == std::abs ( xPk ) ) && ( xPk < 0.0 ) ) ||
                     ( ( std::max ( std::abs ( xPk ), std::max( std::abs ( yPk ), std::abs ( zPk ) ) ) == std::abs ( yPk ) ) && ( yPk < 0.0 ) ) ||
                     ( ( std::max ( std::abs ( xPk ), std::max( std::abs ( yPk ), std::abs ( zPk ) ) ) == std::abs ( zPk ) ) && ( zPk < 0.0 ) ) )
                {
                    xPk                              *= -1.0;
                    yPk                              *= -1.0;
                    zPk                              *= -1.0;
                    anglPk                           *= -1.0;
                }
                
                //==================================== Does the peak match the required axis?
                if ( ( ( ( xPk - axErr ) < xVal ) && ( ( xPk + axErr ) > xVal ) ) &&
                     ( ( ( yPk - axErr ) < yVal ) && ( ( yPk + axErr ) > yVal ) ) &&
                     ( ( ( zPk - axErr ) < zVal ) && ( ( zPk + axErr ) > zVal ) ) )
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
        if ( axVec->at(symIt)[0] == hlpSym[0] )
        {
            if ( ( ( ( axVec->at(symIt)[1] + axErr ) > hlpSym[1] ) && ( ( axVec->at(symIt)[1] - axErr ) < hlpSym[1] ) ) &&
                 ( ( ( axVec->at(symIt)[2] + axErr ) > hlpSym[2] ) && ( ( axVec->at(symIt)[2] - axErr ) < hlpSym[2] ) ) &&
                 ( ( ( axVec->at(symIt)[3] + axErr ) > hlpSym[3] ) && ( ( axVec->at(symIt)[3] - axErr ) < hlpSym[3] ) ) )
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
    proshade_double* symHlp                           = new proshade_double[6];
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
            if ( ( ( std::max ( std::abs ( solVec.at(0) ), std::max( std::abs ( solVec.at(1) ), std::abs ( solVec.at(2) ) ) ) == std::abs ( solVec.at(0) ) ) && ( solVec.at(0) < 0.0 ) ) ||
                 ( ( std::max ( std::abs ( solVec.at(0) ), std::max( std::abs ( solVec.at(1) ), std::abs ( solVec.at(2) ) ) ) == std::abs ( solVec.at(1) ) ) && ( solVec.at(1) < 0.0 ) ) ||
                 ( ( std::max ( std::abs ( solVec.at(0) ), std::max( std::abs ( solVec.at(1) ), std::abs ( solVec.at(2) ) ) ) == std::abs ( solVec.at(2) ) ) && ( solVec.at(2) < 0.0 ) ) )
            {
                solVec.at(0)                         *= -1.0;
                solVec.at(1)                         *= -1.0;
                solVec.at(2)                         *= -1.0;
            }
            
            //======================================== Does the solution fit the whole group?
            symHlp[1] = solVec.at(0); symHlp[2] = solVec.at(1); symHlp[3] = solVec.at(2);
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
            if ( ( ( std::max ( std::abs ( solVec.at(0) ), std::max( std::abs ( solVec.at(1) ), std::abs ( solVec.at(2) ) ) ) == std::abs ( solVec.at(0) ) ) && ( solVec.at(0) < 0.0 ) ) ||
                 ( ( std::max ( std::abs ( solVec.at(0) ), std::max( std::abs ( solVec.at(1) ), std::abs ( solVec.at(2) ) ) ) == std::abs ( solVec.at(1) ) ) && ( solVec.at(1) < 0.0 ) ) ||
                 ( ( std::max ( std::abs ( solVec.at(0) ), std::max( std::abs ( solVec.at(1) ), std::abs ( solVec.at(2) ) ) ) == std::abs ( solVec.at(2) ) ) && ( solVec.at(2) < 0.0 ) ) )
            {
                solVec.at(0)                         *= -1.0;
                solVec.at(1)                         *= -1.0;
                solVec.at(2)                         *= -1.0;
            }
            
            //======================================== Does the solution fit the whole group?
            symHlp[1] = solVec.at(0); symHlp[2] = solVec.at(1); symHlp[3] = solVec.at(2);
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
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
 */
void ProSHADE_internal_symmetry::findTetra3C2s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_unsign verbose, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C3s, prospectiveC2s, C2PossibilitiesHlp;
    std::vector< std::vector< proshade_unsign > > C2Possibilities;
    proshade_double dotProd;
    bool groupMatched;
    for ( proshade_unsign iter = 0; iter < 4; iter++ ) { ProSHADE_internal_misc::addToUnsignVector ( &C3s, iter ); }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of three C2 axes." );
    
    //================================================ For each C3
    for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( ret->size() ); rIt++ )
    {
        //============================================ For each C2, check it has angle ( acos(0.5) ) to the tested C3
        for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
        {
            //======================================== Search only using C2s
            if ( CSymList->at(cIt)[0] != 2.0 || CSymList->at(cIt)[5] < minPeakHeight ) { continue; }
            
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
            ProSHADE_internal_messages::printProgressMessage ( verbose, 3, "Detection of three C2 axes successfull." );
            
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
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting O symmetry detection." );
    
    //================================================ Are the basic requirements for tetrahedral symmetry met?
    if ( ProSHADE_internal_symmetry::detectOctahedralSymmetry ( CSymList, settings->axisErrTolerance, settings->minSymPeak ) )
    {
        //============================================ Search for all the symmetry axes
        ProSHADE_internal_symmetry::findOcta3C4s ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->minSymPeak );
        if ( ret.size() != 3 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        
        ProSHADE_internal_symmetry::findOcta4C3s ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->minSymPeak );
        if ( ret.size() != 7 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        
        ProSHADE_internal_symmetry::findOcta6C2s ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->minSymPeak );
        if ( ret.size() != 13 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        else
        {
            for ( proshade_unsign csIt = 0; csIt < static_cast<proshade_unsign> ( CSymList->size() ); csIt++ )
            {
                for ( proshade_unsign retIt = 0; retIt < static_cast<proshade_unsign> ( ret.size() ); retIt++ )
                {
                    if ( ( CSymList->at(csIt)[0] == ret.at(retIt)[0] ) &&
                         ( CSymList->at(csIt)[1] == ret.at(retIt)[1] ) &&
                         ( CSymList->at(csIt)[2] == ret.at(retIt)[2] ) &&
                         ( CSymList->at(csIt)[3] == ret.at(retIt)[3] ) &&
                         ( CSymList->at(csIt)[4] == ret.at(retIt)[4] ) &&
                         ( CSymList->at(csIt)[5] == ret.at(retIt)[5] ) )
                    {
                        ProSHADE_internal_misc::addToUnsignVector ( &settings->allDetectedOAxes, csIt );
                    }
                }
            }
        }
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "O symmetry detection complete." );
    
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
        if ( CSymList->at(cSym)[0] == 4 && CSymList->at(cSym)[5] >= minPeakHeight ) { ProSHADE_internal_misc::addToUnsignVector ( &C4List, cSym ); }
    }
    
    //================================================ For each unique pair of C3s
    for ( proshade_unsign c4 = 0; c4 < static_cast<proshade_unsign> ( C4List.size() ); c4++ )
    {
        for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ )
        {
            //======================================== Compare only C3s to the C3List
            if ( CSymList->at(cSym)[0] != 3 ) { continue; }
            
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
    \param[in] verobse How loud the announcments should be?
 */
void ProSHADE_internal_symmetry::findOcta3C4s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_unsign verbose, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C4PossibilitiesHlp;
    std::vector< std::vector< proshade_unsign > > C4Possibilities;
    bool groupMatched;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of three C4 axes." );
    
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
        ProSHADE_internal_messages::printProgressMessage ( verbose, 3, "Detection of three C4 axes successfull." );
        
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
    \param[in] verobse How loud the announcments should be?
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
    \param[in] TetraSymList A vector containing the already detected tetrahedral symmetries - this is to avoid the same search for four C3 symmetry axes.
 */
void ProSHADE_internal_symmetry::findOcta4C3s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_unsign verbose, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C4s, prospectiveC3s, C3PossibilitiesHlp;
    std::vector< std::vector< proshade_unsign > > C3Possibilities;
    proshade_double dotProd;
    bool groupMatched;
    for ( proshade_unsign iter = 0; iter < 3; iter++ ) { ProSHADE_internal_misc::addToUnsignVector ( &C4s, iter ); }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of four C3 axes." );
    
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
            ProSHADE_internal_messages::printProgressMessage ( verbose, 3, "Detection of four C3 axes successfull." );

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
    \param[in] verobse How loud the announcments should be?
    \param[in] minPeakHeight The minimum average peak height for axis to be considered.
 */
void ProSHADE_internal_symmetry::findOcta6C2s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_unsign verbose, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > prospectiveC2s, retGrp;
    proshade_double dotProd;
    proshade_unsign noPerpendicular, noSqrtTwo;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of six C2 axes." );
    
    //================================================ For each C2
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
    {
        //============================================ Use only C2s
        if ( CSymList->at(cIt)[0] != 2.0 || CSymList->at(cIt)[5] < minPeakHeight ) { continue; }
        
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
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "Detection of six C2 axes successfull." );
    
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
        ProSHADE_internal_symmetry::addAxisUnlessSame ( CSymList->at(possibilities->at(prIt))[0],
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
            ProSHADE_internal_misc::addToUnsignVector ( possibilities, CSymList->size() );
            ProSHADE_internal_misc::addToDblPtrVector ( CSymList, prosp.at(iter) );
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
    \param[in] prosp The vector to which the axis is to be saved.
    \param[in] axErr The error tolerance on angle matching.
 */
void ProSHADE_internal_symmetry::addAxisUnlessSame ( proshade_unsign fold, proshade_double axX, proshade_double axY, proshade_double axZ, proshade_double axHeight, std::vector< proshade_double* >* prosp, proshade_double axErr )
{
    //================================================ Initialise variables
    proshade_double* symHlp                           = new proshade_double[6];
    ProSHADE_internal_misc::checkMemoryAllocation     ( symHlp, __FILE__, __LINE__, __func__ );
    
    //================================================ Fill in the prospective axis
    symHlp[0]                                         = static_cast<proshade_double> ( fold );
    symHlp[1]                                         = axX;
    symHlp[2]                                         = axY;
    symHlp[3]                                         = axZ;
    symHlp[4]                                         = 2.0 * M_PI / symHlp[0];
    symHlp[5]                                         = axHeight;
    
    //================================================ If it is not the same as already saved axes
    if ( !ProSHADE_internal_symmetry::isSymmetrySame ( prosp, symHlp, axErr ) )
    {
        ProSHADE_internal_misc::addToDblPtrVector     ( prosp, symHlp );
    }
    else
    {
        delete[] symHlp;
    }
    
    //================================================ Done
    return ;
    
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
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting I symmetry detection." );
    
    //================================================ Are the basic requirements for tetrahedral symmetry met?
    if ( ProSHADE_internal_symmetry::detectIcosahedralSymmetry ( CSymList, settings->axisErrTolerance, settings->minSymPeak ) )
    {
        //============================================ Search for all the symmetry axes
        ProSHADE_internal_symmetry::findIcos6C5s      ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->minSymPeak );
        if ( ret.size() != 6 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }

        ProSHADE_internal_symmetry::findIcos10C3s     ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->minSymPeak );
        if ( ret.size() != 16 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }

        ProSHADE_internal_symmetry::findIcos15C2s     ( CSymList, &ret, settings->axisErrTolerance, this, settings->verbose, settings->minSymPeak );
        if ( ret.size() != 31 ) { ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Failed to detect some of the polyhedral symmetries, while detecting the correct dihedral angles.", "WS00031" ); return ( ret ); }
        else
        {
            for ( proshade_unsign csIt = 0; csIt < static_cast<proshade_unsign> ( CSymList->size() ); csIt++ )
            {
                for ( proshade_unsign retIt = 0; retIt < static_cast<proshade_unsign> ( ret.size() ); retIt++ )
                {
                    if ( ( CSymList->at(csIt)[0] == ret.at(retIt)[0] ) &&
                         ( CSymList->at(csIt)[1] == ret.at(retIt)[1] ) &&
                         ( CSymList->at(csIt)[2] == ret.at(retIt)[2] ) &&
                         ( CSymList->at(csIt)[3] == ret.at(retIt)[3] ) &&
                         ( CSymList->at(csIt)[4] == ret.at(retIt)[4] ) &&
                         ( CSymList->at(csIt)[5] == ret.at(retIt)[5] ) )
                    {
                        ProSHADE_internal_misc::addToUnsignVector ( &settings->allDetectedIAxes, csIt );
                    }
                }
            }
        }
    }

    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "I symmetry detection complete." );

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
        if ( CSymList->at(cSym)[0] == 5 && CSymList->at(cSym)[5] >= minPeakHeight ) { ProSHADE_internal_misc::addToUnsignVector ( &C5List, cSym ); }
    }
    
    //================================================ For each unique pair of C5 and C3
    for ( proshade_unsign c5 = 0; c5 < static_cast<proshade_unsign> ( C5List.size() ); c5++ )
    {
        for ( proshade_unsign cSym = 0; cSym < static_cast<proshade_unsign> ( CSymList->size() ); cSym++ )
        {
            //======================================== Compare only C3s to the C5List
            if ( CSymList->at(cSym)[0] != 3 ) { continue; }
            
            //========================================  Check the angle between the C5 and C3 axes
            dotProduct                                = ProSHADE_internal_maths::computeDotProduct ( &CSymList->at(C5List.at(c5))[1],
                                                                                                     &CSymList->at(C5List.at(c5))[2],
                                                                                                     &CSymList->at(C5List.at(c5))[3],
                                                                                                     &CSymList->at(cSym)[1],
                                                                                                     &CSymList->at(cSym)[2],
                                                                                                     &CSymList->at(cSym)[3] );
            
            //======================================== Is the angle approximately the dihedral angle
            if ( ( ( sqrt ( 5.0 ) / 3.0 ) > ( dotProduct - axErr ) ) && ( ( sqrt ( 5.0 ) / 3.0 ) < ( dotProduct + axErr ) ) )
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
    \param[in] verobse How loud the announcments should be?
 */
void ProSHADE_internal_symmetry::findIcos6C5s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_unsign verbose, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > C5PossibilitiesHlp;
    std::vector< std::vector< proshade_unsign > > C5Possibilities;
    bool groupMatched;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of six C5 axes." );

    //================================================ For all symmetries in the C symmetries list
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
    {
        //============================================ Search only using C5s and check peak height
        if ( CSymList->at(cIt)[0] != 5.0 || CSymList->at(cIt)[5] < minPeakHeight ) { continue; }

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

    //======================================== ========Any group has 6 entries? If more such groups, take the one with highest average height.
    proshade_double maxHeight = 0.0; proshade_unsign maxGrp = 0;
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( C5Possibilities.size() ); iter++ ) { if ( C5Possibilities.at(iter).size() == 6 ) { if ( ( ( CSymList->at(C5Possibilities.at(iter).at(0))[5] + CSymList->at(C5Possibilities.at(iter).at(1))[5] + CSymList->at(C5Possibilities.at(iter).at(2))[5] + CSymList->at(C5Possibilities.at(iter).at(3))[5] + CSymList->at(C5Possibilities.at(iter).at(4))[5] + CSymList->at(C5Possibilities.at(iter).at(5))[5] ) / 6.0 ) > maxHeight ) { maxHeight = ( ( CSymList->at(C5Possibilities.at(iter).at(0))[5] + CSymList->at(C5Possibilities.at(iter).at(1))[5] + CSymList->at(C5Possibilities.at(iter).at(2))[5] + CSymList->at(C5Possibilities.at(iter).at(3))[5] + CSymList->at(C5Possibilities.at(iter).at(4))[5] + CSymList->at(C5Possibilities.at(iter).at(5))[5] ) / 6.0 ); maxGrp = iter; } } }

    if ( C5Possibilities.at(maxGrp).size() == 6 )
    {
        //============================================ Success! Save and exit
        for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( C5Possibilities.at(maxGrp).size() ); it++ ) { ProSHADE_internal_misc::addToDblPtrVector ( ret, CSymList->at(C5Possibilities.at(maxGrp).at(it)) ); }

        //============================================ Report progress
        ProSHADE_internal_messages::printProgressMessage ( verbose, 3, "Detection of six C5 axes successfull." );

        //============================================ Done
        return ;
    }
    
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
    \param[in] verobse How loud the announcments should be?
 */
void ProSHADE_internal_symmetry::findIcos10C3s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_unsign verbose, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > prospectiveC3s, retGrp;
    proshade_double dotProd;
    proshade_unsign noClose, noAway;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of ten C3 axes." );

    //================================================ For each C3
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
    {
        //============================================ Use only C3s with hight enough average
        if ( CSymList->at(cIt)[0] != 3.0 || CSymList->at(cIt)[0] < minPeakHeight ) { continue; }

        //============================================ Check the C3 has acos ( sqrt ( 5.0 ) / 3.0 ) to 3 C5s and acos ( 1.0 - ( sqrt ( 5.0 ) / 3.0 ) ) to the other three C5s
        noClose = 0; noAway = 0;
        for ( proshade_unsign rIt = 0; rIt < 6; rIt++ )
        {
            dotProd                                   = ProSHADE_internal_maths::computeDotProduct ( &ret->at(rIt)[1],
                                                                                                     &ret->at(rIt)[2],
                                                                                                     &ret->at(rIt)[3],
                                                                                                     &CSymList->at(cIt)[1],
                                                                                                     &CSymList->at(cIt)[2],
                                                                                                     &CSymList->at(cIt)[3] );

            if ( ( std::abs ( dotProd ) > (       ( sqrt ( 5.0 ) / 3.0 ) - axErr ) ) && ( std::abs ( dotProd ) < (       ( sqrt ( 5.0 ) / 3.0 ) + axErr ) ) ) { noClose += 1; continue; }
            if ( ( std::abs ( dotProd ) > ( 1.0 - ( sqrt ( 5.0 ) / 3.0 ) - axErr ) ) && ( std::abs ( dotProd ) < ( 1.0 - ( sqrt ( 5.0 ) / 3.0 ) + axErr ) ) ) { noAway  += 1; continue; }
        }

        //============================================ If correct angles distribution is found, save the axis
        if ( ( noClose == 3 ) && ( noAway == 3 ) )
        {
            ProSHADE_internal_misc::addToUnsignVector ( &prospectiveC3s, cIt );
        }
    }

    //================================================ Search for missing axes
    for ( proshade_unsign iter = 0; iter < 6; iter++ ) { ProSHADE_internal_misc::addToUnsignVector ( &retGrp, iter ); }
    if ( !ProSHADE_internal_symmetry::findMissingAxesDual ( &prospectiveC3s, CSymList, ret, &retGrp, 10, axErr, 3, sqrt ( 5.0 ) / 3.0, 3, 1.0 - ( sqrt ( 5.0 ) / 3.0 ), 3, dataObj ) )
    {
        return ;
    }

    //================================================ Found correct number of axes! Now save the
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( prospectiveC3s.size() ); iter++ )
    {
        ProSHADE_internal_misc::addToDblPtrVector     ( ret, CSymList->at(prospectiveC3s.at(iter)) );
    }

    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "Detection of ten C3 axes successfull." );
    
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
    \param[in] verobse How loud the announcments should be?
 */
void ProSHADE_internal_symmetry::findIcos15C2s ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_unsign verbose, proshade_double minPeakHeight )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > prospectiveC2s, retGrp;
    proshade_double dotProd;
    proshade_unsign noClose, noMidway, noAway;
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 2, "Starting detection of fifteen C2 axes." );

    //================================================ For each C2
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSymList->size() ); cIt++ )
    {
        //============================================ Use only C2s
        if ( CSymList->at(cIt)[0] != 2.0 || CSymList->at(cIt)[0] < minPeakHeight ) { continue; }

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
    ProSHADE_internal_messages::printProgressMessage  ( verbose, 3, "Detection of fifteen C2 axes successfull." );
    
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
        ProSHADE_internal_symmetry::addAxisUnlessSame ( CSymList->at(possibilities->at(prIt))[0],
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
            //======================================== Add
            ProSHADE_internal_misc::addToDblPtrVector ( CSymList, prosp.at(axIt) );
            ProSHADE_internal_misc::addToUnsignVector ( possibilities, static_cast<proshade_unsign> ( CSymList->size()-1 ) );
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
