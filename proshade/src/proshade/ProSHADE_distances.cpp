/*! \file ProSHADE_distances.cpp
    \brief This is the source file containing functions required for computation of shape distances.

    This source file contains the functions required to compute shape distances (specifically energy levels distances, trace sigma distances and the
    full rotation function distances) between two shapes. The user should not need to access these functions directly, as there are automated functions
    available in the higher levels of the ProSHADE organisation, which will call these in correct order and parse the results properly.

    Copyright by Michal Tykac and individual contributors. All rights reserved.
    
    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.
    
    This software is provided by the copyright holders and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In     no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data     or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility     of such damage.
    
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.6.2
    \date      DEC 2021
 */

//==================================================== ProSHADE
#include "ProSHADE_distances.hpp"

/*! \brief This function allocates the required memory for the RRP matrices.
 
    This function belongs to the ProSHADE_data class and its role is to allocate the require memory
    for the RRP matrices, given the already determined bandwidths and shell count.
 */
void ProSHADE_internal_data::ProSHADE_data::allocateRRPMemory ( )
{
    //================================================ Allocate the required memory
    this->rrpMatrices                                 = new proshade_double** [this->maxShellBand];
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->rrpMatrices, __FILE__, __LINE__, __func__ );
    
    for ( proshade_unsign bwIt = 0; bwIt < this->maxShellBand; bwIt++ )
    {
        //============================================ For rach sphere
        this->rrpMatrices[bwIt]                       = new proshade_double*  [this->noSpheres];
        ProSHADE_internal_misc::checkMemoryAllocation ( this->rrpMatrices[bwIt], __FILE__, __LINE__, __func__ );
        
        for ( proshade_unsign shIt = 0; shIt < this->noSpheres; shIt++ )
        {
            this->rrpMatrices[bwIt][shIt]             = new double [this->noSpheres];
            ProSHADE_internal_misc::checkMemoryAllocation ( this->rrpMatrices[bwIt][shIt], __FILE__, __LINE__, __func__ );
        }
    }
}

/*! \brief This function pre-computes the RRP matrices for a data object.
 
    This function belongs to the ProSHADE_data class and its role is to set the objects internal
    variables properly and provide all the required calculations, so that the object will in the
    end have all the RRP matrices computed and be ready for the energy levels calculation.
 
    \param[in] settings A pointer to settings class containing all the information required for the task.
 */
void ProSHADE_internal_data::ProSHADE_data::computeRRPMatrices ( ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Computing RRP matrices for structure " + this->fileName, settings->messageShift );
    
    //================================================ Allocate the memory
    this->allocateRRPMemory                           ( );
    
    //================================================ Start computation: For each band (l)
    proshade_double descValR                          = 0.0;
    proshade_unsign arrPos1, arrPos2;
    for ( proshade_unsign band = 0; band < this->maxShellBand; band++ )
    {
        //============================================ For each unique shell couple
        for ( proshade_unsign shell1 = 0; shell1 < this->noSpheres; shell1++ )
        {
            //======================================== Does the band exist for this shell1?
            if ( !ProSHADE_internal_distances::isBandWithinShell ( band, shell1, this->spheres ) )
            {
                for ( proshade_unsign shell2 = 0; shell2 < this->noSpheres; shell2++ )
                {
                    this->rrpMatrices[band][shell1][shell2] = 0.0;
                    this->rrpMatrices[band][shell2][shell1] = 0.0;
                }
                continue;
            }
            
            for ( proshade_unsign shell2 = 0; shell2 < this->noSpheres; shell2++ )
            {
                //==================================== Compute each values only once
                if ( shell1 > shell2 ) { continue; }
                
                //==================================== Check if band exists for this shell2?
                if ( !ProSHADE_internal_distances::isBandWithinShell ( band, shell2, this->spheres ) )
                {
                    this->rrpMatrices[band][shell1][shell2] = 0.0;
                    this->rrpMatrices[band][shell2][shell1] = 0.0;
                    continue;
                }

                //==================================== Initialise
                descValR                              = 0.0;

                //==================================== Sum over order (m)
                for ( proshade_unsign order = 0; order < static_cast< proshade_unsign >  ( ( 2 * band ) + 1 ); order++ )
                {
                    arrPos1                           = static_cast< proshade_unsign > ( seanindex ( static_cast< int > ( order ) - static_cast<int > ( band ),
                                                                                                   static_cast< int > ( band ), static_cast< int > ( this->spheres[shell1]->getLocalBandwidth() ) ) );
                    arrPos2                           = static_cast< proshade_unsign > ( seanindex ( static_cast< int > ( order ) - static_cast< int > ( band ),
                                                                                                   static_cast< int > ( band ), static_cast< int > ( this->spheres[shell2]->getLocalBandwidth() ) ) );
                    descValR                         += ProSHADE_internal_maths::complexMultiplicationConjugRealOnly ( &this->sphericalHarmonics[shell1][arrPos1][0],
                                                                                                                       &this->sphericalHarmonics[shell1][arrPos1][1],
                                                                                                                       &this->sphericalHarmonics[shell2][arrPos2][0],
                                                                                                                       &this->sphericalHarmonics[shell2][arrPos2][1]  );
                }

                //==================================== Save the matrices
                this->rrpMatrices[band][shell1][shell2] = descValR;
                this->rrpMatrices[band][shell2][shell1] = descValR;
            }
        }
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, "RRP matrices successfully computed.", settings->messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function checks if a band is available for a given shell.
 
    This function simply checks if a particular sphere bandwidth limit is higher than a requested
    band, returning true if it is and false if not.
 
    \param[in] bandInQuestion The value of the band existence of which is to be checked.
    \param[in] shellInQuestion The index of the shell for which the band existence is checked.
    \param[in] spheres The ProSHADE structure holding all the shells and their information.
 */
bool ProSHADE_internal_distances::isBandWithinShell ( proshade_unsign bandInQuestion, proshade_unsign shellInQuestion, ProSHADE_internal_spheres::ProSHADE_sphere** spheres )
{
    if ( bandInQuestion < spheres[shellInQuestion]->getLocalBandwidth() )
    {
        return                                        ( true );
    }
    else
    {
        return                                        ( false );
    }
}

/*! \brief This function computes the energy levels descriptor value between two objects.
 
    This function is where the enery levels descriptor computation is controlled and done from. It starts by
    making sure that both input data objects have the RRP matrices computed and then it proceeds to compute the
    Pearson's coefficients for each band, finally averaging the band values and returning the descriptos.
 
    \param[in] obj1 The first ProSHADE_data object against which comparison is done.
    \param[in] obj2 The second ProSHADE_data object which is compared to the first.
    \param[in] settings A pointer to settings class containing all the information required for the task.
 */
proshade_double ProSHADE_internal_distances::computeEnergyLevelsDescriptor ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2, ProSHADE_settings* settings )
{
    //================================================ Report starting the task
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting energy levels distance computation.", settings->messageShift );
    
    //================================================ Initialise local variables
    proshade_double ret                               = 0.0;
    std::vector<proshade_double> bandDists;
    
    //================================================ Sanity check
    if ( !settings->computeEnergyLevelsDesc )
    {
        throw ProSHADE_exception ( "Attempted computing energy levels descriptors when it was not required.", "ED00017", __FILE__, __LINE__, __func__, "Attempted to pre-compute the RRP matrices, when the user\n                    : has specifically stated that these should not be computed.\n                    : Unless you manipulated the code, this error should never\n                    : occur; if you see this, I made a large blunder. Please let\n                    : me know!" );
    }
    
    //================================================ Get the RRP matrices for both objects
    obj1->computeRRPMatrices                          ( settings );
    obj2->computeRRPMatrices                          ( settings );
    
    //================================================ Find the minimium comparable shells and bands
    proshade_unsign minCommonShells                   = std::min ( obj1->getMaxSpheres(), obj2->getMaxSpheres() );
    proshade_unsign minCommonBands                    = std::min ( obj1->getMaxBand(),    obj2->getMaxBand()    );

    //================================================ Get the Pearson's coefficients for each common band
    computeRRPPearsonCoefficients                     ( obj1, obj2, settings, minCommonBands, minCommonShells, &bandDists );
    
    //================================================ Get distance (by averaging Patterson's coefficients)
    ret                                               = static_cast<proshade_double> ( std::accumulate ( bandDists.begin(), bandDists.end(), 0.0 ) ) /
                                                                                       static_cast<proshade_double> ( bandDists.size() );
    
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Energy levels distance computation complete.", settings->messageShift );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function gets the Pearson's coefficients or all bands between two objects.
 
    This function takes two data objects with their RRP matrices computed and proceeds to compute the
    Pearson's correlation coefficient for each band, saving it into the supplied vector.
 
    \param[in] obj1 The first ProSHADE_data object against which comparison is done.
    \param[in] obj2 The second ProSHADE_data object which is compared to the first.
    \param[in] settings A pointer to settings class containing all the information required for the task.
    \param[in] minCommonBands The number of common bands between the two objects.
    \param[in] minCommonShells The index of highest common shell in both objects.
    \param[in] bandDists Empty vector of proshade_doubles to which the Pearson's Coefficients will be saved for each band.
 */
void ProSHADE_internal_distances::computeRRPPearsonCoefficients ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2, ProSHADE_settings* settings, proshade_unsign minCommonBands, proshade_unsign minCommonShells, std::vector<proshade_double>* bandDists )
{
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Correlating RRP matrices.", settings->messageShift );
    
    //================================================ Initialise local variables
    proshade_double *str1Vals                         = new proshade_double[minCommonShells * minCommonShells];
    proshade_double *str2Vals                         = new proshade_double[minCommonShells * minCommonShells];
    ProSHADE_internal_misc::checkMemoryAllocation     ( str1Vals, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( str2Vals, __FILE__, __LINE__, __func__ );
    proshade_unsign arrIter                           = 0;

    //================================================ Start computation: For each band (l)
    for ( proshade_unsign band = 0; band < minCommonBands; band++ )
    {
        //============================================ Reset local counter
        arrIter                                       = 0;

        //============================================ For each shell pair
        for ( proshade_unsign shell1 = 0; shell1 < minCommonShells; shell1++ )
        {
            //======================================== Check if band exists (progressive only)
            if ( settings->progressiveSphereMapping ) { if ( !obj1->shellBandExists( shell1, band ) || !obj2->shellBandExists( shell1, band ) ) { continue; } }
            
            for ( proshade_unsign shell2 = 0; shell2 < minCommonShells; shell2++ )
            {
                //============================ Check the other shell as well
                if ( !obj1->shellBandExists( shell2, band ) || !obj2->shellBandExists( shell2, band ) ) { continue; }
                
                //==================================== Set values between which the Person's correlation coefficient should be computed
                str1Vals[arrIter]                     = obj1->getRRPValue ( band, shell1, shell2 ) *
                                                        pow ( static_cast<proshade_double> ( shell1 ), settings->enLevMatrixPowerWeight ) *
                                                        pow ( static_cast<proshade_double> ( shell2 ), settings->enLevMatrixPowerWeight );
                str2Vals[arrIter]                     = obj2->getRRPValue ( band, shell1, shell2 ) *
                                                        pow ( static_cast<proshade_double> ( shell1 ), settings->enLevMatrixPowerWeight ) *
                                                        pow ( static_cast<proshade_double> ( shell2 ), settings->enLevMatrixPowerWeight );
        
                arrIter                              += 1;
            }
        }

        //============================================ Get Pearson's Correlation Coefficient
        ProSHADE_internal_misc::addToDoubleVector     ( bandDists, ProSHADE_internal_maths::pearsonCorrCoeff ( str1Vals, str2Vals, arrIter ) );
    }

    //================================================ Clean up
    delete[] str1Vals;
    delete[] str2Vals;
    
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, "RRP matrices correlation computed.", settings->messageShift );
    
    //================================================ Done
    return ;
}

/*! \brief This function allocates the required memory for the E matrices.
 
    This function belongs to the ProSHADE_data class and its role is to allocate the require memory
    for the E matrices required by both, the Trace Sigma and Full Rotational descriptors, as well as
    symmetry and rotation tasks.
 
    \param[in] band The minimal band of the comparison for which E matrices are computed.
 */
void ProSHADE_internal_data::ProSHADE_data::allocateEMatrices ( proshade_unsign  band )
{
    //================================================ Save the maximum band to the object
    this->maxCompBand                                 = band;
    
    //================================================ Allocate the required memory
    this->eMatrices                                   = new proshade_complex** [this->maxCompBand];
    ProSHADE_internal_misc::checkMemoryAllocation ( this->eMatrices, __FILE__, __LINE__, __func__ );
    
    for ( proshade_unsign bandIter = 0; bandIter < this->maxCompBand; bandIter++ )
    {
        //============================================ Allocate the data structure
        this->eMatrices[bandIter]                     = new proshade_complex*  [static_cast<proshade_unsign> ( ( bandIter * 2 ) + 1 )];
        ProSHADE_internal_misc::checkMemoryAllocation ( this->eMatrices[bandIter], __FILE__, __LINE__, __func__ );
        
        for ( proshade_unsign band2Iter = 0; band2Iter < static_cast<proshade_unsign> ( ( bandIter * 2 ) + 1 ); band2Iter++ )
        {
            this->eMatrices[bandIter][band2Iter]      = new proshade_complex [static_cast<proshade_unsign> ( ( bandIter * 2 ) + 1 )];
            ProSHADE_internal_misc::checkMemoryAllocation ( this->eMatrices[bandIter][band2Iter], __FILE__, __LINE__, __func__ );
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This helper function is responsible for allocating the workspace memory required for trace sigma descriptor computation.
 
    \param[in] minSpheres The minima of the number of spheres available in the compared objects.
    \param[in] intOrder The integration order for the computation.
    \param[in] obj1Vals Array to hold the shell values for the first object integgration.
    \param[in] obj2Vals Array to hold the shell values for the second object integgration.
    \param[in] GLabscissas An array to hold the pre-computed anscissas for the Gauss-Legendre integration.
    \param[in] glWeights An array to hold the pre-computed weights for the Gauss-Legendre integration.
    \param[in] radiiVals A complex array to hold the results of combining spherical harmonics coefficients of the two objects for each shell.
 */
void ProSHADE_internal_distances::allocateTrSigmaWorkspace ( proshade_unsign minSpheres, proshade_unsign intOrder, proshade_double*& obj1Vals, proshade_double*& obj2Vals, proshade_double*& GLabscissas, proshade_double*& GLweights, proshade_complex*& radiiVals )
{
    //================================================ Allocate the memory
    obj1Vals                                          = new proshade_double [minSpheres];
    obj2Vals                                          = new proshade_double [minSpheres];
    radiiVals                                         = new proshade_complex[minSpheres];
    GLabscissas                                       = new proshade_double [intOrder];
    GLweights                                         = new proshade_double [intOrder];
    
    //================================================ Check the memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( obj1Vals, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( obj2Vals, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( radiiVals, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( GLabscissas, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( GLweights, __FILE__, __LINE__, __func__ );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the magnitude of a particular spherical harmonics position for a given object, weighting it by the radius^2 (for integration).
 
    \param[in] obj The ProSHADE_data object for which the computation is to be done.
    \param[in] band The bandwidth of the SH value for which this should be done.
    \param[in] order The order of the SH value for which this should be done.
    \param[in] radius The shell of the SH value for which this should be done.
    \param[in] result The location where the result is to be saved.
 */
void ProSHADE_internal_distances::computeSphericalHarmonicsMagnitude ( ProSHADE_internal_data::ProSHADE_data* obj, proshade_unsign band, proshade_unsign order, proshade_unsign radius, proshade_double* result )
{
    //================================================ Find the magnitude
   *result                                            = ProSHADE_internal_maths::complexMultiplicationConjugRealOnly ( obj->getRealSphHarmValue ( band, order, radius ),
                                                                                                                       obj->getImagSphHarmValue ( band, order, radius ),
                                                                                                                       obj->getRealSphHarmValue ( band, order, radius ),
                                                                                                                       obj->getImagSphHarmValue ( band, order, radius ) );
    
    //================================================ Weight by radius^2 for the integration that will follow
   *result                                           *= pow ( static_cast<proshade_double> ( obj->getAnySphereRadius( radius ) ), 2.0 );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the E matrix un-weighted values for a given band and order and saves these into the obj2 parameter.
 
    \param[in] obj1 The ProSHADE_data object for which the comparison is done in regards to.
    \param[in] obj2 The ProSHADE_data object for which the comparison is done in regards from - the E matrices will be saved into this object.
    \param[in] band The bandwidth of the SH value for which this should be done.
    \param[in] order The order of the SH value for which this should be done.
    \param[in] radiiVals Already allocated array of proshade_complex to which integrated values will be saved. It must have size equal to minimum of spheres in the two compared objects.
    \param[in] integOrder The Gauss-Legendre integration order to be used.
    \param[in] abscissas The pre-computed abscissas for the Gauss-Legendre integration.
    \param[in] weights The pre-computed weights for the Gauss-Legendre integration.
    \param[in] integRange The range in angstroms between the smalleds and largest shell which are integrated over (might not be 0 to max for progressive shell sampling).
    \param[in] sphereDist The distance between any two spheres.
 */
void ProSHADE_internal_distances::computeEMatricesForLM ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2, proshade_unsign bandIter, proshade_unsign orderIter, proshade_complex* radiiVals, proshade_unsign integOrder, proshade_double* abscissas, proshade_double* weights, proshade_double integRange, proshade_double sphereDist )
{
    //================================================ Initialise local variables
    proshade_unsign objCombValsIter                   = 0;
    proshade_double hlpReal, hlpImag;
    proshade_complex arrVal;
    
    //================================================ For each combination of m and m' for E matrices
    for ( proshade_unsign order2Iter = 0; order2Iter < ( ( bandIter * 2 ) + 1 ); order2Iter++ )
    {
        //============================================ Reset loop
        objCombValsIter                               = 0;
        
        //============================================ Find the c*conj(c) values for different radii
        for ( proshade_unsign radiusIter = 0; radiusIter < std::min( obj1->getMaxSpheres(), obj2->getMaxSpheres() ); radiusIter++ )
        {
    
            //======================================== Get only values where the shell has the band
            if ( std::min ( obj1->getShellBandwidth ( radiusIter ), obj2->getShellBandwidth ( radiusIter ) ) <= bandIter ) { continue; }
            
            //======================================== Multiply coeffs
            ProSHADE_internal_maths::complexMultiplicationConjug ( obj1->getRealSphHarmValue ( bandIter, orderIter,  radiusIter ),
                                                                   obj1->getImagSphHarmValue ( bandIter, orderIter,  radiusIter ),
                                                                   obj2->getRealSphHarmValue ( bandIter, order2Iter, radiusIter ),
                                                                   obj2->getImagSphHarmValue ( bandIter, order2Iter, radiusIter ),
                                                                  &hlpReal, &hlpImag );
  
            //======================================== Apply r^2 integral weight
            radiiVals[objCombValsIter][0]             = hlpReal *  pow ( ( static_cast<proshade_double> ( obj1->getAnySphereRadius( radiusIter ) ) ), 2.0 );
            radiiVals[objCombValsIter][1]             = hlpImag *  pow ( ( static_cast<proshade_double> ( obj1->getAnySphereRadius( radiusIter ) ) ), 2.0 );
        
            objCombValsIter                          += 1;
        }
        
        //============================================ Integrate over all radii using n-point Gauss-Legendre integration
        ProSHADE_internal_maths::gaussLegendreIntegration ( radiiVals, objCombValsIter, integOrder, abscissas, weights, integRange, sphereDist, &hlpReal, &hlpImag );
  
        //============================================ Save the result into E matrices
        arrVal[0]                                     = hlpReal;
        arrVal[1]                                     = hlpImag;
        obj2->setEMatrixValue                         ( bandIter, orderIter, order2Iter, arrVal );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the E matrix weight values for a given band and order and saves these into the appropriate objects.
 
    \param[in] obj1 The ProSHADE_data object for which the comparison is done in regards to.
    \param[in] obj2 The ProSHADE_data object for which the comparison is done in regards from - the E matrices will be saved into this object.
    \param[in] band The bandwidth of the SH value for which this should be done.
    \param[in] order The order of the SH value for which this should be done.
    \param[in] obj1Vals Already allocated array of proshade_double to which integrated values will be saved. It must have size equal to minimum of spheres in the two compared objects.
    \param[in] obj2Vals Already allocated array of proshade_double to which integrated values will be saved. It must have size equal to minimum of spheres in the two compared objects.
    \param[in] integOrder The Gauss-Legendre integration order to be used.
    \param[in] abscissas The pre-computed abscissas for the Gauss-Legendre integration.
    \param[in] weights The pre-computed weights for the Gauss-Legendre integration.
    \param[in] integRange The range in angstroms between the smalleds and largest shell which are integrated over (might not be 0 to max for progressive shell sampling).
    \param[in] sphereDist The distance between any two spheres.
    \param[out] sphereRange The distance between the smallest and largest usable sphere (usable as in having the required band).
 */
proshade_double ProSHADE_internal_distances::computeWeightsForEMatricesForLM ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2, proshade_unsign bandIter, proshade_unsign orderIter, proshade_double* obj1Vals, proshade_double* obj2Vals, proshade_unsign integOrder, proshade_double* abscissas, proshade_double* weights, proshade_single sphereDist )
{
    //================================================ Initialise local values
    proshade_unsign obj1ValsIter                      = 0;
    proshade_unsign obj2ValsIter                      = 0;
    
    //================================================ Set sphere counters
    proshade_unsign minSphere                         = std::min( obj1->getMaxSpheres(), obj2->getMaxSpheres() );
    proshade_unsign maxSphere                         = 0;
    
    //================================================ For each radius, deal with weights
    for ( proshade_unsign radiusIter  = 0; radiusIter < std::min( obj1->getMaxSpheres(), obj2->getMaxSpheres() ); radiusIter++ )
    {
        //============================================ Get only values where the shell has the band
        if ( std::min ( obj1->getShellBandwidth ( radiusIter ), obj2->getShellBandwidth ( radiusIter ) ) <= bandIter ) { continue; }
        minSphere                                     = std::min ( radiusIter, minSphere );
        maxSphere                                     = std::max ( radiusIter, maxSphere );
        
        //============================================ Get the magnitudes for weighting
        computeSphericalHarmonicsMagnitude            ( obj1, bandIter, orderIter, radiusIter, &(obj1Vals[obj1ValsIter]) );
        computeSphericalHarmonicsMagnitude            ( obj2, bandIter, orderIter, radiusIter, &(obj2Vals[obj2ValsIter]) );
        obj1ValsIter                                 += 1;
        obj2ValsIter                                 += 1;
    }
    
    //================================================ Integrate weights
    proshade_single minSphereRad                      = obj1->getSpherePosValue ( minSphere ) - ( sphereDist * 0.5f );
    proshade_single maxSphereRad                      = obj1->getSpherePosValue ( maxSphere ) + ( sphereDist * 0.5f );
            
    obj1->setIntegrationWeightCumul                   ( ProSHADE_internal_maths::gaussLegendreIntegrationReal ( obj1Vals, obj1ValsIter, integOrder, abscissas, weights, static_cast< proshade_double > ( maxSphereRad - minSphereRad ), static_cast< proshade_double > ( sphereDist ) ) );
    obj2->setIntegrationWeightCumul                   ( ProSHADE_internal_maths::gaussLegendreIntegrationReal ( obj2Vals, obj2ValsIter, integOrder, abscissas, weights, static_cast< proshade_double > ( maxSphereRad - minSphereRad ), static_cast< proshade_double > ( sphereDist ) ) );
    
    //================================================ Done
    return                                            ( static_cast< proshade_double > ( maxSphereRad - minSphereRad ) );
    
}

/*! \brief This helper function is responsible for deleting the workspace memory required for trace sigma descriptor computation.
 
    \param[in] obj1Vals Array to hold the shell values for the first object integgration.
    \param[in] obj2Vals Array to hold the shell values for the second object integgration.
    \param[in] GLabscissas An array to hold the pre-computed anscissas for the Gauss-Legendre integration.
    \param[in] glWeights An array to hold the pre-computed weights for the Gauss-Legendre integration.
    \param[in] radiiVals A complex array to hold the results of combining spherical harmonics coefficients of the two objects for each shell.
 */
void ProSHADE_internal_distances::releaseTrSigmaWorkspace ( proshade_double*& obj1Vals, proshade_double*& obj2Vals, proshade_double*& GLabscissas, proshade_double*& GLweights, proshade_complex*& radiiVals )
{
    //================================================ Release memory
    delete[] obj1Vals;
    delete[] obj2Vals;
    delete[] radiiVals;
    delete[] GLabscissas;
    delete[] GLweights;
    
    //================================================ Set to NULL
    obj1Vals                                          = nullptr;
    obj2Vals                                          = nullptr;
    radiiVals                                         = nullptr;
    GLabscissas                                       = nullptr;
    GLweights                                         = nullptr;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the complete E matrices and their weights between any two objects.
 
    This function allocates the space required for storing the E matrices, allocates all the workspace requierd for the computation
    and proceeds to compute the values for all band (l), order1(m) and order2(m') E matrix values. It then proceeds to release all
    non required memory and terminates, leaving all its results in the second ProSHADE data object supplied. This function does NOT
    apply the weights to the matrices, it needs to be done subsequently!
 
    \param[in] obj1 The first ProSHADE_data object for which the computation is done.
    \param[in] obj2 The second ProSHADE_data object for which the computation is done.
    \param[in] settings A pointer to settings class containing all the information required for the task.
 */
void ProSHADE_internal_distances::computeEMatrices    ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2, ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Starting computation of E matrices.", settings->messageShift );
    
    //================================================ Allocatre memory for E matrices in the second object (first may be compared to more structures and therefore its data would be written over)
    obj2->allocateEMatrices                           ( std::min ( obj1->getMaxBand(), obj2->getMaxBand() ) );
    
    //================================================ Initialise local variables
    proshade_double *obj1Vals, *obj2Vals, *GLAbscissas, *GLWeights;
    proshade_complex* radiiVals;
    proshade_double integRange;
    
    //================================================ Allocate workspace memory
    allocateTrSigmaWorkspace                          ( std::min( obj1->getMaxSpheres(), obj2->getMaxSpheres() ), settings->integOrder, obj1Vals, obj2Vals, GLAbscissas, GLWeights,  radiiVals);
    
    //================================================ Initialise abscissas and weights for integration
    ProSHADE_internal_maths::getLegendreAbscAndWeights ( settings->integOrder, GLAbscissas, GLWeights, settings->taylorSeriesCap );
    
    //================================================ For each band (l), compute the E matrix integrals
    for ( proshade_unsign bandIter = 0; bandIter < std::min ( obj1->getMaxBand(), obj2->getMaxBand() ); bandIter++ )
    {
        //============================================ For each order (m)
        for ( proshade_unsign orderIter = 0; orderIter < ( ( bandIter * 2 ) + 1 ); orderIter++ )
        {
            //======================================== Get weights for the required band(l) and order (m)
            integRange                                = computeWeightsForEMatricesForLM ( obj1, obj2, bandIter, orderIter, obj1Vals, obj2Vals, settings->integOrder, GLAbscissas, GLWeights, settings->maxSphereDists );

            //======================================== Compute E matrices value for given band (l) and order(m)
            computeEMatricesForLM                     ( obj1, obj2, bandIter, orderIter, radiiVals, settings->integOrder, GLAbscissas, GLWeights, integRange, static_cast< proshade_double > ( settings->maxSphereDists ) );
        }
        
        //============================================ Report progress
        if ( settings->verbose > 3 )
        {
            std::stringstream hlpSS;
            hlpSS << "E matrices computed for band " << bandIter;
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 4, hlpSS.str(), settings->messageShift );
        }
    }
    
    //================================================ Release the workspace memory
    releaseTrSigmaWorkspace                           ( obj1Vals, obj2Vals, GLAbscissas, GLWeights, radiiVals );
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, "E matrices computed.", settings->messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function normalises the E matrices.
 
    This function assumes that the E matrices and integration magnitude weights were already computed. It now proceeds to compute the weighting
    factor (sqrt of the product of the magnitudes of the two objects) and apply it to the E matrices. This normalisation is similar in formula
    and meaning to the Pearson's correlation coefficient normalisation.
 
    \param[in] obj1 The first ProSHADE_data object for which the computation is done.
    \param[in] obj2 The second ProSHADE_data object for which the computation is done.
    \param[in] settings A pointer to settings class containing all the information required for the task.
 */
void ProSHADE_internal_distances::normaliseEMatrices ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2, ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, "Starting E matrices normalisation.", settings->messageShift );
    
    //================================================ Normalise by the Pearson's c.c. like formula
    proshade_double eMatNormFactor                    = std::sqrt ( obj1->getIntegrationWeight() * obj2->getIntegrationWeight() );
    
    //================================================ If this is self-correlation (i.e. obj1 == obj2), then divide normalisation factor by 2 as the weight was applied cumulatively!
    if ( obj1->inputOrder == obj2->inputOrder ) { eMatNormFactor /= 2.0; }
    
    for ( proshade_unsign bandIter = 0; bandIter < std::min ( obj1->getMaxBand(), obj2->getMaxBand() ); bandIter++ )
    {
        //============================================ For each combination of m and m' for E matrices
        for ( proshade_unsign orderIter = 0; orderIter < ( ( bandIter * 2 ) + 1 ); orderIter++ )
        {
            for ( proshade_unsign order2Iter = 0; order2Iter < ( ( bandIter * 2 ) + 1 ); order2Iter++ )
            {
                obj2->normaliseEMatrixValue           ( bandIter, orderIter, order2Iter, eMatNormFactor );
            }
        }
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 4, "E matrices normalised.", settings->messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the trace sigma descriptor value between two objects.
 
    This function starts by checking if the trace sigma descriptor was requested and if so, proceeds to compute the E matrices.
    These are 3D matrices with each l,m,m' value being the combination of the c_{l,m} and c*_{l,m'} spherical harmonics coefficients.
    Once computed, the E matrices are normalised by the magnitudes of the objects spherical harmonics coefficients and the SVD is computed
    for each l (i.e. on each m x m' matrix). The sum of the trace of the sigmas of the SVD is then the trace sigma descriptor, whose value
    is returned.
 
    \param[in] obj1 The first ProSHADE_data object against which comparison is done.
    \param[in] obj2 The second ProSHADE_data object which is compared to the first.
    \param[in] settings A pointer to settings class containing all the information required for the task.
    \param[out] ret The final normalised value of the trace sigma descriptor for the two objects.
 */
proshade_double ProSHADE_internal_distances::computeTraceSigmaDescriptor ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2, ProSHADE_settings* settings )
{
    //================================================ Report starting the task
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting trace sigma distance computation.", settings->messageShift );
    
    //================================================ Initialise return variable
    proshade_double ret                               = 0.0;
    
    //================================================ Sanity check
    if ( !settings->computeTraceSigmaDesc )
    {
        throw ProSHADE_exception ( "Attempted computing trace sigma descriptors when it was\n                    : not required.", "ED00018", __FILE__, __LINE__, __func__, "Attempted to pre-compute the E matrices, when the user\n                    : has specifically stated that these should not be computed.\n                    : Unless you manipulated the code, this error should never\n                    : occur; if you see this, I made a large blunder. Please let\n                    : me know!" );
    }
    
    //================================================ Empty the cumulative weights back to 0.0 for each structure
    obj1->setIntegrationWeight                        ( 0.0 );
    obj1->setIntegrationWeight                        ( 0.0 );

    //================================================ Compute un-weighted E matrices and their weights
    computeEMatrices                                  ( obj1, obj2, settings );
    
    //================================================ Normalise E matrices by the magnitudes
    normaliseEMatrices                                ( obj1, obj2, settings );
    
    //================================================ Allocate the required memory
    double* singularValues                            = new double[( ( std::min ( obj1->getMaxBand(), obj2->getMaxBand() ) * 2 ) + 1 )];
    ProSHADE_internal_misc::checkMemoryAllocation     ( singularValues, __FILE__, __LINE__, __func__ );
    
    //================================================ Compute the distance
    for ( proshade_unsign lIter = 0; lIter < std::min ( obj1->getMaxBand(), obj2->getMaxBand() ); lIter++ )
    {
        //============================================ Find the complex matrix SVD singular values
        ProSHADE_internal_maths::complexMatrixSVDSigmasOnly ( obj2->getEMatrixByBand ( lIter ), static_cast<int> ( ( lIter * 2 ) + 1 ), singularValues );

        //============================================ Now sum the trace
        for ( proshade_unsign iter = 0; iter < ( ( lIter * 2 ) + 1 ); iter++  )
        {
            ret                                      += singularValues[iter];
        }
    }
    
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, "E matrices decomposed to singular values.", settings->messageShift );
    
    //================================================ Release the memory
    delete[] singularValues;
    
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Trace sigma distance computation complete.", settings->messageShift );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function allocates the memory for the SO(3) coefficients and the inverse for that calling object.
 
    \param[in] band The bandwidth to which the computation will be done.
 */
void ProSHADE_internal_data::ProSHADE_data::allocateSO3CoeffsSpace ( proshade_unsign band )
{
    //================================================ Allocate the memory
    this->so3Coeffs                                   = new fftw_complex [static_cast<proshade_unsign>( ( 4 * pow( static_cast<proshade_double> ( band ), 3.0 ) - static_cast<proshade_double> ( band ) ) / 3.0 )];
    this->so3CoeffsInverse                            = new fftw_complex [static_cast<proshade_unsign>( pow( static_cast<proshade_double> ( band ) * 2.0, 3.0 ) )];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->so3Coeffs,        __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->so3CoeffsInverse, __FILE__, __LINE__, __func__ );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function converts the E matrices to SO(3) coefficients.
 
    This function starts by allocating the memory for the SO(3) coefficients and their inverse. It then
    proceeds to convert the E matrix values into the SO(3) transform coefficients by applying the Wigner
    normalisation factor and changing the sign as required by SOFT library. Upon termination, the coeffs
    will be saved in the obj2 class.
 
    \param[in] obj1 The first ProSHADE_data object against which comparison is done.
    \param[in] obj2 The second ProSHADE_data object which is compared to the first.
    \param[in] settings A pointer to settings class containing all the information required for the task.
 */
void ProSHADE_internal_distances::generateSO3CoeffsFromEMatrices ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2, ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Converting E matrices to SO(3) coefficients.", settings->messageShift );
    
    //================================================ Allocate memory for the coefficients
    obj2->allocateSO3CoeffsSpace                      ( std::min ( obj1->getMaxBand(), obj2->getMaxBand() ) );
    
    //================================================ Initialise local variables
    proshade_double wigNorm, hlpValReal, hlpValImag;
    proshade_double signValue                         = 1.0;
    proshade_unsign indexO;
    proshade_complex hlpVal;
    
    //================================================ For each band (l)
    for ( proshade_signed bandIter = 0; bandIter < static_cast<proshade_signed> ( std::min ( obj1->getMaxBand(), obj2->getMaxBand() ) ); bandIter++ )
    {
        //============================================ Get wigner normalisation factor
        wigNorm                                       = 2.0 * M_PI * sqrt ( 2.0 / (2.0 * static_cast< proshade_double > ( bandIter ) + 1.0 ) );
        
        //============================================ For each order (m)
        for ( proshade_signed orderIter = 0; orderIter < ( ( bandIter * 2 ) + 1 ); orderIter++ )
        {
            //======================================== Set the sign
            if ( ( orderIter - bandIter + bandIter ) % 2 ) { signValue = -1.0 ; }
            else                                           { signValue =  1.0 ; }
            
            //======================================== For each order2 (m')
            for ( proshade_signed order2Iter = 0; order2Iter < ( ( bandIter * 2 ) + 1 ); order2Iter++ )
            {
                //==================================== Find output index
                indexO                                = static_cast< proshade_unsign > ( so3CoefLoc ( static_cast< int > ( orderIter - bandIter ), static_cast< int > ( order2Iter - bandIter ), static_cast< int > ( bandIter ), static_cast< int > ( std::min ( obj1->getMaxBand(), obj2->getMaxBand() ) ) ) );
                
                //==================================== Compute and save the SO(3) coefficients
                obj2->getEMatrixValue                 ( static_cast< proshade_unsign > ( bandIter ), static_cast< proshade_unsign > ( orderIter ), static_cast< proshade_unsign > ( order2Iter ), &hlpValReal, &hlpValImag );
                hlpVal[0]                             = hlpValReal * wigNorm * signValue;
                hlpVal[1]                             = hlpValImag * wigNorm * signValue;
                obj2->setSO3CoeffValue                ( indexO, hlpVal );
                
                //==================================== Switch the sign value
                signValue                            *= -1.0;
            }
        }
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, "SO(3) coefficients obtained.", settings->messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allocates the workspaces required to compute the inverse SOFT transform.
 
    \param[in] work1 The first workspace pointer to be allocated.
    \param[in] work2 The second workspace pointer to be allocated.
    \param[in] work3 The third workspace pointer to be allocated.
    \param[in] band The bandwidth of the computations (this determines how much workspace will be required).
 */
void ProSHADE_internal_distances::allocateInvSOFTWorkspaces ( proshade_complex*& work1, proshade_complex*& work2, proshade_double*& work3, proshade_unsign band )
{
    //================================================ Allocate memory
    work1                                             = new proshade_complex[8  * static_cast<proshade_unsign> ( pow( static_cast<double> ( band ), 3.0 ) )];
    work2                                             = new proshade_complex[14 * static_cast<proshade_unsign> ( pow( static_cast<double> ( band ), 2.0 ) ) + (48 * band)];
    work3                                             = new proshade_double [2  * static_cast<proshade_unsign> ( pow( static_cast<double> ( band ), 2.0 ) ) + (24 * band)];
    
    //================================================ Check the memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( work1, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( work2, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( work3, __FILE__, __LINE__, __func__ );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function prepares the FFTW plan for the inverse SO(3) transform.
 
    \param[in] inverseSO3 The FFTW_PLAN pointer where the result will be saved.
    \param[in] band The bandwidth of the computations.
    \param[in] work1 The workspace to be used for the computation.
    \param[in] invCoeffs The pointer to where the inverse SOFT transform results will be saved.
 */
void ProSHADE_internal_distances::prepareInvSOFTPlan ( fftw_plan* inverseSO3, int band, fftw_complex* work1, proshade_complex* invCoeffs )
{
    //================================================ Prepare the plan describing variables
    int howmany                                       = 4 * band * band;
    int idist                                         = 2 * band;
    int odist                                         = 2 * band;
    int rank                                          = 2;
            
    int inembed[2], onembed[2];
    inembed[0]                                        = 2 * band;
    inembed[1]                                        = 4 * band * band;
    onembed[0]                                        = 2 * band;
    onembed[1]                                        = 4 * band * band;
            
    int istride                                       = 1;
    int ostride                                       = 1;
            
    int na[2];
    na[0]                                             = 1;
    na[1]                                             = 2 * band;
    
    //================================================ Create the plan
   *inverseSO3                                        = fftw_plan_many_dft ( rank,
                                                                             na,
                                                                             howmany,
                                                                             work1,
                                                                             inembed,
                                                                             istride,
                                                                             idist,
                                                                             invCoeffs,
                                                                             onembed,
                                                                             ostride,
                                                                             odist,
                                                                             FFTW_FORWARD,
                                                                             FFTW_ESTIMATE );
            
    //================================================ Done
    return ;
    
}

/*! \brief This function releases the memory used for computation of the inverse SOFT transform.
 
    \param[in] work1 The first workspace pointer to be released.
    \param[in] work2 The second workspace pointer to be released.
    \param[in] work3 The third workspace pointer to be released.
 */
void ProSHADE_internal_distances::releaseInvSOFTMemory ( proshade_complex*& work1, proshade_complex*& work2, proshade_double*& work3 )
{
    //================================================ Release memory
    delete[] work1;
    delete[] work2;
    delete[] work3;
    
    //================================================ Done
    return ;
}

/*! \brief This function computes the inverse SO(3) transform.
 
    This function firstly allocates all the required workspaces for the inverse SO(3) Fourier Transform, then it
    prepares the FFTW plans for performing the FFTW inverse Fourier transform in the SO(3) space using FFTW and
    finally it subjects the SO(3) coeffficients available at this point to the computation. The results are saved
    into the second object, memory is released and function terminates.
 
    \param[in] obj1 The first ProSHADE_data object against which comparison is done.
    \param[in] obj2 The second ProSHADE_data object which is compared to the first.
    \param[in] settings A pointer to settings class containing all the information required for the task.
 */
void ProSHADE_internal_distances::computeInverseSOFTTransform ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2, ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Computing inverse SO(3) Fourier transform.", settings->messageShift );
    
    //================================================ Initialise local variables
    proshade_complex *workspace1, *workspace2;
    proshade_double *workspace3;
    fftw_plan inverseSO3;
    
    //================================================ Allocate memory for the workspaces
    allocateInvSOFTWorkspaces                         ( workspace1, workspace2, workspace3, std::min ( obj1->getMaxBand(), obj2->getMaxBand() ) );
    
    //================================================ Prepare the FFTW plan
    prepareInvSOFTPlan                                ( &inverseSO3, static_cast< int > ( std::min ( obj1->getMaxBand(), obj2->getMaxBand() ) ), workspace1, obj2->getInvSO3Coeffs ( ) );
    
    //================================================ Compute the transform
    Inverse_SO3_Naive_fftw                            ( static_cast< int > ( std::min ( obj1->getMaxBand(), obj2->getMaxBand() ) ),
                                                        obj2->getSO3Coeffs ( ),
                                                        obj2->getInvSO3Coeffs ( ),
                                                        workspace1,
                                                        workspace2,
                                                        workspace3,
                                                       &inverseSO3,
                                                        0 );
    
    //================================================ Release memory
    releaseInvSOFTMemory                              ( workspace1, workspace2, workspace3 );
    fftw_destroy_plan                                 ( inverseSO3 );
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, "Inverse SO(3) Fourier transform computed.", settings->messageShift );
    
    
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the rotation function descriptor value between two objects.
 
    This function starts by sanity checks and computation and normalisation of the E matrices used for the trace sigma
    descriptor (these only need to be computed if trace sigma descriptor is NOT computed, otherwise the stored values are
    used to save time). It then converts the E matrix values to SO(3) transform coefficients and invers these using the
    SOFT library. From the resulting rotation function map, it selects the highest peak and applies its rotation to the
    E matrix values. The resulting values are simply summed, as this sum can be proven to be the argmin of the distance
    between the two objects in their spherical harmonics decomposition space.
 
    \param[in] obj1 The first ProSHADE_data object against which comparison is done.
    \param[in] obj2 The second ProSHADE_data object which is compared to the first.
    \param[in] settings A pointer to settings class containing all the information required for the task.
    \param[out] ret The final normalised value of the rotation function descriptor for the two objects.
 */
proshade_double ProSHADE_internal_distances::computeRotationFunctionDescriptor ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2, ProSHADE_settings* settings )
{
    //================================================ Report starting the task
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting rotation function distance computation.", settings->messageShift );
    
    //================================================ Initialise return variable
    proshade_double ret                               = 0.0;
    proshade_double eulA, eulB, eulG, EMatR, EMatI, WigDR, WigDI;
    
    //================================================ Sanity check
    if ( !settings->computeRotationFuncDesc )
    {
        throw ProSHADE_exception ( "Attempted computing rotation function descriptors when it\n                    : was not required.", "ED00023", __FILE__, __LINE__, __func__, "Attempted to compute the SO(3) transform and the rotation \n                    : function descriptor when the user did not request this. \n                    : Unless you manipulated the code, this error should never \n                    : occur; if you see this, I made a large blunder. \n                    : Please let me know!" );
    }
    
    //================================================ Compute weighted E matrices if not already present
    if ( !settings->computeTraceSigmaDesc )
    {
        computeEMatrices                              ( obj1, obj2, settings );
        normaliseEMatrices                            ( obj1, obj2, settings );
    }
    
    //================================================ Generate SO(3) coefficients
    generateSO3CoeffsFromEMatrices                    ( obj1, obj2, settings );
    
    //================================================ Compute the inverse SO(3) Fourier Transform (SOFT) on the newly computed coefficients
    computeInverseSOFTTransform                       ( obj1, obj2, settings );

    //================================================ Get inverse SO(3) map top peak Euler angle values
    ProSHADE_internal_peakSearch::getBestPeakEulerAngsNaive ( obj2->getInvSO3Coeffs (),
                                                              std::min ( obj1->getMaxBand(), obj2->getMaxBand() ) * 2,
                                                             &eulA, &eulB, &eulG, settings );

    //================================================ Compute the Wigner D matrices for the Euler angles
    ProSHADE_internal_wigner::computeWignerMatricesForRotation ( settings, obj2, eulA, eulB, eulG );

    //================================================ Compute the distance
    for ( proshade_unsign bandIter = 0; bandIter < obj2->getComparisonBand(); bandIter++ )
    {
        //============================================ For each order1
        for ( proshade_unsign order1 = 0; order1 < ( ( bandIter * 2 ) + 1 ); order1++ )
        {
            //======================================== For each order2
            for ( proshade_unsign order2 = 0; order2 < ( ( bandIter * 2 ) + 1 ); order2++ )
            {
                //==================================== Multiply D_{l} * E_{l} and get sum over l of traces (i.e. just sum all together)
                obj2->getEMatrixValue                 ( bandIter, order1, order2, &EMatR, &EMatI );
                obj2->getWignerMatrixValue            ( bandIter, order2, order1, &WigDR, &WigDI );
                ret                                  += ProSHADE_internal_maths::complexMultiplicationRealOnly ( &WigDR, &WigDI, &EMatR, &EMatI );
            }
        }
    }
    
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Rotation function distance computation complete.", settings->messageShift );
    
    //================================================ Done
    return                                            ( ret );
    
}
