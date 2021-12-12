/*! \file ProSHADE_sphericalHarmonics.cpp
    \brief This source file contains the function required to compute the spherical harmonics decompostion in ProSHADE.
 
    The functions in this source file are all required to run the spherical harmonic decomposition of a ProSHADE_sphere objet. They make use of the SOFT2.0 library in order to get the
    compuation done. This file also contains the support, memory allocation and similar functions, all related to the main task.
 
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
#include "ProSHADE_sphericalHarmonics.hpp"

/*! \brief This function determines the integration order for the between spheres integration.
 
    This function simply takes all pointer variables required for the spherical harmonics computation and allocates the required
    amount of memory for them. It also does the memory checks in case memory allocation fails.
 
    \param[in] band The bandwidth to which the computation will be done.
    \param[in] inputReal The real part of the input will be copied here.
    \param[in] inputImag The immaginary part of the input will be copied here.
    \param[in] outputReal The real part of the output will be saved here.
    \param[in] outputImag The immaginary part of the output will be saved here.
    \param[in] shWeights The weights for spherical harmonics computation will be stored here.
    \param[in] tableSpaceHelper This space is required by SOFT for pre-computing values into this table.
    \param[in] workspace The space where multiple minor results are saved by SOFT.
 */
void ProSHADE_internal_sphericalHarmonics::allocateComputationMemory ( proshade_unsign band, proshade_double*& inputReal, proshade_double*& inputImag, proshade_double*& outputReal, proshade_double*& outputImag, proshade_double*& shWeights, proshade_double*& tableSpaceHelper, fftw_complex*& workspace )
{
    //================================================ Initialise local variables
    proshade_unsign oneDimmension                     = 2 * band;
    
    //================================================ Allocate Input Memory
    inputReal                                         = new proshade_double [oneDimmension * oneDimmension];
    inputImag                                         = new proshade_double [oneDimmension * oneDimmension];
    
    //================================================ Allocate Output Memory
    outputReal                                        = new proshade_double [oneDimmension * oneDimmension];
    outputImag                                        = new proshade_double [oneDimmension * oneDimmension];
    
    //================================================ Allocate Working Memory
    shWeights                                         = new proshade_double [band * 4];
    workspace                                         = new fftw_complex    [(  8 * band * band ) +  ( 10 * band )];
    
    //================================================ Allocate table
    tableSpaceHelper                                  = new proshade_double [static_cast<proshade_unsign> ( Reduced_Naive_TableSize ( static_cast< int > ( band ), static_cast< int > ( band ) ) +
                                                                                                            Reduced_SpharmonicTableSize ( static_cast< int > ( band ), static_cast< int > ( band ) ) )];
    
    //================================================ Check memory allocation success
    ProSHADE_internal_misc::checkMemoryAllocation     ( inputReal,        __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( inputImag,        __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( outputReal,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( outputImag,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( shWeights,        __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( tableSpaceHelper, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( workspace,        __FILE__, __LINE__, __func__ );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes the workspace pointer and correctly places the other internal pointers.
 
    This is a simple helper function, which places the internal workspace pointers to the correct addresses as
    required by the SOFT2.0 dependency.
 
    \param[in] workspace Pointer to the allocated workspace, within which the other pointers will be placed.
    \param[in] oDim The size of the single transform dimension (twice the transform bandwidth).
    \param[in] rres Pointer to where the real part of the results will be temporarily saved.
    \param[in] ires Pointer to where the imaginary part of the results will be temporarily saved.
    \param[in] fltres Pointer to where the temporary bandwise results will be saved.
    \param[in] scratchpad Pointer to extra space which is internally used (but not allocated) by SOFT2.0.
 */
void ProSHADE_internal_sphericalHarmonics::placeWithinWorkspacePointers ( fftw_complex*& workspace, proshade_unsign oDim, proshade_double*& rres, proshade_double*& ires, proshade_double*& fltres, proshade_double*& scratchpad )
{
    //================================================ Place pointers as required by SOFT2.0
    rres                                              = reinterpret_cast<proshade_double*> ( workspace );
    ires                                              = rres + ( oDim * oDim );
    fltres                                            = ires + ( oDim * oDim );
    scratchpad                                        = fltres + ( oDim / 2 );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function initialises the FFTW plans.
 
    This function initialises the FFTW plans for the spherical harmonics computations as required by the SOFT2.0
    library.
 
    \param[in] band The bandwidth to which the computation will be done.
    \param[in] fftPlan pointer to the variable where the Fourier transform should be set.
    \param[in] dctPlan pointer to the variable where the 1D r2r Fourier transform should be set.
    \param[in] inputReal pointer to the array containing (or which will contain) the input real values.
    \param[in] inputImag pointer to the array containing (or which will contain) the input imaginary values.
    \param[in] rres pointer to the array where the real values of result should be saved.
    \param[in] ires pointer to the array where the imaginary values of result should be saved.
    \param[in] scratchpad pointer to the array where temporary results will be saved.
 */
void ProSHADE_internal_sphericalHarmonics::initialiseFFTWPlans ( proshade_unsign band, fftw_plan& fftPlan, fftw_plan& dctPlan, proshade_double*& inputReal, proshade_double*& inputImag, proshade_double*& rres, proshade_double*& ires, proshade_double*& scratchpad )
{
    //================================================ Initialize fft plan along phi angles
    fftw_iodim dims[1];
    fftw_iodim howmany_dims[1];
    
    int rank                                          = 1;
    int howmany_rank                                  = 1;
            
    dims[0].n                                         = static_cast<int> ( band * 2 );
    dims[0].is                                        = 1;
    dims[0].os                                        = static_cast<int> ( band * 2 );
            
    howmany_dims[0].n                                 = static_cast<int> ( band * 2 );
    howmany_dims[0].is                                = static_cast<int> ( band * 2 );
    howmany_dims[0].os                                = 1;
    
    //================================================ Plan fft transform
    fftPlan                                           = fftw_plan_guru_split_dft ( rank,
                                                                                   dims,
                                                                                   howmany_rank,
                                                                                   howmany_dims,
                                                                                   inputReal,
                                                                                   inputImag,
                                                                                   rres,
                                                                                   ires,
                                                                                   FFTW_ESTIMATE  );
    
    //================================================ Initialize dct plan for SHT
    dctPlan                                           = fftw_plan_r2r_1d ( static_cast<int> ( band * 2 ),
                                                                           scratchpad,
                                                                           scratchpad + static_cast<int> ( band * 2 ),
                                                                           FFTW_REDFT10,
                                                                           FFTW_ESTIMATE ) ;
    
    //================================================ Done
    return ;
    
}


/*! \brief This function computes the spherical harmonics of a aingle shell, saving them in supplied pointer.
 
    This function takes all the pointers required for the spherical harmonics computation and releases all the
    memory.
 
    \param[in] inputReal pointer to the array that contained the input real values to be freed.
    \param[in] inputImag pointer to the array that contained the input imaginary values to be freed.
    \param[in] outputReal pointer to the array that contained the output real values to be freed.
    \param[in] outputImag pointer to the array that contained the output imaginary values to be freed.
    \param[in] tableSpaceHelper pointer to the helper array for Legendre polynomials values to be freed.
    \param[in] tableSpace pointer to the array of Legendre polynomials values to be freed.
    \param[in] shWeights pointer to the array spherical harmonics weighhts to be freed.
    \param[in] workspace pointer to the array for miscellaneous temporary results to be freed.
    \param[in] fftPlan pointer to the variable where the Fourier transform was done to be freed.
    \param[in] dctPlan pointer to the variable where the 1D r2r Fourier transform was done to be freed.
 */
void ProSHADE_internal_sphericalHarmonics::releaseSphericalMemory ( proshade_double*& inputReal, proshade_double*& inputImag, proshade_double*& outputReal, proshade_double*& outputImag, double*& tableSpaceHelper, double**& tableSpace, double*& shWeights, fftw_complex*& workspace, fftw_plan& fftPlan, fftw_plan& dctPlan )
{
    //================================================ Release all memory related to SH
    delete[] inputReal;
    delete[] inputImag;
    delete[] outputReal;
    delete[] outputImag;
    delete[] tableSpaceHelper;
    delete[] shWeights;
    delete[] workspace;
            
    //================================================ Set pointers to NULL
    tableSpaceHelper                                  = nullptr;
    tableSpace                                        = nullptr;
    shWeights                                         = nullptr;
    workspace                                         = nullptr;
          
    //================================================ Delete fftw plans
    fftw_destroy_plan                                 ( dctPlan );
    fftw_destroy_plan                                 ( fftPlan );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function initialises all the memory required for spherical harmonics computation.
 
    This function takes on all the memory allocation and filling in all data required later for the
    spherical harmonics computation by the SOFT2.0 library.
 
    \param[in] band The bandwidth to which the computation will be done.
    \param[in] inputReal The real part of the input will be copied here.
    \param[in] inputImag The immaginary part of the input will be copied here.
    \param[in] outputReal The real part of the output will be saved here.
    \param[in] outputImag The immaginary part of the output will be saved here.
    \param[in] shWeights The weights for spherical harmonics computation will be stored here.
    \param[in] tableSpace This space is required by SOFT for pre-computing values into this table.
    \param[in] tableSpaceHelper This space is required by SOFT for proper computation of the table space.
    \param[in] workspace The space where multiple minor results are saved by SOFT2.0.
    \param[in] rres Pointer to where the real part of the results will be temporarily saved.
    \param[in] ires Pointer to where the imaginary part of the results will be temporarily saved.
    \param[in] fltres Pointer to where the temporary bandwise results will be saved.
    \param[in] scratchpad Pointer to extra space which is internally used (but not allocated) by SOFT2.0.
    \param[in] fftPlan pointer to the variable where the Fourier transform should be set.
    \param[in] dctPlan pointer to the variable where the 1D r2r Fourier transform should be set.
 */
void ProSHADE_internal_sphericalHarmonics::initialiseAllMemory ( proshade_unsign band, proshade_double*& inputReal, proshade_double*& inputImag, proshade_double*& outputReal, proshade_double*& outputImag, double*& shWeights, double**& tableSpace, double*& tableSpaceHelper, fftw_complex*& workspace, proshade_double*& rres, proshade_double*& ires, proshade_double*& fltres, proshade_double*& scratchpad, fftw_plan& fftPlan, fftw_plan& dctPlan )
{
    //================================================ Initialise local variables
    proshade_unsign oneDim                            = band * 2;
    
    //================================================ Allocate memory for local pointers
    allocateComputationMemory                         ( band, inputReal, inputImag, outputReal, outputImag, shWeights, tableSpaceHelper, workspace );
    
    //================================================ Within workspace pointers
    placeWithinWorkspacePointers                      ( workspace, oneDim, rres, ires, fltres, scratchpad );
    
    //================================================ Generate Seminaive and naive tables for Legendre Polynomials
    tableSpace                                        = SemiNaive_Naive_Pml_Table ( static_cast< int > ( band ), static_cast< int > ( band ), tableSpaceHelper, reinterpret_cast<double*> ( workspace ) );
    
    //================================================ Make weights for spherical transform
    makeweights                                       ( static_cast< int > ( band ), shWeights );
    
    //================================================ Initialize FFTW Plans
    initialiseFFTWPlans                               ( band, fftPlan, dctPlan, inputReal, inputImag, rres, ires, scratchpad );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the spherical harmonics of a aingle shell, saving them in supplied pointer.
 
    This function takes the already initialised and prepared pointers and values and proceeds to load the data
    into the proper places and compute the split discrete Fourier transform, thus preparing for the spherical
    transform to be done.
 
    \param[in] oneDim This is the size of any dimension of the transform (2 * bandwidth).
    \param[in] inputReal Pointer to array which should be subjected to the transform (real part).
    \param[in] inputReal Pointer to array which should be subjected to the transform (imaginary part).
    \param[in] rres Pointer to array where the transform results will be saved (real part).
    \param[in] ires Pointer to array where the transform results will be saved (imaginary part).
    \param[in] mappedData Pointer to the data which should be decomposed.
    \param[in] fftPlan The prepared plan which states how the transform will be done as set by the initialiseFFTWPlans() function.
    \param[in] normCoeff The transform normalisation factor.
 */
void ProSHADE_internal_sphericalHarmonics::initialSplitDiscreteTransform ( proshade_unsign oneDim, proshade_double*& inputReal, proshade_double*& inputImag, proshade_double*& rres, proshade_double*& ires, proshade_double* mappedData, fftw_plan& fftPlan, proshade_double normCoeff )
{
    //================================================ Load mapped data to decomposition array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( oneDim * oneDim ); iter++ )
    {
        inputReal[iter]                               = mappedData[iter];
        inputImag[iter]                               = 0.0;
    }
    
    //================================================ Execute fft plan along phi
    fftw_execute_split_dft                            ( fftPlan, inputReal, inputImag, rres, ires ) ;
    
    //================================================ Normalize
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( oneDim * oneDim ); iter++ )
    {
        rres[iter]                                   *= normCoeff;
        ires[iter]                                   *= normCoeff;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes the split discrete transform and proceeds to complete the spherical harmonics decomposition.
 
    This function takes the results of the initial split discrete transform and proceeds to compute the spherical harmonics
    coefficients for all applicable bands using all the pre-computed valus (i.e. the Legendre polynomials table and the weights).
    To do this, the SOFT2.0 function is called and for more details, see SOFT2.0.
 
    \param[in] band The bandwidth to which the computation will be done.
    \param[in] rdataptr Pointer to be used as iterative locator in the large results array for real results.
    \param[in] idataptr Pointer to be used as iterative locator in the large results array for imaginary results.
    \param[in] outputReal An array for storing the real results of the completed transform.
    \param[in] outputImag An array for storing the imaginary results of the completed transform.
    \param[in] rres Array containing the real results of the initial split discrete transform.
    \param[in] ires Array containing the imaginary results of the initial split discrete transform.
    \param[in] fltres Helper array for transform computations.
    \param[in] scratchpad Array for keeping temporary results of the transform computations.
    \param[in] tablePml Pre-computed array of the Legendre polynomials as done by SOFT2.0.
    \param[in] shWeights The weights for the spherical harmonics.
    \param[in] dctPlan The FFTW plan for the final spherical harmonics transform.
 */
void ProSHADE_internal_sphericalHarmonics::computeSphericalTransformCoeffs ( proshade_unsign band, proshade_double*& rdataptr, proshade_double*& idataptr, proshade_double*& outputReal, proshade_double*& outputImag, proshade_double*& rres, proshade_double*& ires, proshade_double*& fltres, proshade_double*& scratchpad, double**& tablePml, double*& shWeights, fftw_plan& dctPlan )
{
    //================================================ Calculate the coefficients for each band
    rdataptr                                          = outputReal;
    idataptr                                          = outputImag;
    for ( proshade_unsign bandIter = 0; bandIter < band; bandIter++ )
    {
        //============================================ Real part calculation
        SemiNaiveReduced                              ( rres + ( bandIter * ( band * 2 ) ),
                                                        static_cast< int > ( band ),
                                                        static_cast< int > ( bandIter ),
                                                        fltres,
                                                        scratchpad,
                                                        tablePml[bandIter],
                                                        shWeights,
                                                       &dctPlan);
        
        //============================================ Save the real results to temporary holder
        memcpy                                        ( rdataptr, fltres, sizeof(proshade_double) * ( band - bandIter ) );
        rdataptr                                     += band - bandIter;
        
        //============================================ Imaginary part calculation
        SemiNaiveReduced                              ( ires + ( bandIter * ( band * 2 ) ),
                                                        static_cast< int > ( band ),
                                                        static_cast< int > ( bandIter ),
                                                        fltres,
                                                        scratchpad,
                                                        tablePml[bandIter],
                                                        shWeights,
                                                       &dctPlan );
        
        //============================================ Save the imaginary results
        memcpy                                        ( idataptr, fltres, sizeof(proshade_double) * ( band - bandIter ) );
        idataptr                                     += band - bandIter;
    }
    
    //================================================ DONE
    return ;
    
}

/*! \brief This is the final step in computing the full spherical harmonics decomposition of the input data.
 
    This is the final function that is needed for the complete spherical harmonics decomposition. It applied the
    Condon-Shortley phase as well as computing the negative orders, then saving the output into the final results
    array for further processing.
 
    \param[in] band The bandwidth to which the computation will be done.
    \param[in] outputReal The real results of the complete transform as done by the initialSplitDiscreteTransform() and computeSphericalTransformCoeffs() functions.
    \param[in] outputImag The imaginary results of the complete transform as done by the initialSplitDiscreteTransform() and computeSphericalTransformCoeffs() functions.
    \param[in] shArray An array of complex numbers to which the results of the spherical harmonics decomposition are to be saved.
 */
void ProSHADE_internal_sphericalHarmonics::applyCondonShortleyPhase ( proshade_unsign band, proshade_double* outputReal, proshade_double* outputImag, proshade_complex*& shArray )
{
    //================================================ Copy the results into the final holder
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( (band * 2) * (band * 2) ); iter++ )
    {
        shArray[iter][0]                              = outputReal[iter];
        shArray[iter][1]                              = outputImag[iter];
    }
    
    //================================================ Apply the Condon-Shortley phase sign
    proshade_double powerOne                          = 1.0;
    proshade_unsign hlp1                              = 0;
    proshade_unsign hlp2                              = 0;
    for ( proshade_signed order = 1; order < static_cast<proshade_signed> ( band ); order++)
    {
        powerOne                                     *= -1.0;
        for ( proshade_signed bandIter = order; bandIter < static_cast<proshade_signed> ( band ); bandIter++)
        {
            hlp1                                      = static_cast< proshade_unsign > ( seanindex ( static_cast< int > (  order ), static_cast< int > ( bandIter ), static_cast< int > ( band ) ) );
            hlp2                                      = static_cast< proshade_unsign > ( seanindex ( static_cast< int > ( -order ), static_cast< int > ( bandIter ), static_cast< int > ( band ) ) );
                    
            shArray[hlp2][0]                          =  powerOne * static_cast<proshade_double> ( outputReal[hlp1] );
            shArray[hlp2][1]                          = -powerOne * static_cast<proshade_double> ( outputImag[hlp1] );
        }
    }
    
    //================================================ DONE
    return ;
    
}

/*! \brief This function computes the spherical harmonics of a aingle shell, saving them in supplied pointer.
 
    This function does all the spherical harmonics computations for a single shell, including the memory allocation and
    releasing and the FFTW transforms. Because the shells can have different resolutions, the memory management is left
    until here, but possible speed-up could be gained from having the same resolution on all shells as in the older ProSHADE
    versions.
 
    \param[in] band The bandwidth to which the computation will be done.
    \param[in] sphereMappedData An array of doubles containing the mapped data onto a sphere for the sphere to be decomposed.
    \param[in] shArray An array of complex numbers to which the results of the spherical harmonics decomposition are to be saved.
 */
void ProSHADE_internal_sphericalHarmonics::computeSphericalHarmonics ( proshade_unsign band, proshade_double* sphereMappedData, proshade_complex*& shArray )
{
    //================================================ Initialise local variables
    proshade_double *inputReal = nullptr, *inputImag = nullptr, *outputReal = nullptr, *outputImag = nullptr;
    double *shWeights = nullptr, *tableSpaceHelper = nullptr;
    double** tablePml                                 = nullptr;
    fftw_complex* workspace                           = nullptr;
    proshade_unsign oneDim                            = static_cast<proshade_unsign> ( band * 2 );
    proshade_double normCoeff                         = ( 1.0 / ( static_cast<proshade_double> ( band * 2 ) ) ) * sqrt( 2.0 * M_PI );
    
    //================================================ Set output to zeroes (so that all unfilled data are not random)
    for ( proshade_unsign i = 0; i < ( 2 * band * 2 * band); i++ )
    {
        shArray[i][0]                                 = 0.0;
        shArray[i][1]                                 = 0.0;
    }
    
    //================================================ Within workspace pointers
    proshade_double *rres = nullptr, *ires = nullptr, *fltres = nullptr, *scratchpad = nullptr, *rdataptr = nullptr, *idataptr = nullptr;
    
    //================================================ FFTW Plans
    fftw_plan fftPlan                                 = nullptr;
    fftw_plan dctPlan                                 = nullptr;
    
    //================================================ Initialise all memory
    initialiseAllMemory                               ( band, inputReal, inputImag, outputReal, outputImag, shWeights, tablePml, tableSpaceHelper, workspace,
                                                        rres, ires, fltres, scratchpad, fftPlan, dctPlan );
    
    //================================================ Do the initial discrete split transform
    initialSplitDiscreteTransform                     ( oneDim, inputReal, inputImag, rres, ires, sphereMappedData, fftPlan, normCoeff );
    
    //================================================ Complete the spherical harmonics transform
    computeSphericalTransformCoeffs                   ( band, rdataptr, idataptr, outputReal, outputImag, rres, ires, fltres, scratchpad, tablePml, shWeights, dctPlan );
    
    //================================================ Apply the Condon-Shortley phase and save result to the final array
    applyCondonShortleyPhase                          ( band, outputReal, outputImag, shArray );
    
    //================================================ Free memory
    releaseSphericalMemory                            ( inputReal, inputImag, outputReal, outputImag, tableSpaceHelper, tablePml, shWeights, workspace, fftPlan, dctPlan );
    
    //================================================ Done
    return ;
    
}
