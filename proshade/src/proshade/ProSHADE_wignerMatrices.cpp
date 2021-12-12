/*! \file ProSHADE_wignerMatrices.cpp
    \brief This source file contains all the functions required to compute the Wigner D matrices
 
    The functions in this source file all serve to compute the Wigner D matrices, which are used to achieve two things, firstly to allow the inverse SOFT transform, which uses the Wigner D matrices as its
    basis functions and secondly, when a rotation of the internal map is required, these serve to compute the rotation in the spherical harmonics coefficient space.
 
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
#include "ProSHADE_wignerMatrices.hpp"

/*! \brief This function allocates the memory for the Wigner matrices for the calling object.
 */
void ProSHADE_internal_data::ProSHADE_data::allocateWignerMatricesSpace ( )
{
    //================================================ Sanity check
    if ( this->maxCompBand == 0 )
    {
        throw ProSHADE_exception ( "Attempted allocating Wigner D matrices before\n                    : allocating E matrices memory.", "EW00024", __FILE__, __LINE__, __func__, "The E matrices and Wigner matrices both require to know\n                    : the bandwidth of the comparison (which may differ from the\n                    : object bandwidth). This is set when allocating E matrices\n                    : and therefore if it is 0 now, E matrices were not  yet\n                    : allocated." );
    }
    
    //================================================ Allocate bands
    this->wignerMatrices                              = new proshade_complex** [this->maxCompBand];
    ProSHADE_internal_misc::checkMemoryAllocation ( this->wignerMatrices, __FILE__, __LINE__, __func__ );
    
    //================================================ Allocate the arrays
    for ( proshade_unsign bandIter = 0; bandIter < this->maxCompBand; bandIter++ )
    {
        //============================================ Allocate order 1
        this->wignerMatrices[bandIter]                = new proshade_complex* [(bandIter * 2) + 1];
        ProSHADE_internal_misc::checkMemoryAllocation ( this->wignerMatrices[bandIter], __FILE__, __LINE__, __func__ );
        
        //============================================ Allocate order 2
        for ( proshade_unsign order1Iter = 0; order1Iter < ( (bandIter * 2) + 1 ); order1Iter++ )
        {
            this->wignerMatrices[bandIter][order1Iter] = new proshade_complex [(bandIter * 2) + 1];
            ProSHADE_internal_misc::checkMemoryAllocation ( this->wignerMatrices[bandIter], __FILE__, __LINE__, __func__ );
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allocates the memory for the Wigner matrices computation.
 
    This function allocates all the internal computation space, which will be required only for the duration of the computation
    and not after.
 
    \param[in] matIn A pointer to the array holding the inputs for Wigner d matrix computation by SOFT.
    \param[in] matOut A pointer to the array holding the outputs for Wigner d matrix computation by SOFT.
    \param[in] sqrts Array of square roots of the bands.
    \param[in] workspace Array required by SOFT for internal computation space.
    \param[in] alphaExponentReal An array of exponents of the alpha angle for the real part of the Wigner matrices.
    \param[in] alphaExponentImag An array of exponents of the alpha angle for the imaginary part of the Wigner matrices.
    \param[in] gammaExponentReal An array of exponents of the gamma angle for the real part of the Wigner matrices.
    \param[in] gammaExponentImag An array of exponents of the gamma angle for the imaginary part of the Wigner matrices.
    \param[in] trigs Array of 2 doubles holding the  trigonometric results for the beta angle.
    \param[in] compBand The bendwidth of the computation.
 */
void ProSHADE_internal_wigner::allocateWignerWorkspace ( proshade_double*& matIn, proshade_double*& matOut, proshade_double*& sqrts, proshade_double*& workspace, proshade_double*& alphaExponentReal, proshade_double*& alphaExponentImag, proshade_double*& gammaExponentReal, proshade_double*& gammaExponentImag, proshade_double*& trigs, proshade_unsign compBand )
{
    //================================================ Allocate the memory
    matIn                                             = new proshade_double[static_cast<proshade_unsign> ( 4 * pow( compBand, 2.0 ) ) - 4 * compBand + 1];
            
    matOut                                            = new proshade_double[static_cast<proshade_unsign> ( 4 * pow( compBand, 2.0 ) ) - 4 * compBand + 1];
    sqrts                                             = new proshade_double[static_cast<proshade_unsign> ( 2 * compBand )];
    workspace                                         = new proshade_double[static_cast<proshade_unsign> ( 4 * pow( compBand, 2.0 ) )];
    alphaExponentReal                                 = new proshade_double[static_cast<proshade_unsign> ( 2 * compBand - 1 )];
    alphaExponentImag                                 = new proshade_double[static_cast<proshade_unsign> ( 2 * compBand - 1 )];
    gammaExponentReal                                 = new proshade_double[static_cast<proshade_unsign> ( 2 * compBand - 1 )];
    gammaExponentImag                                 = new proshade_double[static_cast<proshade_unsign> ( 2 * compBand - 1 )];
    trigs                                             = new proshade_double[2];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( matIn,             __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( matOut,            __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( sqrts,             __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( workspace,         __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( alphaExponentReal, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( alphaExponentImag, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( gammaExponentReal, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( gammaExponentImag, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( trigs,             __FILE__, __LINE__, __func__ );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function releases the memory for the Wigner matrices computation.
 
    \param[in] matIn A pointer to the array holding the inputs for Wigner d matrix computation by SOFT.
    \param[in] matOut A pointer to the array holding the outputs for Wigner d matrix computation by SOFT.
    \param[in] sqrts Array of square roots of the bands.
    \param[in] workspace Array required by SOFT for internal computation space.
    \param[in] alphaExponentReal An array of exponents of the alpha angle for the real part of the Wigner matrices.
    \param[in] alphaExponentImag An array of exponents of the alpha angle for the imaginary part of the Wigner matrices.
    \param[in] gammaExponentReal An array of exponents of the gamma angle for the real part of the Wigner matrices.
    \param[in] gammaExponentImag An array of exponents of the gamma angle for the imaginary part of the Wigner matrices.
 */
void ProSHADE_internal_wigner::releaseWignerWorkspace ( proshade_double*& matIn, proshade_double*& matOut, proshade_double*& sqrts, proshade_double*& workspace, proshade_double*& alphaExponentReal, proshade_double*& alphaExponentImag, proshade_double*& gammaExponentReal, proshade_double*& gammaExponentImag, proshade_double*& trigs )
{
    //================================================ Allocate the memory
    if ( matIn             != nullptr ) { delete[] matIn;             }
    if ( matOut            != nullptr ) { delete[] matOut;            }
    if ( sqrts             != nullptr ) { delete[] sqrts;             }
    if ( workspace         != nullptr ) { delete[] workspace;         }
    if ( trigs             != nullptr ) { delete[] trigs;             }
    if ( alphaExponentReal != nullptr ) { delete[] alphaExponentReal; }
    if ( alphaExponentImag != nullptr ) { delete[] alphaExponentImag; }
    if ( gammaExponentReal != nullptr ) { delete[] gammaExponentReal; }
    if ( gammaExponentImag != nullptr ) { delete[] gammaExponentImag; }

    //================================================ Done
    return ;
    
}

/*! \brief This function sets all the values repeatedly required for the computation.
 
    \param[in] sqrts Array of square roots of the bands.
    \param[in] alphaExponentReal An array of exponents of the alpha angle for the real part of the Wigner matrices.
    \param[in] alphaExponentImag An array of exponents of the alpha angle for the imaginary part of the Wigner matrices.
    \param[in] gammaExponentReal An array of exponents of the gamma angle for the real part of the Wigner matrices.
    \param[in] gammaExponentImag An array of exponents of the gamma angle for the imaginary part of the Wigner matrices.
    \param[in] trigs proshade_complex pointer which will hold the beta angle trigonometric function results.
    \param[in] compBand The bendwidth of the computation.
    \param[in] angAlpha The Euler alpha angle by which to rotate.
    \param[in] angBeta The Euler beta angle by which to rotate.
    \param[in] angGamma The Euler gamma angle by which to rotate.
 */
void ProSHADE_internal_wigner::prepareTrigsSqrtsAndExponents ( proshade_double* sqrts, proshade_double* alphaExponentReal, proshade_double* alphaExponentImag, proshade_double* gammaExponentReal, proshade_double* gammaExponentImag, proshade_double* trigs, proshade_unsign compBand, proshade_double angAlpha, proshade_double angBeta, proshade_double angGamma )
{
    //================================================ Compute the square roots
    for ( proshade_unsign iter = 0; iter < ( 2 * compBand ); iter++ )
    {
        sqrts[iter]                                   = static_cast<proshade_double> ( sqrt ( static_cast<proshade_double> ( iter ) ) );
    }
    
    //================================================ Compute the trig values
    trigs[0]                                          = static_cast<proshade_double> ( cos ( 0.5 * -angBeta ) );
    trigs[1]                                          = static_cast<proshade_double> ( sin ( 0.5 * -angBeta ) );
    
    //================================================ Get alpha and gamma exponents
    genExp                                            ( static_cast< int > ( compBand ), angAlpha, alphaExponentReal, alphaExponentImag );
    genExp                                            ( static_cast< int > ( compBand ), angGamma, gammaExponentReal, gammaExponentImag );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function does the actual computation of the Wigner D matrices.
 
    This function is the workhorse of the Wigner D matrices computation. It does iterate throught all the appropriate bands and order combinations
    and it computes the Wigner D matrix (l,m,m') values using the Wigner d matrices (which only take into account the beta Euler ZXZ angle) and the
    exponents obtained using the Euler ZXZ alpha and gamma angles. It also deals with the signs, but it does not deal with the memory allocation,
    release and general value set-up.
 
    \param[in] settings A pointer to settings class containing all the information required for the task.
    \param[in] obj A ProSHADE_data class object for which the Wigner matrices should be computed.
    \param[in] alphaExponentReal An array of exponents of the alpha angle for the real part of the Wigner matrices.
    \param[in] alphaExponentImag An array of exponents of the alpha angle for the imaginary part of the Wigner matrices.
    \param[in] gammaExponentReal An array of exponents of the gamma angle for the real part of the Wigner matrices.
    \param[in] gammaExponentImag An array of exponents of the gamma angle for the imaginary part of the Wigner matrices.
    \param[in] matIn A pointer to the array holding the inputs for Wigner d matrix computation by SOFT.
    \param[in] matOut A pointer to the array holding the outputs for Wigner d matrix computation by SOFT.
    \param[in] trigs proshade_complex pointer which holds the beta angle trigonometric function results.
    \param[in] sqrts Array of square roots of the bands.
    \param[in] workspace Array required by SOFT for internal computation space.
 */
void ProSHADE_internal_wigner::computeWignerMatrices ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* obj, proshade_double* alphaExponentReal, proshade_double* alphaExponentImag, proshade_double* gammaExponentReal, proshade_double* gammaExponentImag, proshade_double* matIn, proshade_double* matOut, proshade_double* trigs, proshade_double* sqrts, proshade_double* workspace )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 4, "Start Wigner D matrix computation.", settings->messageShift );
    
    //================================================ For each band, find the Wigned d matrix
    proshade_double *expARStart, *expAIStart, *expGRStart, *expGIStart;
    proshade_double Dij, eARi, eAIi, eGRj, eGIj, iSign, rSign;
    proshade_complex hlpVal;
    proshade_unsign noOrders, arrConvIter;
    for ( proshade_unsign bandIter = 0; bandIter < obj->getComparisonBand(); bandIter++ )
    {
        //============================================ Initialise loop
        noOrders                                      = 2 * bandIter + 1;
        arrConvIter                                   = 0;
        expARStart                                    = &alphaExponentReal[ (obj->getComparisonBand() - 1) - bandIter ];
        expAIStart                                    = &alphaExponentImag[ (obj->getComparisonBand() - 1) - bandIter ];
        expGRStart                                    = &gammaExponentReal[ (obj->getComparisonBand() - 1) - bandIter ];
        expGIStart                                    = &gammaExponentImag[ (obj->getComparisonBand() - 1) - bandIter ];
        iSign                                         = 1.0;
        rSign                                         = 1.0;
        
        //============================================ Get wigner d matrix values using beta angles only
        wignerdmat                                    ( static_cast< int > ( bandIter ), matIn, matOut, trigs, sqrts, workspace );
        
        //============================================ Multiply the wigner d matrix by alpha and gamma values and save the wigner D matrix to output array
        for ( proshade_unsign d1Iter = 0; d1Iter < noOrders; d1Iter++ )
        {
            eARi                                      = expARStart[d1Iter];
            eAIi                                      = expAIStart[d1Iter];
            
            for ( proshade_unsign d2Iter = 0; d2Iter < noOrders; d2Iter++ )
            {
                Dij                                   = matOut[arrConvIter];
                eGRj                                  = expGRStart[d2Iter];
                eGIj                                  = expGIStart[d2Iter];
                        
                hlpVal[0]                             = ( Dij * eGRj * eARi - Dij * eGIj * eAIi ) * rSign;
                hlpVal[1]                             = ( Dij * eGRj * eAIi + Dij * eGIj * eARi ) * iSign;
                obj->setWignerMatrixValue             ( hlpVal, bandIter, d1Iter, d2Iter );
                        
                arrConvIter                          += 1;
                iSign                                *= -1.0;
                rSign                                *= -1.0;
            }
        }
        
        //============================================ Get ready for next wigner matrix calculation
        memcpy                                        ( matIn, matOut, sizeof ( proshade_double ) * ( noOrders * noOrders ) );
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 5, "Wigner D matrices obtained.", settings->messageShift );

    //================================================ Done
    return ;
    
}

/*! \brief This function computes the Wigner D matrices for a particular set of Euler angles.
 
    This function starts by allocating the required memory to store the Wigner D matrices for a particular object and then it
    proceeds to allocate all the computation memory as well. Subsequently, it prepares all the values needed for the computation
    and it calls the workhhorse computing function. Once complete, it deletes all the now redundant memory and exits.
 
    \param[in] settings A pointer to settings class containing all the information required for the task.
    \param[in] obj A ProSHADE_data class object for which the Wigner matrices should be computed.
    \param[in] eulerAlpha The Euler ZXZ convention alpha angle value for the rotation in SO(3) space.
    \param[in] eulerBeta The Euler ZXZ convention beta angle value for the rotation in SO(3) space.
    \param[in] eulerGamma The Euler ZXZ convention gamma angle value for the rotation in SO(3) space.
 */
void ProSHADE_internal_wigner::computeWignerMatricesForRotation ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* obj, proshade_double eulerAlpha, proshade_double eulerBeta, proshade_double eulerGamma )
{
    //================================================ Initialise local variables
    proshade_double *matIn, *matOut, *sqrts, *workspace, *alphaExponentReal, *alphaExponentImag, *gammaExponentReal, *gammaExponentImag, *trigs;
    
    //================================================ Allocate memory for Wigner matrices
    obj->allocateWignerMatricesSpace                  ( );
    
    //================================================ Allocate the workspace memory
    allocateWignerWorkspace                           ( matIn, matOut, sqrts, workspace, alphaExponentReal, alphaExponentImag,
                                                        gammaExponentReal, gammaExponentImag, trigs, obj->getComparisonBand() );
    
    //================================================ Prepare all values for the computation
    prepareTrigsSqrtsAndExponents                     ( sqrts, alphaExponentReal, alphaExponentImag, gammaExponentReal, gammaExponentImag,
                                                        trigs, obj->getComparisonBand(), eulerAlpha, eulerBeta, eulerGamma );
    
    //================================================ Compute the values
    computeWignerMatrices                             ( settings, obj, alphaExponentReal, alphaExponentImag, gammaExponentReal, gammaExponentImag,
                                                        matIn, matOut, trigs, sqrts, workspace );
    
    //================================================ Release the workspace memory
    releaseWignerWorkspace                            ( matIn, matOut, sqrts, workspace, alphaExponentReal, alphaExponentImag,
                                                        gammaExponentReal, gammaExponentImag, trigs );

    //================================================ Done
    return ;
    
}
