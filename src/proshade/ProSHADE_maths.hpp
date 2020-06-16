/*! \file ProSHADE_maths.hpp
    \brief This header file declares all the functions required for computing various information from the ProSHADE data.
 
    This header file declares the ProSHADE_internal_maths namespace, which groups all the functions required to computed various information from
    the specific ProSHADE data and its organisation. The functionalities available here include complex number computations, rotation representation
    conversions as well as Gauss-Legendre integration or Taylor series approximation.
 
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
#include "ProSHADE_misc.hpp"

//==================================================== Overinclusion protection
#ifndef __PROSHADE_MATHS__
#define __PROSHADE_MATHS__

//==================================================== Declare LAPACK functions
extern "C"
{
    // ... The complex matrix singular value decomposition function
    extern void zgesdd_( char* jobz, int* m, int* n, std::complex<double>* a, int* lda, double* s, std::complex<double>* u, int* ldu, std::complex<double>* vt, int* ldvt, std::complex<double>* work, int* lwork, double* rwork, int* iwork, int* info );
}

//==================================================== ProSHADE_internal_io Namespace
/*! \namespace ProSHADE_internal_maths
 \brief This namespace contains the internal functions for common mathematical operations.
 
 The ProSHADE_internal_maths namespace contains a set of common mathematical operations used in many places by ProSHADE. These
 typically include complex number operations and angle conversions.
 */
namespace ProSHADE_internal_maths
{
    void complexMultiplication                        ( proshade_double* r1, proshade_double* i1, proshade_double* r2, proshade_double* i2,
                                                        proshade_double* retReal, proshade_double* retImag );
    void complexMultiplicationConjug                  ( proshade_double* r1, proshade_double* i1, proshade_double* r2, proshade_double* i2,
                                                        proshade_double* retReal, proshade_double* retImag );
    proshade_double complexMultiplicationRealOnly     ( proshade_double* r1, proshade_double* i1, proshade_double* r2, proshade_double* i2 );
    proshade_double complexMultiplicationConjugRealOnly ( proshade_double* r1, proshade_double* i1, proshade_double* r2, proshade_double* i2 );
    void vectorMeanAndSD                              ( std::vector<proshade_double>* vec, proshade_double*& ret );
    void vectorMedianAndIQR                           ( std::vector<proshade_double>* vec, proshade_double*& ret );
    void arrayMedianAndIQR                            ( proshade_double* vec, proshade_unsign vecSize, proshade_double*& ret );
    proshade_double pearsonCorrCoeff                  ( proshade_double* valSet1, proshade_double* valSet2, proshade_unsign length );
    void getLegendreAbscAndWeights                    ( proshade_unsign order, proshade_double* abscissas, proshade_double* weights,
                                                        proshade_unsign taylorSeriesCap );
    void getGLPolyAtZero                              ( proshade_unsign order, proshade_double *polyValue, proshade_double *deriValue );
    void getGLFirstEvenRoot                           ( proshade_double polyAtZero, proshade_unsign order, proshade_double *abscAtZero,
                                                        proshade_double *weighAtZero, proshade_unsign taylorSeriesCap );
    proshade_double evaluateGLSeries                  ( proshade_double *series, proshade_double target, proshade_unsign terms );
    proshade_double advanceGLPolyValue                ( proshade_double from, proshade_double to, proshade_double valAtFrom,
                                                        proshade_unsign noSteps, proshade_unsign taylorSeriesCap );
    void completeLegendreSeries                       ( proshade_unsign order, proshade_double* abscissa, proshade_double* weights,
                                                        proshade_unsign taylorSeriesCap );
    proshade_double gaussLegendreIntegrationReal      ( proshade_double* vals, proshade_unsign valsSize, proshade_unsign order,
                                                        proshade_double* abscissas, proshade_double* weights, proshade_double integralOverRange,
                                                        proshade_double maxSphereDists );
    void gaussLegendreIntegration                     ( proshade_complex* vals, proshade_unsign valsSize, proshade_unsign order,
                                                        proshade_double* abscissas, proshade_double* weights, proshade_double integralOverRange,
                                                        proshade_double maxSphereDists, proshade_double* retReal, proshade_double* retImag );
    void complexMatrixSVDSigmasOnly                   ( proshade_complex** mat, int dim, double*& singularValues );
    void complexMatrixSVDUandVOnly                    ( proshade_double* mat, int dim, proshade_double* uAndV, bool fail = true );
    void getEulerZXZFromSOFTPosition                  ( proshade_signed band, proshade_signed x, proshade_signed y, proshade_signed z, proshade_double* eulerAlpha,
                                                        proshade_double* eulerBeta, proshade_double* eulerGamma );
    void getRotationMatrixFromEulerZXZAngles          ( proshade_double eulerAlpha, proshade_double eulerBeta, proshade_double eulerGamma, proshade_double* matrix );
    void getAxisAngleFromRotationMatrix               ( proshade_double* rotMat, proshade_double* x, proshade_double* y, proshade_double* z, proshade_double* ang );
    void getRotationMatrixFromAngleAxis               ( proshade_double* rotMat, proshade_double x, proshade_double y, proshade_double z, proshade_double ang );
    void getEulerZXZFromRotMatrix                     ( proshade_double* rotMat, proshade_double* eA, proshade_double* eB, proshade_double* eG );
    void multiplyTwoSquareMatrices                    ( proshade_double* A, proshade_double* B, proshade_double* res, proshade_unsign dim );
    std::vector < proshade_signed > primeFactorsDecomp ( proshade_signed number );
    proshade_double normalDistributionValue           ( proshade_double mean, proshade_double standardDev, proshade_double value );
    proshade_double computeDotProduct                 ( proshade_double* x1, proshade_double* y1, proshade_double* z1, proshade_double* x2, proshade_double* y2,
                                                        proshade_double* z2 );
}

#endif
