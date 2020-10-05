/*! \file ProSHADE_maths.cpp
    \brief This source file contains all the mathematical functions not simply available from elsewhere or modified to work with ProSHADE specific data formats.
 
    The functions in this source file provide the computational power for the rest of the code. The functions here are specifically written to work with the ProSHADE internal data organisation and
    are used throughtout the rest of the code. The functionalities implemented here range from complex number maths, Taylor series and Gauss-Legendre integration to conversion betwee different
    rotation representations.
 
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
#include "ProSHADE_maths.hpp"

/*! \brief Function to multiply two complex numbers.
 
    This function takes pointers to the real and imaginary parts of two complex numbers and
    returns the result of their multiplication.
 
    \param[in] r1 Pointer to the real value of number 1.
    \param[in] i1 Pointer to the imaginary value of number 1.
    \param[in] r2 Pointer to the real value of number 2.
    \param[in] i2 Pointer to the imaginary value of number 2.
    \param[in] retReal Pointer to the real part of the complex variable to which the result will be saved.
    \param[in] retImag Pointer to the imaginary part of the complex variable to which the result will be saved.
 */
void ProSHADE_internal_maths::complexMultiplication ( proshade_double* r1, proshade_double* i1, proshade_double* r2, proshade_double* i2, proshade_double* retReal, proshade_double* retImag )
{
    //================================================ Multiplication
    *retReal                                          = (*r1)*(*r2) - (*i1)*(*i2);
    *retImag                                          = (*r1)*(*i2) + (*i1)*(*r2);
    
    //================================================ Return
    return ;
    
}

/*! \brief Function to multiply two complex numbers by using the second number's conjugate.
 
    This function takes pointers to the real and imaginary parts of two complex numbers and
    returns the result of their multiplication, while using the conjugate of the second complex
    number.
 
    \param[in] r1 Pointer to the real value of number 1.
    \param[in] i1 Pointer to the imaginary value of number 1.
    \param[in] r2 Pointer to the real value of number 2.
    \param[in] i2 Pointer to the imaginary value of number 2.
    \param[in] retReal Pointer to the real part of the complex variable to which the result will be saved.
    \param[in] retImag Pointer to the imaginary part of the complex variable to which the result will be saved.
 */
void ProSHADE_internal_maths::complexMultiplicationConjug ( proshade_double* r1, proshade_double* i1, proshade_double* r2, proshade_double* i2, proshade_double* retReal, proshade_double* retImag )
{
    //================================================ Multiplication
    *retReal                                          =  (*r1)*(*r2) + (*i1)*(*i2);
    *retImag                                          = -(*r1)*(*i2) + (*i1)*(*r2);
    
    //================================================ Return
    return ;
    
}

/*! \brief Function to multiply two complex numbers and return the real part only.
 
    This function takes pointers to the real and imaginary parts of two complex numbers and
    returns the real part of the result of their multiplication.
 
    \param[in] r1 Pointer to the real value of number 1.
    \param[in] i1 Pointer to the imaginary value of number 1.
    \param[in] r2 Pointer to the real value of number 2.
    \param[in] i2 Pointer to the imaginary value of number 2.
 */
proshade_double ProSHADE_internal_maths::complexMultiplicationRealOnly ( proshade_double* r1, proshade_double* i1, proshade_double* r2, proshade_double* i2 )
{
    //================================================ Multiplication
    proshade_double ret                               = (*r1)*(*r2) - (*i1)*(*i2);
    
    //================================================ Return
    return                                            ( ret );
    
}

/*! \brief Function to conjuggate multiply two complex numbers and return the real part only.
 
    This function takes pointers to the real and imaginary parts of two complex numbers and
    returns the real part of the result of their conjugate multiplication.
 
    \param[in] r1 Pointer to the real value of number 1.
    \param[in] i1 Pointer to the imaginary value of number 1.
    \param[in] r2 Pointer to the real value of number 2.
    \param[in] i2 Pointer to the imaginary value of number 2.
 */
proshade_double ProSHADE_internal_maths::complexMultiplicationConjugRealOnly ( proshade_double* r1, proshade_double* i1, proshade_double* r2, proshade_double* i2 )
{
    //================================================ Multiplication
    proshade_double ret                               = (*r1)*(*r2) + (*i1)*(*i2);
    
    //================================================ Return
    return                                            ( ret );
    
}

/*! \brief Function to get vector mean and standard deviation.
 
    This function takes a pointer to a vector of proshade_double's and returns the mean and standard deviation of
    such vector.
 
    \param[in] vec Pointer to a vector of proshade_double's for which mean and sd should be obtained.
    \param[in] ret Pointer to array of 2 proshade_double's, which will be the return values - first mean and second sd.
 */
void ProSHADE_internal_maths::vectorMeanAndSD ( std::vector<proshade_double>* vec, proshade_double*& ret )
{
    //================================================ Get mean
    ret[0]                                            = std::accumulate ( vec->begin(), vec->end(), 0.0 ) / static_cast<proshade_double> ( vec->size() );
    
    //================================================ Get standard deviation
    proshade_double squaredSum                        = std::inner_product ( vec->begin(), vec->end(), vec->begin(), 0.0 );
    ret[1]                                            = std::sqrt ( ( squaredSum / static_cast<proshade_double> ( vec->size() ) ) - std::pow ( ret[0], 2.0 ) );
    
    //================================================ Return
    return ;
    
}

/*! \brief Function to get vector median and inter-quartile range.
 
    This function takes a pointer to a vector of proshade_double's and returns the median and the inter-quartile range of
    such vector.
 
    \param[in] vec Pointer to a vector of proshade_double's for which median and IQR should be obtained.
    \param[in] ret Pointer to array of 2 proshade_double's, which will be the return values - first median and second IQR.
 */
void ProSHADE_internal_maths::vectorMedianAndIQR ( std::vector<proshade_double>* vec, proshade_double*& ret )
{
    //================================================ Sanity check
    if ( vec->size() < 3 ) { ret[0] = 0.0; ret[1] = 0.0; return; }
    
    //================================================ Sort the vector
    std::sort                                         ( vec->begin(), vec->end() );
    
    //================================================ Get median
    if ( static_cast<proshade_unsign> ( vec->size() ) % 2 == 0)
    {
        ret[0]                                        = ( vec->at( ( static_cast<proshade_unsign> ( vec->size() ) / 2 ) - 1 ) +
                                                        vec->at(   static_cast<proshade_unsign> ( vec->size() ) / 2 ) ) / 2.0;
    }
    else
    {
        ret[0]                                        = vec->at( static_cast<proshade_unsign> ( vec->size() ) / 2 );
    }
    
    //================================================ Get first and third quartile
    proshade_double Q1, Q3;
    if ( static_cast<proshade_unsign> ( vec->size() ) % 2 == 0)
    {
        Q1                                            = ( vec->at( ( static_cast<proshade_unsign> ( vec->size() ) / 4 ) - 1 ) +
                                                          vec->at(   static_cast<proshade_unsign> ( vec->size() ) / 4 ) ) / 2.0;
        Q3                                            = ( vec->at( ( ( static_cast<proshade_unsign> ( vec->size() ) / 4 ) * 3 ) - 1 ) +
                                                          vec->at(   ( static_cast<proshade_unsign> ( vec->size() ) / 4 ) * 3 ) ) / 2.0;
    }
    else
    {
        Q1                                            = vec->at( static_cast<proshade_unsign> ( vec->size() ) / 4 );
        Q3                                            = vec->at( ( static_cast<proshade_unsign> ( vec->size() ) / 4 ) * 3 );
    }
    
    //================================================ And now save the IQR
    ret[1]                                            = Q3 -  Q1;
    
    //================================================ Return
    return ;
    
}

/*! \brief Function to get array median and inter-quartile range.
 
    This function takes a pointer to a array of proshade_double's and returns the median and the inter-quartile range of
    such vector.
 
    \param[in] vec Pointer to an array of proshade_double's for which median and IQR should be obtained.
    \param[in] vecSize The length of the array.
    \param[in] ret Pointer to array of 2 proshade_double's, which will be the return values - first median and second IQR.
 */
void ProSHADE_internal_maths::arrayMedianAndIQR ( proshade_double* vec, proshade_unsign vecSize, proshade_double*& ret )
{
    //================================================ Sort the vector
    std::sort                                         ( vec, vec + vecSize );
    
    //================================================ Get median
    if ( vecSize % 2 == 0)
    {
        ret[0]                                        = ( vec[ ( vecSize / 2 ) - 1 ] + vec[ vecSize / 2 ] ) / 2.0;
    }
    else
    {
        ret[0]                                        = vec[ vecSize / 2 ];
    }
    
    //================================================ Get first and third quartile
    proshade_double Q1, Q3;
    if ( vecSize % 2 == 0)
    {
        Q1                                            = ( vec[ ( vecSize / 4 ) - 1 ] + vec[ vecSize / 4 ] ) / 2.0;
        Q3                                            = ( vec[ ( ( vecSize / 4 ) * 3 ) - 1 ] + vec[ ( vecSize / 4 ) * 3 ] ) / 2.0;
    }
    else
    {
        Q1                                            = vec[ vecSize / 4 ];
        Q3                                            = vec[ ( vecSize / 4 ) * 3 ];
    }
    
    //================================================ And now save the IQR
    ret[1]                                            = Q3 -  Q1;
    
    //================================================ Return
    return ;
    
}

/*! \brief Function for computing the Pearson's correlation coefficient.
 
    This function takes two numerical arrays of same length and proceeds to compute the Pearson's
    correlation coefficient, which it then returns.
 
    \param[in] valSet1 This is the set of x-values.
    \param[in] valSet2 This is the set of y-values.
    \param[in] length The length of both arrays (both arrays have to have the same length).
    \param[out] X The Pearson's correlation coefficient value.
 */
proshade_double ProSHADE_internal_maths::pearsonCorrCoeff ( proshade_double* valSet1, proshade_double* valSet2, proshade_unsign length )
{
    //================================================ Find vector means
    proshade_double xMean                             = 0.0;
    proshade_double yMean                             = 0.0;
    proshade_double zeroCount                         = 0.0;
    for ( proshade_unsign iter = 0; iter < length; iter++ )
    {
        xMean                                        += valSet1[iter];
        yMean                                        += valSet2[iter];
    }
    xMean                                            /= static_cast<proshade_double> ( length - zeroCount );
    yMean                                            /= static_cast<proshade_double> ( length - zeroCount );
    
    //================================================ Get Pearson's correlation coefficient
    proshade_double xmmymm                            = 0.0;
    proshade_double xmmsq                             = 0.0;
    proshade_double ymmsq                             = 0.0;
    for ( proshade_unsign iter = 0; iter < length; iter++ )
    {
        xmmymm                                       += ( valSet1[iter] - xMean ) * ( valSet2[iter] - yMean );
        xmmsq                                        += pow( valSet1[iter] - xMean, 2.0 );
        ymmsq                                        += pow( valSet2[iter] - yMean, 2.0 );
    }
    
    proshade_double ret                               = xmmymm / ( sqrt(xmmsq) * sqrt(ymmsq) );
    
    //================================================ Done
    if ( std::isnan ( ret ) ) { return ( 0.0 ); }
    return                                            ( ret );
    
}

/*! \brief Function to prepare abscissas and weights for Gauss-Legendre integration.
 
    This function fills in the Gauss-Legendre interpolation points positions (abscissas) and their weights vectors, which will then be used for computing the
    Gauss-Legendre interpolation.
 
    \param[in] order The order to which the abscissas and weights should be prepared.
    \param[in] abscissas The array holding the abscissa values.
    \param[in] weights The array holding the weight values.
    \param[in] taylorSeriesCap The limit on the Taylor series.
 */
void ProSHADE_internal_maths::getLegendreAbscAndWeights ( proshade_unsign order, proshade_double* abscissas, proshade_double* weights, proshade_unsign taylorSeriesCap )
{
    //================================================ Sanity check
    if ( order < 2 )
    {
        throw ProSHADE_exception ( "The integration order is too low.", "EI00019", __FILE__, __LINE__, __func__, "The Gauss-Legendre integration order is less than 2. This\n                    : seems very low; if you have a very small structure or very\n                    : low resolution, please manually increase the integration\n                    : order. Otherwise, please report this as a bug." );
    }
    
    //================================================ Initialise
    proshade_double polyValue                         = 0.0;
    proshade_double deriValue                         = 0.0;
    proshade_double weightSum                         = 0.0;
    
    //================================================ Find the polynomial and derivative values at 0
    getGLPolyAtZero                                   ( order,
                                                       &polyValue,
                                                       &deriValue );
    
    //================================================ If the order is odd, then 0 is a root ...
    if ( order % 2 == 1 )
    {
        abscissas[((order-1)/2)]                      = polyValue;
        weights[((order-1)/2)]                        = deriValue;
    }
    else
    {
        // ... and if order is even, find the first root
        getGLFirstEvenRoot                            ( polyValue, order, &abscissas[(order/2)], &weights[(order/2)], taylorSeriesCap );
    }

    //================================================ Now, having computed the first roots, complete the series
    completeLegendreSeries                            ( order, abscissas, weights, taylorSeriesCap );

    //================================================ Correct weights by anscissa values
    for ( proshade_unsign iter = 0; iter < order; iter++ )
    {
        weights[iter]                                 = 2.0 / ( 1.0 - abscissas[iter] ) / ( 1.0 + abscissas[iter] ) / weights[iter] / weights[iter];
        weightSum                                     = weightSum + weights[iter];
    }

    //================================================ Normalise weights
    for ( proshade_unsign iter = 0; iter < order; iter++ )
    {
        weights[iter]                                 = 2.0 * weights[iter] / weightSum;
    }
    
    //================================================ Done
    return ;
}

/*! \brief This function obtains the Legendre polynomial values and its derivative at zero for any positive integer order polynomial.
 
    This function takes the positive integer order of the Legendre polynomial and uses the recursive
    properties of the polynomials to work up to the order, computing the value at zero and its derivative
    for all lesser orders. It then returns the final values.
 
    \param[in] order Positive integer order of the Legendre polynomial which value at zero we want.
    \param[in] polyValue Pointer to variable which will store the resulting polynomial value at zero.
    \param[in] deriValue Pointer to variable which will store the derivative of the zero value.
 */
void ProSHADE_internal_maths::getGLPolyAtZero ( proshade_unsign order, proshade_double *polyValue, proshade_double *deriValue )
{
    //================================================ Initialise
    proshade_double hlpVal                            = 0.0;
    proshade_double prevPoly                          = 1.0;
    proshade_double prevPrevPoly                      = 0.0;
    proshade_double prevDeri                          = 0.0;
    proshade_double prevPrevDeri                      = 0.0;
    
    for ( proshade_unsign ordIt = 0; ordIt < order; ordIt++ )
    {
        hlpVal                                        = static_cast<proshade_double> ( ordIt );
        *polyValue                                    = -hlpVal * prevPrevPoly / ( hlpVal + 1.0 );
        *deriValue                                    = ( ( 2.0 * hlpVal + 1.0 ) * prevPoly - hlpVal * prevPrevDeri ) / ( hlpVal + 1.0 );
        prevPrevPoly                                  = prevPoly;
        prevPoly                                      = *polyValue;
        prevPrevDeri                                  = prevDeri;
        prevDeri                                      = *deriValue;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds the first root for Legendre polynomials of odd order.
 
    The Legendre polynomials with odd order have zero as the first root, but the even oder polenomials
    have different value and this function serves the purpose of finding this value (i.e. the first
    root of the polynomial if the order is even).
 
    \param[in] polyAtZero The value of the polynomial at zero.
    \param[in] order The positive integer value of the polynomial order.
    \param[in] abscAtZero Pointer to variable storing the abscissa value at zero.
    \param[in] weightAtZero Pointer to variable storing the weight value at zero.
    \param[in] taylorSeriesCap The limit on the Taylor series.
 */
void ProSHADE_internal_maths::getGLFirstEvenRoot ( proshade_double polyAtZero, proshade_unsign order, proshade_double *abscAtZero, proshade_double *weighAtZero, proshade_unsign taylorSeriesCap )
{
    //================================================ Sanity check
    if ( taylorSeriesCap < 2 )
    {
        throw ProSHADE_exception ( "The Taylor series cap is too low.", "EI00020", __FILE__, __LINE__, __func__, "The Taylor series expansion limit is less than 2. This\n                    : seems very low; if you have a very small structure or very\n                    : low resolution, please manually increase the integration\n                    : order. Otherwise, please report this as a bug." );
    }
    
    //================================================ Initialise variables
   *abscAtZero                                        = advanceGLPolyValue ( 0.0, -M_PI / 2.0, 0.0, order, taylorSeriesCap );
    proshade_double hlp                               = 0.0;
    proshade_double hlpVal                            = static_cast<proshade_double> ( order );
    proshade_double *abscSteps;
    proshade_double *weightSteps;

    //================================================ Allocate memory
    abscSteps                                         = new proshade_double [taylorSeriesCap+2];
    weightSteps                                       = new proshade_double [taylorSeriesCap+1];

    //================================================ Pre-set values
    abscSteps[0]                                      = 0.0;
    abscSteps[1]                                      = polyAtZero;
    weightSteps[0]                                    = 0.0;

    //================================================ Fill in abscissa and weight steps
    for ( proshade_unsign iter = 0; iter <= taylorSeriesCap - 2; iter = iter + 2 )
    {
        hlp                                           = static_cast<proshade_double> ( iter );
        
        abscSteps[iter+2]                             = 0.0;
        abscSteps[iter+3]                             = ( hlp * ( hlp + 1.0 ) - hlpVal * ( hlpVal + 1.0 ) ) * abscSteps[iter+1] / (hlp + 1.0) / (hlp + 2.0 );
        
        weightSteps[iter+1]                           = 0.0;
        weightSteps[iter+2]                           = ( hlp + 2.0 ) * abscSteps[iter+3];
    }

    //================================================ Find abscissa and weights
    for ( proshade_double iter = 0; iter < 5; iter++ )
    {
        *abscAtZero                                   = *abscAtZero - evaluateGLSeries ( abscSteps, *abscAtZero, taylorSeriesCap ) / evaluateGLSeries ( weightSteps, *abscAtZero, taylorSeriesCap-1 );
    }
    *weighAtZero                                      = evaluateGLSeries ( weightSteps, *abscAtZero, taylorSeriesCap-1 );

    //================================================ Free memory
    delete abscSteps;
    delete weightSteps;

    //================================================ Done
    return ;
    
}

/*! \brief This function evaluates the Taylor expansion.
 
    This function takes the series array, the target value and the cap on Taylor expansion and proceeds to
    evaluate the series. The main use of this is to evaluate the series twice, one where the series evaluation
    'overshoots' and once where it 'undershoots' and taking value in between those, thus adding accuracy.
 
    \param[in] series Pointer to array with the series values.
    \param[in] target The target location on the series value.
    \param[in] terms The Taylor expansion cap.
    \param[out] X The value of the series at the target location.
 */
proshade_double ProSHADE_internal_maths::evaluateGLSeries ( proshade_double *series, proshade_double target, proshade_unsign terms )
{
    //================================================ Initalise
    proshade_double factorialValue                    = 1.0;
    proshade_double value                             = 0.0;
    
    //================================================ Compute
    for ( proshade_unsign iter = 1; iter <= terms; iter++ )
    {
        value                                         = value + series[iter] * factorialValue;
        factorialValue                                = factorialValue * target;
    }
    
    //================================================ Done
    return ( value );
    
}

/*! \brief This function finds the next value of the polynomial.
 
    Given the previous value of the polynomial, the distance to proceed and the number of steps to
    take, this function finds the next value of the polynomial using the Taylor series.
 
    \param[in] from Current polynomial position.
    \param[in] to Polynomial position to move to.
    \param[in] valAtFrom The current value of the polynomial at the <from> position.
    \param[in] noSteps Number of steps in which to reach the <to> position.
    \param[in] taylorSeriesCap The limit on the Taylor series.
    \param[out] X The polynomial value at the <to> position.
 */
proshade_double ProSHADE_internal_maths::advanceGLPolyValue ( proshade_double from, proshade_double to, proshade_double valAtFrom, proshade_unsign noSteps, proshade_unsign taylorSeriesCap )
{
    //================================================ Initialise variables
    proshade_double hlpVal                            = 0.0;
    proshade_double stepSize                          = 0.0;
    proshade_double valChange                         = 0.0;
    proshade_double valSecChange                      = 0.0;
    proshade_double squareSteps                       = 0.0;
    proshade_double curVal                            = 0.0;
    
    //================================================ Set initial values
    stepSize                                          = ( to - from ) / static_cast<proshade_double> ( taylorSeriesCap );
    squareSteps                                       = sqrt ( static_cast<proshade_double> ( noSteps * ( noSteps + 1 ) ) );
    curVal                                            = from;
    
    //================================================ Go through the series and iteratively improve the estimate
    for ( proshade_unsign iter = 0; iter < taylorSeriesCap; iter++ )
    {
        hlpVal                                        = ( 1.0 - valAtFrom ) * ( 1.0 + valAtFrom );
        valChange                                     = - stepSize * hlpVal / ( squareSteps * sqrt ( hlpVal ) - 0.5 * valAtFrom * sin ( 2.0 * curVal ) );
        valAtFrom                                     = valAtFrom + valChange;
                
        curVal                                        = curVal + stepSize;
                
        hlpVal                                        = ( 1.0 - valAtFrom ) * ( 1.0 + valAtFrom );
        valSecChange                                  = - stepSize * hlpVal / ( squareSteps * sqrt ( hlpVal ) - 0.5 * valAtFrom * sin ( 2.0 * curVal ) );
        valAtFrom                                     = valAtFrom + 0.5 * ( valSecChange - valChange );
    }
    
    //================================================ Done
    return valAtFrom;
    
}

/*! \brief This function completes the Legendre polynomial series assuming you have obtained the first values.
 
    Given that the polynomial value at zero is known, this function will complete the Legendre polynomial and with it
    the absicassas and weights for the Gauss-Legendre integration using the other functions defined above.
 
    \param[in] order The positive integer value of the polynomial order.
    \param[in] abscissas Pointer to an array of abscissas containing the first value.
    \param[in] weights Pointer to an array of weights containing the first value.
    \param[in] taylorSeriesCap The limit on the Taylor series.
 */
void ProSHADE_internal_maths::completeLegendreSeries ( proshade_unsign order, proshade_double* abscissas, proshade_double* weights, proshade_unsign taylorSeriesCap )
{
    //================================================ Initialise internal variables
    proshade_double hlpTaylorVal                      = 0.0;
    proshade_double hlpOrderVal                       = static_cast<proshade_double> ( order );
    proshade_double abscValueChange                   = 0.0;
    proshade_double prevAbsc                          = 0.0;
    proshade_double *hlpAbscSeries;
    proshade_double *hlpWeightSeries;
    proshade_unsign noSeriesElems                     = 0;
    proshade_unsign oddEvenSwitch                     = 0;
    
    //================================================ Pre-set internal values
    if ( order % 2 == 1 )
    {
        noSeriesElems                                 = ( order - 1 ) / 2 - 1;
        oddEvenSwitch                                 = 1;
    }
    else
    {
        noSeriesElems                                 = order / 2 - 1;
        oddEvenSwitch                                 = 0;
    }
    
    //================================================ Allocate memory
    hlpAbscSeries                                     = new proshade_double[taylorSeriesCap+2];
    hlpWeightSeries                                   = new proshade_double[taylorSeriesCap+1];
    
    //================================================ For each series element
    for ( proshade_unsign serIt = noSeriesElems + 1; serIt < order - 1; serIt++ )
    {
        //============================================ Init loop
        prevAbsc                                      = abscissas[serIt];
        abscValueChange                               = advanceGLPolyValue ( M_PI/2.0, -M_PI/2.0, prevAbsc, order, taylorSeriesCap ) - prevAbsc;
        
        //============================================ Init abscissas
        hlpAbscSeries[0]                              = 0.0;
        hlpAbscSeries[1]                              = 0.0;
        hlpAbscSeries[2]                              = weights[serIt];
        
        //============================================ Init weights
        hlpWeightSeries[0]                            = 0.0;
        hlpWeightSeries[1]                            = hlpAbscSeries[2];
        
        //============================================ Taylor expansion
        for ( proshade_unsign tayIt = 0; tayIt <= taylorSeriesCap - 2; tayIt++ )
        {
            hlpTaylorVal                              = static_cast<proshade_double> ( tayIt );
                    
            hlpAbscSeries[tayIt+3]                    = ( 2.0 * prevAbsc * ( hlpTaylorVal + 1.0 ) * hlpAbscSeries[tayIt+2] + ( hlpTaylorVal * ( hlpTaylorVal + 1.0 ) - hlpOrderVal *
                                                        ( hlpOrderVal + 1.0 ) ) * hlpAbscSeries[tayIt+1] / ( hlpTaylorVal + 1.0 ) ) / ( 1.0 - prevAbsc ) / ( 1.0 + prevAbsc ) /
                                                        ( hlpTaylorVal + 2.0 );
            
            hlpWeightSeries[tayIt+2]                  = ( hlpTaylorVal + 2.0 ) * hlpAbscSeries[tayIt+3];
        }
        
        //============================================ Sum over results
        for ( proshade_unsign iter = 0; iter < 5; iter++ )
        {
            abscValueChange                           = abscValueChange - evaluateGLSeries ( hlpAbscSeries,   abscValueChange, taylorSeriesCap   ) /
                                                                          evaluateGLSeries ( hlpWeightSeries, abscValueChange, taylorSeriesCap-1 );
        }
        
        //============================================ Save results
        abscissas[serIt+1]                            = prevAbsc + abscValueChange;
        weights[serIt+1]                              = evaluateGLSeries ( hlpWeightSeries, abscValueChange, taylorSeriesCap - 1 );
    }
    
    for ( proshade_unsign serIt = 0; serIt <= noSeriesElems + oddEvenSwitch; serIt++ )
    {
        abscissas[serIt]                              = -abscissas[order-serIt-1];
        weights[serIt]                                = weights[order-serIt-1];
    }
    
    //================================================ Free memory
    delete hlpAbscSeries;
    delete hlpWeightSeries;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function to compute real part of the Gauss-Legendre integration over spherical harmonic values in different shells.
 
    This function takes the real parts of the spherical harmonics value in different shells and proceeds to compute
    the Gauss-Legendre integration over them. It uses the shell positions to appropriately place abscissas and their
    weights, which it assumes were pre-computed by the getLegendreAbscAndWeights() function.
 
    \param[in] vals Pointer to an array of values over which the integration to be done.
    \param[in] valsSize The length of the input array.
    \param[in] order The integration order value.
    \param[in] abscissas The allocated array for holding the abscissa values.
    \param[in] weights The allocated array for holding the weight values.
    \param[in] integralOverRange The range of the intgral. If progressive shell mapping is used, this will not be max shell radius.
    \param[in] maxSphereDists Distance between two shells.
    \param[out] X The real part of Gauss-Legendre integration over the shperical harmonics values.
 */
proshade_double ProSHADE_internal_maths::gaussLegendreIntegrationReal ( proshade_double* vals, proshade_unsign valsSize, proshade_unsign order, proshade_double* abscissas, proshade_double* weights, proshade_double integralOverRange, proshade_double maxSphereDists )
{
    //================================================ Initialise local variables
    proshade_double ret                               = 0.0;
    proshade_complex* intData                         = new proshade_complex [order];
    ProSHADE_internal_misc::checkMemoryAllocation ( intData, __FILE__, __LINE__, __func__ );
    proshade_complex posVals;
    proshade_unsign lesserPos                         = 0;
    proshade_unsign upperPos                          = 0;
    proshade_double lesserWeight                      = 0.0;
    proshade_double upperWeight                       = 0.0;
    
    //================================================ Rescale to <order> points
    for ( proshade_unsign absIter = 0; absIter < order; absIter++ )
    {
        //============================================ Init loop
        posVals[0]                                    = 0.0;
        posVals[1]                                    = 0.0;
        
        //============================================ Find real position of abscissas
        posVals[0]                                    = ( ( abscissas[absIter] + 1.0 ) / 2.0 ) * integralOverRange;
        

        //============================================ Find lesser and upper bounds
        for ( proshade_unsign valIt = 0; valIt < valsSize; valIt++ )
        {
            if ( ( ( valIt * maxSphereDists ) <=  posVals[0] ) && ( ( ( valIt + 1 ) * maxSphereDists ) > posVals[0] ) )
            {
                lesserPos                             = static_cast<proshade_unsign> ( valIt );
                upperPos                              = static_cast<proshade_unsign> ( valIt + 1 );
                break;
            }
        }
        
        //============================================ Linear Interpolation
        lesserWeight                                  = 0.0;
        upperWeight                                   = 0.0;
        if ( lesserPos != 0 )
        {
            //======================================== Here we realise that the lesser and upper bounds were determined on scale 1 ... N, while our values are on scale 0 ... N-1 and therefore after determining the linear interpolation weights, we subtract 1 from both lesserPos and upperPos; however ...
            lesserWeight                              = upperPos - ( posVals[0] / maxSphereDists );
            upperWeight                               = 1.0 - lesserWeight;
                    
            posVals[1]                                = ( lesserWeight * vals[lesserPos-1] ) + ( upperWeight * vals[upperPos-1] );
        }
        else
        {
            //======================================== ... this then means that we would require position -1 for when the integration value is between 0 and the first shell. To resolve this, we assume that the values are 0 below the first shell and proceed as follows:
            upperWeight                               = 1.0 - ( upperPos - ( posVals[0] / maxSphereDists ) );
                    
            posVals[1]                                = ( upperWeight * vals[upperPos-1] );
        }
                
        intData[absIter][0]                           = posVals[0];
        intData[absIter][1]                           = posVals[1];
    }

    //================================================ Integrate
    for ( proshade_unsign absPoint = 0; absPoint < order; absPoint++ )
    {
        ret                                          += ( weights[absPoint] * intData[absPoint][1] );
    }
    
    //================================================ Normalise
    ret                                              *= ( integralOverRange / 2.0 );
    
    //================================================ Release memory
    delete[] intData;
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief Function to compute the complete complex Gauss-Legendre integration over spherical harmonic values in different shells.
 
    This function takes the real parts of the spherical harmonics value in different shells and proceeds to compute
    the Gauss-Legendre integration over them. It uses the shell positions to appropriately place abscissas and their
    weights, which it assumes were pre-computed by the getLegendreAbscAndWeights() function.
 
    \param[in] vals Pointer to a complex array of values over which the integration to be done.
    \param[in] valsSize The length of the input array.
    \param[in] order The integration order value.
    \param[in] abscissas The allocated array for holding the abscissa values.
    \param[in] weights The allocated array for holding the weight values.
    \param[in] integralOverRange The range of the intgral. If progressive shell mapping is used, this will not be max shell radius.
    \param[in] maxSphereDists Distance between two shells.
    \param[in] retReal The real part of the complex result of Gauss-Legendre integration over the shperical harmonics values.
    \param[in] retImag The imaginary part of the complex result of Gauss-Legendre integration over the shperical harmonics values.
 */
void ProSHADE_internal_maths::gaussLegendreIntegration ( proshade_complex* vals, proshade_unsign valsSize, proshade_unsign order, proshade_double* abscissas, proshade_double* weights, proshade_double integralOverRange, proshade_double maxSphereDists, proshade_double* retReal, proshade_double* retImag )
{
    //================================================ Initialise local variables
    proshade_triplet* intData                         = new proshade_triplet [order];
    ProSHADE_internal_misc::checkMemoryAllocation ( intData, __FILE__, __LINE__, __func__ );
    proshade_triplet posVals;
    proshade_unsign lesserPos                         = 0;
    proshade_unsign upperPos                          = 0;
    proshade_double lesserWeight                      = 0.0;
    proshade_double upperWeight                       = 0.0;
    
    //================================================ Rescale to <order> points
    for ( proshade_unsign absIter = 0; absIter < order; absIter++ )
    {
        //============================================ Init loop
        posVals[0]                                    = 0.0;
        posVals[1]                                    = 0.0;
        posVals[2]                                    = 0.0;
        
        //============================================ Find real position of abscissas
        posVals[0]                                    = ( ( abscissas[absIter] + 1.0 ) / 2.0 ) * integralOverRange;
        
        
        //============================================ Find lesser and upper bounds
        for ( proshade_unsign valIt = 0; valIt < valsSize; valIt++ )
        {
            if ( ( ( valIt * maxSphereDists ) <=  posVals[0] ) && ( ( ( valIt + 1 ) * maxSphereDists ) > posVals[0] ) )
            {
                lesserPos                             = static_cast<proshade_unsign> ( valIt );
                upperPos                              = static_cast<proshade_unsign> ( valIt + 1 );
                break;
            }
        }
        
        //============================================ Linear Interpolation
        lesserWeight                                  = 0.0;
        upperWeight                                   = 0.0;
        if ( lesserPos != 0 )
        {
            //======================================== Here we realise that the lesser and upper bounds were determined on scale 1 ... N, while our values are on scale 0 ... N-1 and therefore after determining the linear interpolation weights, we subtract 1 from both lesserPos and upperPos; however ...
            lesserWeight                              = upperPos - ( posVals[0] / maxSphereDists );
            upperWeight                               = 1.0 - lesserWeight;
                    
            posVals[1]                                = ( lesserWeight * vals[lesserPos-1][0] ) + ( upperWeight * vals[upperPos-1][0] );
            posVals[2]                                = ( lesserWeight * vals[lesserPos-1][1] ) + ( upperWeight * vals[upperPos-1][1] );
        }
        else
        {
            //======================================== ... this then means that we would require position -1 for when the integration value is between 0 and the first shell. To resolve this, we assume that the values are 0 below the first shell and proceed as follows:
            upperWeight                               = 1.0 - ( upperPos - ( posVals[0] / maxSphereDists ) );
                    
            posVals[1]                                = ( upperWeight * vals[upperPos-1][0] );
            posVals[2]                                = ( upperWeight * vals[upperPos-1][1] );
        }
                
        intData[absIter][0]                           = posVals[0];
        intData[absIter][1]                           = posVals[1];
        intData[absIter][2]                           = posVals[2];
    }
    
    //================================================ Integrate
   *retReal                                           = 0.0;
   *retImag                                           = 0.0;
    for ( proshade_unsign absPoint = 0; absPoint < order; absPoint++ )
    {
       *retReal                                      += ( weights[absPoint] * intData[absPoint][1] );
       *retImag                                      += ( weights[absPoint] * intData[absPoint][2] );
    }
    
    //================================================ Normalise
   *retReal                                          *= ( integralOverRange / 2.0 );
   *retImag                                          *= ( integralOverRange / 2.0 );
    
    //================================================ Release memory
    delete[] intData;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function to compute the complete complex matrix SVD and return only the sigmas.
 
    This function converts the input proshade_complex matrix of dimensions dim onto the LAPACK compatible
    std::complex<double> matrix. It then proceeds to create a dummy variables for the U and V matrices for
    saving the SVD results as well as other required variables. It finally proceeds to call LAPACK ZGESDD
    function to compute the SVD of the complex matrix input, checks the results and terminates. Note that
    this function does not make use of most of the LAPACK capabilities and is limitted onto square matrices.
 
    \param[in] mat Pointer to a complex square matrix with dimensions dim * dim.
    \param[in] dim The dimension of the complex matrix.
    \param[in] singularValues Empty array of size dim where the singular values will be saved.
 */
void ProSHADE_internal_maths::complexMatrixSVDSigmasOnly ( proshade_complex** mat, int dim, double*& singularValues )
{
    //================================================ Initialise local variables
    char job                                          = 'N';                                   // Save computation of parts of U and V matrices, they are not needed here
    std::complex<double> *rotMatU                     = new std::complex<double> [dim*dim];    // The U matrix space
    std::complex<double> *rotMatV                     = new std::complex<double> [dim*dim];    // The V^T matrix space
    std::complex<double> *work                        = new std::complex<double> [( 4 * dim)]; // Workspace, minimum required is 3*dim, using more for performance
    int workDim                                       = ( 4 * dim);                            // Formalism stating just that
    double* rwork                                     = new double[(7 * dim)];                 // Required by LAPACK, from 3.7 requires 7 * dim
    int* iwork                                        = new int[(8 * dim)];                    // Required by LAPACK
    int returnValue                                   = 0;                                     // This will tell if operation succeeded
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMatU, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMatV, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( work,    __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rwork,   __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( iwork,   __FILE__, __LINE__, __func__ );
    
    //================================================ Load input data into array in column-major order
    std::complex<double> *matrixToDecompose           = new std::complex<double>[dim*dim];
    ProSHADE_internal_misc::checkMemoryAllocation     ( matrixToDecompose, __FILE__, __LINE__, __func__ );
    for ( int rowIt = 0; rowIt < dim; rowIt++ )
    {
        for ( int colIt = 0; colIt < dim; colIt++ )
        {
            matrixToDecompose[(colIt*dim)+rowIt]      = std::complex<double> ( mat[rowIt][colIt][0], mat[rowIt][colIt][1] );
        }
    }
    
    //================================================ Run LAPACK ZGESDD
    zgesdd_                                           ( &job, &dim, &dim, matrixToDecompose, &dim, singularValues, rotMatU, &dim, rotMatV, &dim,
                                                        work, &workDim, rwork, iwork, &returnValue );
    
    //================================================ Free memory
    delete[] rotMatU;
    delete[] rotMatV;
    delete[] work;
    delete[] rwork;
    delete[] iwork;
    delete[] matrixToDecompose;
    
    //================================================ Check result
    if ( returnValue != 0 )
    {
        throw ProSHADE_exception ( "The LAPACK complex SVD algorithm did not converge!", "EL00021", __FILE__, __LINE__, __func__, "LAPACK algorithm for computing the singular value\n                    : decomposition of complex matrices did not converge and\n                    : therefore it was not possible to combine SH coefficients\n                    : from multiple shells. Changing the resolution may help,\n                    : contact me if this error persists." );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Function to compute the real matrix SVD and return the U and V matrices.
 
    This function converts the input proshade_double array of dimensions dim*dim onto the LAPACK compatible
    std::complex<double> matrix. It then proceeds to create a dummy variables for the U and V matrices for
    saving the SVD results as well as other required variables. It finally proceeds to call LAPACK ZGESDD
    function to compute the SVD of the real matrix input, checks the results and saves the U and V matrices
    for output. Note that this function does not make use of most of the LAPACK capabilities and is limitted
    onto square matrices.
 
    \param[in] mat Pointer to a real square matrix with dimensions dim * dim.
    \param[in] dim The dimension of the real matrix.
    \param[in] uAndV Empty and allocated array of size dim*6 where the U and V matrices will be saved.
    \param[in] fail If true and an error is encountered (typically algorithm not converging), this function will stop the program (useful for distances computations). However, if false, the function will simply return -777 as the first matrix element and not fail.
 */
void ProSHADE_internal_maths::complexMatrixSVDUandVOnly ( proshade_double* mat, int dim, proshade_double* uAndV, bool fail )
{
    //================================================ Initialise local variables
    char job                                          = 'A';                                 // Save computation of parts of U and V matrices, they are not needed here
    double* singularValues                            = new double[dim];                     // The array of singular values
    std::complex<double> *rotMatU                     = new std::complex<double> [dim*dim];  // The U matrix space
    std::complex<double> *rotMatV                     = new std::complex<double> [dim*dim];  // The V^T matrix space
    std::complex<double> *work                        = new std::complex<double> [static_cast<proshade_unsign>( ( 3 * dim) + pow( dim, 2 ) * dim)]; // Workspace, minimum required is 3*dim, using more for performance
    int workDim                                       = ( 3 * dim) + pow( dim, 2 );          // Formalism stating just that
    double* rwork                                     = new double[static_cast<proshade_unsign>((5 * dim) + 5 * pow(dim,2))]; // Required by LAPACK
    int* iwork                                        = new int[(8 * dim)];                  // Required by LAPACK
    int returnValue                                   = 0;                                   // This will tell if operation succeeded
    ProSHADE_internal_misc::checkMemoryAllocation     ( singularValues, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMatU,        __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMatV,        __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( work,           __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rwork,          __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( iwork,          __FILE__, __LINE__, __func__ );
    
    //================================================ Load input data into array in column-major order
    std::complex<double> *matrixToDecompose           = new std::complex<double>[dim*dim];
    ProSHADE_internal_misc::checkMemoryAllocation     ( matrixToDecompose, __FILE__, __LINE__, __func__ );
    for ( int rowIt = 0; rowIt < dim; rowIt++ )
    {
        for ( int colIt = 0; colIt < dim; colIt++ )
        {
            matrixToDecompose[(colIt*dim)+rowIt]      = std::complex<double> ( mat[(rowIt*dim)+colIt], 0.0 );
        }
    }
    
    //================================================ Run LAPACK ZGESDD
    zgesdd_                                           ( &job, &dim, &dim, matrixToDecompose, &dim, singularValues, rotMatU, &dim, rotMatV, &dim,
                                                        work, &workDim, rwork, iwork, &returnValue );
    
    //================================================ Free memory
    delete[] work;
    delete[] rwork;
    delete[] iwork;
    delete[] matrixToDecompose;
    delete[] singularValues;
    
    //================================================ Check result
    if ( ( returnValue != 0 ) && ( fail ) )
    {
        throw ProSHADE_exception ( "The LAPACK complex SVD algorithm did not converge!", "EL00022", __FILE__, __LINE__, __func__, "LAPACK algorithm for computing the singular value\n                    : decomposition of complex matrices did not converge and\n                    : therefore it was not possible to optimise the peak\n                    : positions in the (self-)rotation function. Changing the\n                    : resolution may help, contact me if this error persists." );
    }
    if ( ( returnValue != 0 ) && ( !fail ) )
    {
        uAndV[0]                                      =  -777.7;
        return ;
    }
    
    //================================================ Save U
    for ( proshade_signed rowIt = 0; rowIt < dim; rowIt++ )
    {
        for ( proshade_signed colIt = 0; colIt < dim; colIt++ )
        {
            uAndV[(rowIt*3)+colIt]                    = rotMatU[( rowIt * 3 ) + colIt].real();
        }
    }
    
    //================================================ Save V
    for ( proshade_signed rowIt = 0; rowIt < dim; rowIt++ )
    {
        for ( proshade_signed colIt = 0; colIt < dim; colIt++ )
        {
            uAndV[(rowIt*3)+colIt+9]                  = rotMatV[( rowIt * 3 ) + colIt].real();
        }
    }
    
    //================================================ Release the rest of the memory
    delete[] rotMatU;
    delete[] rotMatV;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function to find Euler angles (ZXZ convention) from index position in the inverse SOFT map.
 
    This function proceeds to convert the inverse SOFT map x, y and z position to Euler ZXZ convention angles, saving
    these into the inputted pointers. It also changes the Euler angle ranges from (0,2PI> to (-PI,PI> as preferred by
    the functions using the results.
 
    \param[in] band The maximum bandwidth of the computation.
    \param[in] x The x-axis position in the inverse SOFT map.
    \param[in] y The x-axis position in the inverse SOFT map.
    \param[in] z The x-axis position in the inverse SOFT map.
    \param[in] eulerAlpha Pointer to where the Euler alpha angle will be saved.
    \param[in] eulerBeta Pointer to where the Euler beta angle will be saved.
    \param[in] eulerGamma Pointer to where the Euler gamma angle will be saved.
 */
void ProSHADE_internal_maths::getEulerZXZFromSOFTPosition ( proshade_signed band, proshade_signed x, proshade_signed y, proshade_signed z, proshade_double* eulerAlpha, proshade_double* eulerBeta, proshade_double* eulerGamma )
{
    //================================================ Convert index to Euler angles
   *eulerGamma                                        = ( M_PI * y / ( static_cast<proshade_double> ( band ) ) );
   *eulerBeta                                         = ( M_PI * ( 2.0 * x + 1.0 ) / static_cast<proshade_double> ( 4.0 * band ) )  ;
   *eulerAlpha                                        = ( M_PI * z / ( static_cast<proshade_double> ( band ) ) );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function to find the rotation matrix from Euler angles (ZXZ convention).
 
    \param[in] eulerAlpha The Euler alpha angle value.
    \param[in] eulerBeta The Euler beta angle value.
    \param[in] eulerGamma The Euler gamma angle value.
    \param[in] matrix A pointer to array of 9 values to which the results of the function will be saved.
 */
void ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( proshade_double eulerAlpha, proshade_double eulerBeta, proshade_double eulerGamma, proshade_double* matrix )
{
    //================================================ First row
    matrix[0]                                         =  cos ( eulerAlpha ) * cos ( eulerBeta  ) * cos ( eulerGamma ) - sin ( eulerAlpha ) * sin ( eulerGamma );
    matrix[1]                                         =  sin ( eulerAlpha ) * cos ( eulerBeta  ) * cos ( eulerGamma ) + cos ( eulerAlpha ) * sin ( eulerGamma );
    matrix[2]                                         = -sin ( eulerBeta  ) * cos ( eulerGamma );
  
    //================================================ Second row
    matrix[3]                                         = -cos ( eulerAlpha ) * cos ( eulerBeta  ) * sin ( eulerGamma ) - sin ( eulerAlpha ) * cos ( eulerGamma );
    matrix[4]                                         = -sin ( eulerAlpha ) * cos ( eulerBeta  ) * sin ( eulerGamma ) + cos ( eulerAlpha ) * cos ( eulerGamma );
    matrix[5]                                         =  sin ( eulerBeta  ) * sin ( eulerGamma );
  
    //================================================ Third row
    matrix[6]                                         =  cos ( eulerAlpha ) * sin ( eulerBeta  );
    matrix[7]                                         =  sin ( eulerAlpha ) * sin ( eulerBeta  );
    matrix[8]                                         =  cos ( eulerBeta  );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function converts rotation matrix to the axis-angle representation.
 
    This function takes a rotation matrix as an array of 9 numbers and converts it to the Angle-Axis representation,
    which is the main rotation representation used in ProSHADE. This function deals with both the North and South
    pole singularity of the rotation matrices.
 
    \param[in] rotMat Rotation matrix as an array of 9 values.
    \param[in] x Pointer to which the x-axis value of the axis vector will be saved.
    \param[in] y Pointer to which the y-axis value of the axis vector will be saved.
    \param[in] z Pointer to which the z-axis value of the axis vector will be saved.
    \param[in] ang Pointer to which the angle value will be saved.
 */
 void ProSHADE_internal_maths::getAxisAngleFromRotationMatrix ( proshade_double* rotMat, proshade_double* x, proshade_double* y, proshade_double* z, proshade_double* ang )
{
    //================================================ Initialise
    proshade_double singAtPiCheck                     = 0.01;
    proshade_double singAtIdentity                    = 0.05;
    
    //================================================ Check input for singularities
    if ( ( std::abs ( rotMat[1] - rotMat[3] ) < singAtPiCheck ) &&
         ( std::abs ( rotMat[2] - rotMat[6] ) < singAtPiCheck ) &&
         ( std::abs ( rotMat[5] - rotMat[7] ) < singAtPiCheck ) )
    {
        //============================================ Singularity in input! Check for identity matrix
        if ( ( std::abs ( rotMat[1] + rotMat[3] ) < singAtIdentity ) &&
             ( std::abs ( rotMat[2] + rotMat[6] ) < singAtIdentity ) &&
             ( std::abs ( rotMat[5] + rotMat[7] ) < singAtIdentity ) &&
             ( std::abs ( rotMat[0] + rotMat[4] + rotMat[8] - 3.0 ) < singAtIdentity ) )
        {
            //======================================== Identity matrix. Return 0 angle.
           *x                                         = 1.0;
           *y                                         = 0.0;
           *z                                         = 0.0;
           *ang                                       = 0.0;
            
            //======================================== Done
            return ;
        }
        
        //============================================ If we got here, this is the 180deg (pi rad) singularity. Find which axis should the rotation be done along
       *ang                                           = M_PI;
                
        proshade_double xx                            = ( rotMat[0] + 1.0 ) / 2.0;
        proshade_double yy                            = ( rotMat[4] + 1.0 ) / 2.0;
        proshade_double zz                            = ( rotMat[8] + 1.0 ) / 2.0;
        proshade_double xy                            = ( rotMat[1] + rotMat[3] ) / 4.0;
        proshade_double xz                            = ( rotMat[2] + rotMat[6] ) / 4.0;
        proshade_double yz                            = ( rotMat[5] + rotMat[7] ) / 4.0;
        
        if ( ( xx > yy ) && ( xx > zz ) ) // XX is the largest diagonal
        {
            if ( xx < singAtPiCheck ) // and is still 0
            {
               *x                                     = 0.0;
               *y                                     = 1.0 / sqrt(2);
               *z                                     = 1.0 / sqrt(2);
            }
            else
            {
               *x                                     =  sqrt ( xx );
               *y                                     =  xy / sqrt ( xx );
               *z                                     =  xz / sqrt ( xx );
            }
        }
        
        else if ( yy > zz ) // YY is the largest diagonal
        {
            if ( yy < singAtPiCheck ) // and is still 0
            {
               *x                                     =  1.0 / sqrt(2);
               *y                                     =  0.0;
               *z                                     =  1.0 / sqrt(2);
            }
            else
            {
               *y                                     =  sqrt ( yy );
               *x                                     =  xy / sqrt ( yy );
               *z                                     =  yz / sqrt ( yy );
            }
        }
        
        else // ZZ is the largest diagonal
        {
            if ( zz < singAtPiCheck ) // and is still 0
            {
               *x                                     = 1.0 / sqrt(2);
               *y                                     = 1.0 / sqrt(2);
               *z                                     = 0.0;
            }
            else
            {
               *z                                     = sqrt ( zz );
               *x                                     = xz / sqrt ( zz );
               *y                                     = yz / sqrt ( zz );
            }
        }
        
        //============================================ Done
        return ;
    }
    
    //================================================ No singularities! Now get angle
   *ang                                               = std::acos ( ( std::max ( -1.0, std::min ( 3.0, rotMat[0] + rotMat[4] + rotMat[8] ) ) - 1.0 ) / 2.0 );
    
    //================================================ Init return values
   *x                                                 = 1.0;
   *y                                                 = 0.0;
   *z                                                 = 0.0;
    
    //================================================ Is angle 0? This should not happen, but will
    if ( std::abs ( *ang ) < singAtPiCheck )
    {
       *ang                                           = 0.0;
        return ;
    }
    
    //================================================ Axis
   *x                                                 = rotMat[7] - rotMat[5];
   *y                                                 = rotMat[2] - rotMat[6];
   *z                                                 = rotMat[3] - rotMat[1];
    proshade_double normFactor                        = pow ( *x, 2.0 ) + pow ( *y, 2.0 ) + pow ( *z, 2.0 );
            
    if ( normFactor > singAtPiCheck )
    {
        normFactor                                    = sqrt ( normFactor );
       *x                                            /= normFactor;
       *y                                            /= normFactor;
       *z                                            /= normFactor;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function converts the axis-angle representation to the rotation matrix representation.
 
    \param[in] rotMat Rotation matrix as an array of 9 values will be saved to this pointer, must already be allocated.
    \param[in] x The x-axis value of the axis vector.
    \param[in] y The y-axis value of the axis vector.
    \param[in] z The z-axis value of the axis vector.
    \param[in] angThe angle value.
 */
void ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( proshade_double* rotMat, proshade_double x, proshade_double y, proshade_double z, proshade_double ang )
{
    //================================================ If angle is 0 or infinity (anything divided by 0), return identity matrix
    if ( ( ang == 0.0 ) || ( std::isinf ( ang ) ) )
    {
        //============================================ Create identity
        for ( proshade_unsign i = 0; i < 9; i++ ) { rotMat[i] = 0.0; }
        rotMat[0]                                     = 1.0;
        rotMat[4]                                     = 1.0;
        rotMat[8]                                     = 1.0;
        
        //============================================ Done
        return ;
    }
    
    //================================================ Compute the matrix
    proshade_double cAng                              = cos ( ang );
    proshade_double sAng                              = sin ( ang );
    proshade_double tAng                              = 1.0 - cAng;
            
    rotMat[0]                                         = cAng + x * x * tAng;
    rotMat[4]                                         = cAng + y * y * tAng;
    rotMat[8]                                         = cAng + z * z * tAng;
            
    proshade_double tmp1                              = x * y * tAng;
    proshade_double tmp2                              = z * sAng;
    rotMat[3]                                         = tmp1 + tmp2;
    rotMat[1]                                         = tmp1 - tmp2;
            
    tmp1                                              = x * z * tAng;
    tmp2                                              = y * sAng;
    rotMat[6]                                         = tmp1 - tmp2;
    rotMat[2]                                         = tmp1 + tmp2;
            
    tmp1                                              = y * z * tAng;
    tmp2                                              = x * sAng;
    rotMat[7]                                         = tmp1 + tmp2;
    rotMat[5]                                         = tmp1 - tmp2;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function converts rotation matrix to the Euler ZXZ angles representation.
 
    \param[in] rotMat Rotation matrix as an array of 9 values.
    \param[in] eA Pointer to which the Euler angle alpha value will be saved.
    \param[in] eB Pointer to which the Euler angle beta value will be saved.
    \param[in] eG Pointer to which the Euler angle gamma value will be saved.
 */
void ProSHADE_internal_maths::getEulerZXZFromRotMatrix ( proshade_double* rotMat, proshade_double* eA, proshade_double* eB, proshade_double* eG )
{
    //================================================ Get ZXZ Euler from matrix
   *eA                                                = atan2 ( rotMat[7],  rotMat[6] );
   *eB                                                = acos  ( rotMat[8] );
   *eG                                                = atan2 ( rotMat[5], -rotMat[2] );
    
    //================================================ Solve undefined 0,0 inputs (i.e. identity matrix)
    proshade_double errLimit                          = 0.001;
    if ( ( ( rotMat[7] < errLimit ) && ( rotMat[7] > -errLimit ) ) && ( ( rotMat[6] < errLimit ) && ( rotMat[6] > -errLimit ) ) )
    {
        //============================================ atan2 (0,0) is undefined, we want 0.0 here
       *eA                                            = 0.0;
    }
    
    if ( ( ( rotMat[5] < errLimit ) && ( rotMat[5] > -errLimit ) ) && ( ( rotMat[2] < errLimit ) && ( rotMat[2] > -errLimit ) ) )
    {
        //============================================ atan2 (0,0) is undefined, we want 0.0 here
       *eG                                            = 0.0;
    }
    
    //================================================ Get the angles to proper range
    if ( *eA < 0.0 ) { *eA                            = 2.0 * M_PI + *eA; }
    if ( *eB < 0.0 ) { *eB                            =       M_PI + *eB; }
    if ( *eG < 0.0 ) { *eG                            = 2.0 * M_PI + *eG; }
    
    //================================================ Done
    return ;
    
}

/*! \brief Function to compute matrix multiplication.
 
    \param[in] A The left matrix of the matrix multiplication to be solved.
    \param[in] B The right matrix of the matrix multiplication to be solved. (Assuming it already has been transposed).
    \param[in] res Matrix containing the results.
    \param[in] dim The dimension of all the matrices (i.e. assuming square dim*dim matrices).
 
    \warning This function assumes the second matrix has been transposed already!
 */
void ProSHADE_internal_maths::multiplyTwoSquareMatrices ( proshade_double* A, proshade_double* B, proshade_double* res, proshade_unsign dim )
{
    //================================================ Compute the matrix multiplication
    for ( proshade_unsign row = 0; row < dim; row++ )
    {
        for ( proshade_unsign col  = 0; col < dim; col++ )
        {
            for ( proshade_unsign inner = 0; inner < dim; inner++ )
            {
                res[(row*dim)+col]                   += A[(inner*dim)+row] * B[(col*dim)+inner];
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Function to find prime factors of an integer.
 
    \param[in] number A single integer number to be decomposed into its prime factors.
 */
std::vector < proshade_signed > ProSHADE_internal_maths::primeFactorsDecomp ( proshade_signed number )
{
    //================================================ Initialise variables
    std::vector < proshade_signed > ret;
    
    //================================================ Deal with negative numbers
    bool changeSign                                   = false;
    if ( number < 0 ) { changeSign = true; number = -number; }
    
    //================================================ Deal with zero and one
    if ( number == 0 ) { ProSHADE_internal_misc::addToSignedVector ( &ret, 0 ); return ( ret ); }
    if ( number == 1 ) { ProSHADE_internal_misc::addToSignedVector ( &ret, 1 ); return ( ret ); }
    
    //================================================ Divide by 2 as long as you can
    while ( number % 2 == 0 )
    {
        ProSHADE_internal_misc::addToSignedVector     ( &ret, 2 );
        number                                        = number / 2;
    }
    
    //================================================ Check all odd numbers up to the square root
    for ( proshade_unsign posDiv = 3; posDiv <= sqrt ( number ); posDiv += 2)
    {
        // If posDiv is a divisor of the number, save the result
        while ( number % posDiv == 0 )
        {
            ProSHADE_internal_misc::addToSignedVector ( &ret, posDiv );
            number                                    = number / posDiv;
        }
    }
    
    //================================================ If the number was a large prime number, save it as it is
    if ( number > 2 ) { ProSHADE_internal_misc::addToSignedVector ( &ret, number ); }

    //================================================ Finish dealing with negative numbers
    if ( changeSign ) { ret.at(0) = -ret.at(0); }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief Function to the heiht of normal distribution given by mean and standard deviation for a given value.
 
    \param[in] mean The mean of the normal distribution.
    \param[in] standardDev The standard deviation of the normal distribution.
    \param[in] value The value on the axis for which the height of the normal distribution is to be obtained.
    \param[out] X The height of the normal distribution at point given by the value.
 */
proshade_double ProSHADE_internal_maths::normalDistributionValue ( proshade_double mean, proshade_double standardDev, proshade_double value )
{
    //================================================ Compute and return
    return                                            ( ( 1.0 / sqrt ( 2.0 * M_PI * pow(standardDev,2.0) ) ) * std::exp ( - pow( value - mean, 2.0 ) / 2.0 * pow(standardDev,2.0) ) );
    
}

/*! \brief Simple 3D vector dot product computation.
 
    \param[in] x1 The x-axis element of the first vector.
    \param[in] y1 The y-axis element of the first vector.
    \param[in] z1 The z-axis element of the first vector.
    \param[in] x2 The x-axis element of the second vector.
    \param[in] y2 The y-axis element of the second vector.
    \param[in] z2 The z-axis element of the second vector.
 */
proshade_double ProSHADE_internal_maths::computeDotProduct ( proshade_double* x1, proshade_double* y1, proshade_double* z1, proshade_double* x2, proshade_double* y2, proshade_double* z2 )
{
    //================================================ Compute and return
    return                                            ( (*x1 * *x2) + (*y1 * *y2) + (*z1 * *z2) );
}

/*! \brief Function for finding a vector which would have a given two dot products to two other vectors.

    This function takes two vectors and two dot product values. It then basically solves the following set of equations for x, y and z:
    
        solX*x1 + solY*y1 + solZ*z1 = dot1
        solX*x2 + solY*y2 + solZ*z2 = dot2
        sqrt ( solX^2 + solY^2 + solZ^2 ) = 1.0
 
    This should result in a vector, which has the required angles to both input vectors and is normalised (this is done by the third equation, which is required to obtain system of three equations with three unkonwns). The equations are courtesy of https://www.wolframalpha.com/input/?i=Solve%5B%7Ba+x+%2B+b+y+%2B+c+z+%3D%3D+f%2C+k+x+%2B+l+y+%2B+m+z+%3D%3D+g%2C+Sqrt%5Bx%5E2+%2B+y%5E2+%2B+z%5E2%5D+%3D%3D+1%7D%2C+%7Bx%2C+y%2C+z%7D%5D
    webpage of Wolfram Alpha. If in doubt, do not fear to derive yourself :-).
 
    \param[in] x1 The x-axis element of the first vector.
    \param[in] y1 The y-axis element of the first vector.
    \param[in] z1 The z-axis element of the first vector.
    \param[in] x2 The x-axis element of the second vector.
    \param[in] y2 The y-axis element of the second vector.
    \param[in] z2 The z-axis element of the second vector.
    \param[in] dot1 The dot product specifying the angle between the sought vector and the first input vector.
    \param[in] dot2 The dot product specifying the angle between the sought vector and the second input vectors.
    \param[out] vec A std::vector containing the three elements of the sought vector.
*/
std::vector < proshade_double > ProSHADE_internal_maths::findVectorFromTwoVAndTwoD ( proshade_double x1, proshade_double y1, proshade_double z1, proshade_double x2, proshade_double y2, proshade_double z2, proshade_double dot1, proshade_double dot2 )
{
    //================================================ Initialise variables
    std::vector < proshade_double > ret;
    
    //================================================ Solution
    proshade_double solX                              = ( -sqrt ( pow ( 2.0 * x1 * y1 * dot2 * y2 + 2.0 * x1 * z1 * dot2 * z2 - 2.0 * x1 * dot1 * pow ( y2, 2.0 ) - 2.0 * x1 * dot1 * pow ( z2, 2.0 ) - 2.0 * pow ( y1, 2.0 ) * dot2 * x2 + 2.0 * y1 * dot1 * x2 * y2 - 2.0 * pow ( z1, 2.0 ) * dot2 * x2 + 2.0 * z1 * dot1 * x2 * z2, 2.0 ) -
                                                                  4.0 * ( pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * x1 * y1 * x2 * y2 - 2.0 * x1 * z1 * x2 * z2 + pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) ) *
                                                                  ( pow ( y1, 2.0 ) * pow ( dot2, 2.0 ) - pow ( y1, 2.0 ) * pow ( z2, 2.0 ) + 2.0 * y1 * z1 * y2 * z2 - 2.0 * y1 * dot1 * dot2 * y2 + pow ( z1, 2.0 ) * pow ( dot2, 2.0 ) - pow ( z1, 2.0 ) * pow ( y2, 2.0 ) - 2.0 * z1 * dot1 * dot2 * z2 + pow ( dot1, 2.0 ) * pow ( y2, 2.0 ) + pow ( dot1, 2.0 ) * pow ( z2, 2.0 ) ) ) -
                                                          2.0 * x1 * y1 * dot2 * y2 - 2.0 * x1 * z1 * dot2 * z2 + 2.0 * x1 * dot1 * pow ( y2, 2.0 ) + 2.0 * x1 * dot1 * pow ( z2, 2.0 ) + 2.0 * pow ( y1, 2.0 ) * dot2 * x2 - 2.0 * y1 * dot1 * x2 * y2 + 2.0 * pow ( z1, 2.0 ) * dot2 * x2 - 2.0 * z1 * dot1 * x2 * z2 ) /
                                                        ( 2.0 * ( pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * x1 * y1 * x2 * y2 - 2.0 * x1 * z1 * x2 * z2 + pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) ) );
    proshade_double solY                              = ( ( dot2 * pow ( x2, 2.0 ) * pow ( z1, 3.0 ) ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) -
                                                          ( dot1 * pow ( x2, 2.0 ) * z2 * pow ( z1, 2.0 ) ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) -
                                                          ( 2.0 * x1 * dot2 * x2 * z2 * pow ( z1, 2.0 ) ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) - dot2 * z1 -
                                                          ( x2 * sqrt ( pow ( - 2.0 * dot2 * x2 * pow ( y1, 2.0 ) + 2.0 * x1 * dot2 * y2 * y1 + 2.0 * dot1 * x2 * y2 * y1 - 2.0 * x1 * dot1 * pow ( y2, 2.0 ) - 2.0 * x1 * dot1 * pow ( z2, 2.0 ) - 2.0 * pow ( z1, 2.0 ) * dot2 * x2 + 2.0 * x1 * z1 * dot2 * z2 + 2.0 * z1 * dot1 * x2 * z2, 2.0 ) -
                                                                        4.0 * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) *
                                                                        ( pow ( y1, 2.0 ) * pow ( dot2, 2.0 ) + pow ( z1, 2.0 ) * pow ( dot2, 2.0 ) - 2.0 * y1 * dot1 * y2 * dot2 - 2.0 * z1 * dot1 * z2 * dot2 - pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( dot1, 2.0 ) * pow ( y2, 2.0 ) - pow ( y1, 2.0 ) * pow ( z2, 2.0 ) + pow ( dot1, 2.0 ) * pow ( z2, 2.0 ) + 2.0 * y1 * z1 * y2 * z2 ) ) * z1 ) /
                                                          ( 2.0 * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) +
                                                          ( pow ( y1, 2.0 ) * dot2 * pow ( x2, 2.0 ) * z1 ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) +
                                                          ( x1 * dot1 * x2 * pow ( y2, 2.0 ) * z1 ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) +
                                                          ( pow ( x1, 2.0 ) * dot2 * pow ( z2, 2.0 ) * z1 ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) +
                                                          ( 2.0 * x1 * dot1 * x2 * pow ( z2, 2.0 ) * z1 ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) -
                                                          ( y1 * dot1 * pow ( x2, 2.0 ) * y2 * z1 ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) -
                                                          ( x1 * y1 * dot2 * x2 * y2 * z1 ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) + dot1 * z2 +
                                                          ( x1 * z2 * sqrt ( pow ( -2.0 * dot2 * x2 * pow ( y1, 2.0 ) + 2.0 * x1 * dot2 * y2 * y1 + 2.0 * dot1 * x2 * y2 * y1 - 2.0 * x1 * dot1 * pow ( y2, 2.0 ) - 2.0 * x1 * dot1 * pow ( z2, 2.0 ) - 2.0 * pow ( z1, 2.0 ) * dot2 * x2 + 2.0 * x1 * z1 * dot2 * z2 + 2.0 * z1 * dot1 * x2 * z2, 2.0 ) -
                                                                             4.0 * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) *
                                                                             ( pow ( y1, 2.0 ) * pow ( dot2, 2.0 ) + pow ( z1, 2.0 ) * pow ( dot2, 2.0 ) - 2.0 * y1 * dot1 * y2 * dot2 - 2.0 * z1 * dot1 * z2 * dot2 - pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( dot1, 2.0 ) * pow ( y2, 2.0 ) - pow ( y1, 2.0 ) * pow ( z2, 2.0 ) + pow ( dot1, 2.0 ) * pow ( z2, 2.0 ) + 2.0 * y1 * z1 * y2 * z2 ) ) ) /
                                                          ( 2.0 * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) -
                                                          ( pow ( x1, 2.0 ) * dot1 * pow ( z2, 3.0 ) ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) -
                                                          ( pow ( x1, 2.0 ) * dot1 * pow ( y2, 2.0 ) * z2 ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) -
                                                          ( x1 * pow ( y1, 2.0 ) * dot2 * x2 * z2 ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) +
                                                          ( pow ( x1, 2.0 ) * y1 * dot2 * y2 * z2 ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) +
                                                          ( x1 * y1 * dot1 * x2 * y2 * z2 ) /
                                                          ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) / ( y1 * z2 - z1 * y2 );
    proshade_double solZ                              = ( - ( dot2 * pow ( x2, 2.0 ) * y2 * pow ( z1, 3.0 ) ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) -
                                                            ( dot2 * pow ( x2, 2.0 ) * pow ( z1, 2.0 ) ) /
                                                            ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) +
                                                            ( dot1 * pow ( x2, 2.0 ) * y2 * z2 * pow ( z1, 2.0 ) ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) +
                                                            ( 2.0 * x1 * dot2 * x2 * y2 * z2 * pow ( z1, 2.0 ) ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) +
                                                            ( x2 * y2 * sqrt ( pow ( -2.0 * dot2 * x2 * pow ( y1, 2.0 ) + 2.0 * x1 * dot2 * y2 * y1 + 2.0 * dot1 * x2 * y2 * y1 - 2.0 * x1 * dot1 * pow ( y2, 2.0 ) - 2.0 * x1 * dot1 * pow ( z2, 2.0 ) - 2.0 * pow ( z1, 2.0 ) * dot2 * x2 + 2.0 * x1 * z1 * dot2 * z2 + 2.0 * z1 * dot1 * x2 * z2, 2.0 ) -
                                                                              4.0 * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) *
                                                                              ( pow ( y1, 2.0 ) * pow ( dot2, 2.0 ) + pow ( z1, 2.0 ) * pow ( dot2, 2.0 ) - 2.0 * y1 * dot1 * y2 * dot2 - 2.0 * z1 * dot1 * z2 * dot2 - pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( dot1, 2.0 ) * pow ( y2, 2.0 ) - pow ( y1, 2.0 ) * pow ( z2, 2.0 ) + pow ( dot1, 2.0 ) * pow ( z2, 2.0 ) + 2.0 * y1 * z1 * y2 * z2 ) ) * z1 ) /
                                                            ( 2.0 * ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) +
                                                            ( dot2 * y2 * z1 ) / ( y1 * z2 - z1 * y2 ) +
                                                            ( dot1 * pow ( x2, 2.0 ) * z2 * z1 ) /
                                                            ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) +
                                                            ( x1 * dot2 * x2 * z2 * z1 ) /
                                                            ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) -
                                                            ( x1 * dot1 * x2 * pow ( y2, 3.0 ) * z1 ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) +
                                                            ( y1 * dot1 * pow ( x2, 2.0 ) * pow ( y2, 2.0 ) * z1 ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) +
                                                            ( x1 * y1 * dot2 * x2 * pow ( y2, 2.0 ) * z1 ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) -
                                                            ( pow ( x1, 2.0 ) * dot2 * y2 * pow ( z2, 2.0 ) * z1 ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) -
                                                            ( 2.0 * x1 * dot1 * x2 * y2 * pow ( z2, 2.0 ) * z1 ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) -
                                                            ( pow ( y1, 2.0 ) * dot2 * pow ( x2, 2.0 ) * y2 * z1 ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) + dot2 +
                                                            ( x2 * sqrt ( pow ( - 2.0 * dot2 * x2 * pow ( y1, 2.0 ) + 2.0 * x1 * dot2 * y2 * y1 + 2.0 * dot1 * x2 * y2 * y1 - 2.0 * x1 * dot1 * pow ( y2, 2.0 ) - 2.0 * x1 * dot1 * pow ( z2, 2.0 ) - 2.0 * pow ( z1, 2.0 ) * dot2 * x2 + 2.0 * x1 * z1 * dot2 * z2 + 2.0 * z1 * dot1 * x2 * z2, 2.0 ) -
                                                                          4.0 * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) *
                                                                          ( pow ( y1, 2.0 ) * pow ( dot2, 2.0 ) + pow ( z1, 2.0 ) * pow ( dot2, 2.0 ) - 2.0 * y1 * dot1 * y2 * dot2 - 2.0 * z1 * dot1 * z2 * dot2 - pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( dot1, 2.0 ) * pow ( y2, 2.0 ) - pow ( y1, 2.0 ) * pow ( z2, 2.0 ) + pow ( dot1, 2.0 ) * pow ( z2, 2.0 ) + 2.0 * y1 * z1 * y2 * z2 ) ) ) /
                                                            ( 2.0 * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) -
                                                            ( x1 * y2 * z2 * sqrt ( pow ( - 2.0 * dot2 * x2 * pow ( y1, 2.0 ) + 2.0 * x1 * dot2 * y2 * y1 + 2.0 * dot1 * x2 * y2 * y1 - 2.0 * x1 * dot1 * pow ( y2, 2.0 ) - 2.0 * x1 * dot1 * pow ( z2, 2.0 ) - 2.0 * pow ( z1, 2.0 ) * dot2 * x2 + 2.0 * x1 * z1 * dot2 * z2 + 2.0 * z1 * dot1 * x2 * z2, 2.0 ) -
                                                                                    4.0 * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) *
                                                                                    ( pow ( y1, 2.0 ) * pow ( dot2, 2.0 ) + pow ( z1, 2.0 ) * pow ( dot2, 2.0 ) - 2.0 * y1 * dot1 * y2 * dot2 - 2.0 * z1 * dot1 * z2 * dot2 - pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( dot1, 2.0 ) * pow ( y2, 2.0 ) - pow ( y1, 2.0 ) * pow ( z2, 2.0 ) + pow ( dot1, 2.0 ) * pow ( z2, 2.0 ) + 2.0 * y1 * z1 * y2 * z2 ) ) ) /
                                                            ( 2.0 * ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) -
                                                            ( dot1 * y2 * z2 ) / ( y1 * z2 - z1 * y2 ) -
                                                            ( pow ( y1, 2.0 ) * dot2 * pow ( x2, 2.0 ) ) /
                                                            ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) -
                                                            ( x1 * dot1 * x2 * pow ( y2, 2.0 ) ) /
                                                            ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) -
                                                            ( x1 * dot1 * x2 * pow ( z2, 2.0 ) ) /
                                                            ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) +
                                                            ( y1 * dot1 * pow ( x2, 2.0 ) * y2 ) /
                                                            ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) +
                                                            ( x1 * y1 * dot2 * x2 * y2 ) /
                                                            ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) +
                                                            ( pow ( x1, 2.0 ) * dot1 * y2 * pow ( z2, 3.0 ) ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) +
                                                            ( pow ( x1, 2.0 ) * dot1 * pow ( y2, 3.0 ) * z2 ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) -
                                                            ( pow ( x1, 2.0 ) * y1 * dot2 * pow ( y2, 2.0 ) * z2 ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) -
                                                            ( x1 * y1 * dot1 * x2 * pow ( y2, 2.0 ) * z2 ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) +
                                                            ( x1 * pow ( y1, 2.0 ) * dot2 * x2 * y2 * z2 ) /
                                                            ( ( y1 * z2 - z1 * y2 ) * ( pow ( y1, 2.0 ) * pow ( x2, 2.0 ) + pow ( z1, 2.0 ) * pow ( x2, 2.0 ) - 2.0 * x1 * y1 * y2 * x2 - 2.0 * x1 * z1 * z2 * x2 + pow ( x1, 2.0 ) * pow ( y2, 2.0 ) + pow ( z1, 2.0 ) * pow ( y2, 2.0 ) + pow ( x1, 2.0 ) * pow ( z2, 2.0 ) + pow ( y1, 2.0 ) * pow ( z2, 2.0 ) - 2.0 * y1 * z1 * y2 * z2 ) ) ) / z2;

    //================================================ Set largest axis element to positive (ProSHADE standard)
    if ( ( ( std::max ( std::abs ( solX ), std::max( std::abs ( solY ), std::abs ( solZ ) ) ) == std::abs ( solX ) ) && ( solX < 0.0 ) ) ||
         ( ( std::max ( std::abs ( solX ), std::max( std::abs ( solY ), std::abs ( solZ ) ) ) == std::abs ( solY ) ) && ( solY < 0.0 ) ) ||
         ( ( std::max ( std::abs ( solX ), std::max( std::abs ( solY ), std::abs ( solZ ) ) ) == std::abs ( solZ ) ) && ( solZ < 0.0 ) ) ) { solX *= -1.0; solY *= -1.0; solZ *= -1.0; }
    
    //================================================ Save solutions
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, solX  );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, solY  );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, solZ  );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief Function for finding a vector which would have a given three dot products to three other vectors.

    This function takes three vectors and three dot product values. It then basically solves the following set of equations for x, y and z:
    
        solX*x1 + solY*y1 + solZ*z1 = dot1
        solX*x2 + solY*y2 + solZ*z2 = dot2
        solX*x3 + solY*y3 + solZ*z3 = dot3
 
    This should result in a vector, which has the required angles to all input vectors and is normalised (this is done later as part of this function). The equations are courtesy of https://www.wolframalpha.com/input/?i=Solve%5B%7Ba+x+%2B+b+y+%2B+c+z+%3D%3D+f%2C+u+x+%2B+v+y+%2B+w+z+%3D%3D+g%2C+k+x+%2B+l+y+%2B+m+z+%3D%3D+h%7D%2C+%7Bx%2C+y%2C+z%7D%5D
    webpage of Wolfram Alpha. If in doubt, do not fear to derive yourself :-).
 
    \param[in] x1 The x-axis element of the first vector.
    \param[in] y1 The y-axis element of the first vector.
    \param[in] z1 The z-axis element of the first vector.
    \param[in] x2 The x-axis element of the second vector.
    \param[in] y2 The y-axis element of the second vector.
    \param[in] z2 The z-axis element of the second vector.
    \param[in] dot1 The dot product specifying the angle between the sought vector and the first input vector.
    \param[in] dot2 The dot product specifying the angle between the sought vector and the second input vectors.
    \param[out] vec A std::vector containing the three elements of the sought vector.
*/
std::vector < proshade_double > ProSHADE_internal_maths::findVectorFromThreeVAndThreeD ( proshade_double x1, proshade_double y1, proshade_double z1, proshade_double x2, proshade_double y2, proshade_double z2, proshade_double x3, proshade_double y3, proshade_double z3, proshade_double dot1, proshade_double dot2, proshade_double dot3 )
{
    //================================================ Initialise variables
    std::vector < proshade_double > ret;
    
    //================================================ Solution
    proshade_double solX                              = - (  y1 * dot2 * z3 - y1 * dot3 * z2 - z1 * dot2 * y3 + z1 * dot3 * y2 + dot1 * y3 * z2 - dot1 * z3 * y2 ) /
                                                          ( -x1 * y3 * z2 + x1 * z3 * y2 + y1 * x3 * z2 - y1 * z3 * x2 - z1 * x3 * y2 + z1 * y3 * x2 );
    proshade_double solY                              = - (  x1 * dot2 * z3 - x1 * dot3 * z2 - z1 * dot2 * x3 + z1 * dot3 * x2 + dot1 * x3 * z2 - dot1 * z3 * x2 ) /
                                                          (  x1 * y3 * z2 - x1 * z3 * y2 - y1 * x3 * z2 + y1 * z3 * x2 + z1 * x3 * y2 - z1 * y3 * x2 );
    proshade_double solZ                              = - (  x1 * dot2 * y3 - x1 * dot3 * y2 - y1 * dot2 * x3 + y1 * dot3 * x2 + dot1 * x3 * y2 - dot1 * y3 * x2 ) /
                                                          ( -x1 * y3 * z2 + x1 * z3 * y2 + y1 * x3 * z2 - y1 * z3 * x2 - z1 * x3 * y2 + z1 * y3 * x2 );

    //================================================ Normalise the axis to magnitude 1
    proshade_double normFactor                        = sqrt ( pow ( solX, 2.0 ) + pow ( solY, 2.0 ) + pow ( solZ, 2.0 ) );
    solX                                             /= normFactor;
    solY                                             /= normFactor;
    solZ                                             /= normFactor;
    
    //================================================ Set largest axis element to positive (ProSHADE standard)
    if ( ( ( std::max ( std::abs ( solX ), std::max( std::abs ( solY ), std::abs ( solZ ) ) ) == std::abs ( solX ) ) && ( solX < 0.0 ) ) ||
         ( ( std::max ( std::abs ( solX ), std::max( std::abs ( solY ), std::abs ( solZ ) ) ) == std::abs ( solY ) ) && ( solY < 0.0 ) ) ||
         ( ( std::max ( std::abs ( solX ), std::max( std::abs ( solY ), std::abs ( solZ ) ) ) == std::abs ( solZ ) ) && ( solZ < 0.0 ) ) ) { solX *= -1.0; solY *= -1.0; solZ *= -1.0; }
    
    //================================================ Save solutions
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, solX  );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, solY  );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, solZ  );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function computes matrix multiplication using the ProSHADE group element matrix format as input and output.
 
    \param[in] el1 Group element as rotation matrix in the group element matrix format.
    \param[in] el2 Group element as rotation matrix in the group element matrix format.
    \param[out] ret Matrix in the group element format resulting from the input matrices multiplication.
 */
std::vector< proshade_double > ProSHADE_internal_maths::multiplyGroupElementMatrices ( std::vector< proshade_double >* el1, std::vector< proshade_double >* el2 )
{
    //================================================ Initialise variables
    std::vector< proshade_double > ret;
    
    //================================================ Compute
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, ( el1->at(0) * el2->at(0) ) +
                                                              ( el1->at(1) * el2->at(3) ) +
                                                              ( el1->at(2) * el2->at(6) ) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, ( el1->at(0) * el2->at(1) ) +
                                                              ( el1->at(1) * el2->at(4) ) +
                                                              ( el1->at(2) * el2->at(7) ) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, ( el1->at(0) * el2->at(2) ) +
                                                              ( el1->at(1) * el2->at(5) ) +
                                                              ( el1->at(2) * el2->at(8) ) );
    
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, ( el1->at(3) * el2->at(0) ) +
                                                              ( el1->at(4) * el2->at(3) ) +
                                                              ( el1->at(5) * el2->at(6) ) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, ( el1->at(3) * el2->at(1) ) +
                                                              ( el1->at(4) * el2->at(4) ) +
                                                              ( el1->at(5) * el2->at(7) ) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, ( el1->at(3) * el2->at(2) ) +
                                                              ( el1->at(4) * el2->at(5) ) +
                                                              ( el1->at(5) * el2->at(8) ) );
    
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, ( el1->at(6) * el2->at(0) ) +
                                                              ( el1->at(7) * el2->at(3) ) +
                                                              ( el1->at(8) * el2->at(6) ) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, ( el1->at(6) * el2->at(1) ) +
                                                              ( el1->at(7) * el2->at(4) ) +
                                                              ( el1->at(8) * el2->at(7) ) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, ( el1->at(6) * el2->at(2) ) +
                                                              ( el1->at(7) * el2->at(5) ) +
                                                              ( el1->at(8) * el2->at(8) ) );
    
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function compares the distance between two rotation matrices and decides if they are similar using tolerance.
 
    This function computes the distance between two rotation matrices, specifically by computing the trace of (R1 * R2^T). This measure will be
    3.0 if the two matrices are identical and will decrease the more the rotation matrices difference diverges from identity. Therefore, from this trace
    3.0 is subtracted and the absolute value of the result is compared to the tolerance. If the difference is less than the tolerance, true is returned, while
    false is returned otherwise.
 
    \param[in] mat1 Vector of 9 numbers representing first rotation matrix.
    \param[in] mat1 Vector of 9 numbers representing second rotation matrix.
    \param[in] tolerance Double number representing the maximum allowed error on the distance.
    \param[out] res Boolean decision if the two matrices are similar or not.
 */
bool ProSHADE_internal_maths::rotationMatrixSimilarity ( std::vector< proshade_double >* mat1, std::vector< proshade_double >* mat2, proshade_double tolerance )
{
    //================================================ Initialise variables
    bool ret                                          = false;
    
    //================================================ Compute trace of mat1 * mat2^T
    proshade_double trace                             = ( mat1->at(0) * mat2->at(0) ) + ( mat1->at(1) * mat2->at(1) ) + ( mat1->at(2) * mat2->at(2) );
    trace                                            += ( mat1->at(3) * mat2->at(3) ) + ( mat1->at(4) * mat2->at(4) ) + ( mat1->at(5) * mat2->at(5) );
    trace                                            += ( mat1->at(6) * mat2->at(6) ) + ( mat1->at(7) * mat2->at(7) ) + ( mat1->at(8) * mat2->at(8) );
    
    //================================================ Subtract 3 (so that we would have 0 in case of idenity matrix)
    trace                                            -= 3.0;
    
    //================================================ Compare to tolerance
    if ( tolerance >= std::abs ( trace ) ) { ret = true; }
    
    //================================================ Done
    return                                            ( ret );
    
}
