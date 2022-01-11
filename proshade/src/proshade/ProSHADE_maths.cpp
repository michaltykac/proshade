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
    \version   0.7.6.2
    \date      DEC 2021
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
    
    //================================================ Check for NaN's
    const FloatingPoint< proshade_double > lhs1 ( ret[0] );
    const FloatingPoint< proshade_double > lhs2 ( ret[1] );
    if ( !lhs1.AlmostEquals ( lhs1 ) ) { ret[0] = 0.0; }
    if ( !lhs2.AlmostEquals ( lhs2 ) ) { ret[1] = 0.0; }
    
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
    xMean                                            /= static_cast<proshade_double> ( length ) - zeroCount;
    yMean                                            /= static_cast<proshade_double> ( length ) - zeroCount;
    
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
    if ( std::isnan ( ret ) )                         { return ( 0.0 ); }
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
    proshade_complex* intData                         = new proshade_complex[order];
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
            if ( ( ( static_cast< proshade_double > ( valIt ) * maxSphereDists ) <=  posVals[0] ) && ( ( ( static_cast< proshade_double > ( valIt ) + 1.0 ) * maxSphereDists ) > posVals[0] ) )
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
            lesserWeight                              = static_cast< proshade_double > ( upperPos ) - ( posVals[0] / maxSphereDists );
            upperWeight                               = 1.0 - lesserWeight;
                    
            posVals[1]                                = ( lesserWeight * vals[lesserPos-1] ) + ( upperWeight * vals[upperPos-1] );
        }
        else
        {
            //======================================== ... this then means that we would require position -1 for when the integration value is between 0 and the first shell. To resolve this, we assume that the values are 0 below the first shell and proceed as follows:
            upperWeight                               = 1.0 - ( static_cast< proshade_double > ( upperPos ) - ( posVals[0] / maxSphereDists ) );
                    
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
            if ( ( ( static_cast< proshade_double > ( valIt ) * maxSphereDists ) <=  posVals[0] ) && ( ( ( static_cast< proshade_double > ( valIt ) + 1.0 ) * maxSphereDists ) > posVals[0] ) )
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
            lesserWeight                              = static_cast< proshade_double > ( upperPos ) - ( posVals[0] / maxSphereDists );
            upperWeight                               = 1.0 - lesserWeight;
                    
            posVals[1]                                = ( lesserWeight * vals[lesserPos-1][0] ) + ( upperWeight * vals[upperPos-1][0] );
            posVals[2]                                = ( lesserWeight * vals[lesserPos-1][1] ) + ( upperWeight * vals[upperPos-1][1] );
        }
        else
        {
            //======================================== ... this then means that we would require position -1 for when the integration value is between 0 and the first shell. To resolve this, we assume that the values are 0 below the first shell and proceed as follows:
            upperWeight                               = 1.0 - ( static_cast< proshade_double > ( upperPos ) - ( posVals[0] / maxSphereDists ) );
                    
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
void ProSHADE_internal_maths::realMatrixSVDUandVOnly ( proshade_double* mat, int dim, proshade_double* uAndV, bool fail )
{
    //================================================ Initialise local variables
    char job                                          = 'A';                                   // Save computation of parts of U and V matrices, they are not needed here
    double* singularValues                            = new double[dim];                       // The array of singular values
    double *rotMatU                                   = new double [dim*dim];                  // The U matrix space
    double *rotMatV                                   = new double [dim*dim];                  // The V^T matrix space
    double *work                                      = new double [static_cast< proshade_unsign >( ( 3 * dim ) + pow( dim, 2 ) * dim)]; // Workspace, minimum required is 4*dim^2 + 7*dim, using more for performance
    int workDim                                       = static_cast< int > ( 2 * ( ( 4 * dim * dim ) + ( 7 * dim ) ) ); // Formalism stating just that
    double* rwork                                     = new double[static_cast<proshade_unsign>((5 * dim) + 5 * pow(dim,2))]; // Required by LAPACK
    int* iwork                                        = new int[(8 * dim)];                    // Required by LAPACK
    int returnValue                                   = 0;                                     // This will tell if operation succeeded
    ProSHADE_internal_misc::checkMemoryAllocation     ( singularValues, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMatU,        __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMatV,        __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( work,           __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rwork,          __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( iwork,          __FILE__, __LINE__, __func__ );
    
    //================================================ Load input data into array in column-major order
    double *matrixToDecompose                         = new double[dim*dim];
    ProSHADE_internal_misc::checkMemoryAllocation     ( matrixToDecompose, __FILE__, __LINE__, __func__ );
    for ( int rowIt = 0; rowIt < dim; rowIt++ )
    {
        for ( int colIt = 0; colIt < dim; colIt++ )
        {
            matrixToDecompose[(colIt*dim)+rowIt]      = mat[(rowIt*dim)+colIt];
        }
    }
    
    //================================================ Run LAPACK ZGESDD
    dgesdd_                                           ( &job, &dim, &dim, matrixToDecompose, &dim, singularValues, rotMatU, &dim, rotMatV, &dim,
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
            uAndV[(rowIt*3)+colIt]                    = rotMatU[( rowIt * 3 ) + colIt];
        }
    }
    
    //================================================ Save V
    for ( proshade_signed rowIt = 0; rowIt < dim; rowIt++ )
    {
        for ( proshade_signed colIt = 0; colIt < dim; colIt++ )
        {
            uAndV[(rowIt*3)+colIt+9]                  = rotMatV[( rowIt * 3 ) + colIt];
        }
    }
    
    //================================================ Release the rest of the memory
    delete[] rotMatU;
    delete[] rotMatV;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function to find Euler angles (ZYZ convention) from index position in the inverse SOFT map.
 
    This function proceeds to convert the inverse SOFT map x, y and z position to Euler ZYZ convention angles, saving
    these into the supplied pointers.
 
    \param[in] band The maximum bandwidth of the computation.
    \param[in] x The x-axis position in the inverse SOFT map.
    \param[in] y The y-axis position in the inverse SOFT map.
    \param[in] z The z-axis position in the inverse SOFT map.
    \param[in] eulerAlpha Pointer to where the Euler alpha angle will be saved.
    \param[in] eulerBeta Pointer to where the Euler beta angle will be saved.
    \param[in] eulerGamma Pointer to where the Euler gamma angle will be saved.
 */
void ProSHADE_internal_maths::getEulerZYZFromSOFTPosition ( proshade_signed band, proshade_signed x, proshade_signed y, proshade_signed z, proshade_double* eulerAlpha, proshade_double* eulerBeta, proshade_double* eulerGamma )
{
    //================================================ Convert index to Euler angles
   *eulerAlpha                                        = M_PI *         static_cast<proshade_double> ( y )    / (       static_cast<proshade_double> ( band ) ) ;
   *eulerBeta                                         = M_PI * ( 2.0 * static_cast<proshade_double> ( x )  ) / ( 4.0 * static_cast<proshade_double> ( band ) ) ;
   *eulerGamma                                        = M_PI *         static_cast<proshade_double> ( z )    / (       static_cast<proshade_double> ( band ) ) ;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function to find the index position in the inverse SOFT map from given Euler angles (ZYZ convention).
 
    This function does the conversion from Euler angles ZYZ convention to the SOFT map x, y and z position. It is not limitted to
    the SOFT map indices and instead if given Euler agnles between two indices will return a decimal point for the indices.
 
    \param[in] band The maximum bandwidth of the computation.
    \param[in] eulerAlpha The Euler alpha angle value.
    \param[in] eulerBeta The Euler beta angle value.
    \param[in] eulerGamma The Euler gamma angle value.
    \param[in] x Pointer to where the closest x-axis position in the inverse SOFT map will be saved to (position may be decimal!).
    \param[in] y Pointer to where the closest y-axis position in the inverse SOFT map will be saved to (position may be decimal!).
    \param[in] z Pointer to where the closest z-axis position in the inverse SOFT map will be saved to (position may be decimal!).
 */
void ProSHADE_internal_maths::getSOFTPositionFromEulerZYZ ( proshade_signed band, proshade_double eulerAlpha, proshade_double eulerBeta, proshade_double eulerGamma, proshade_double* x, proshade_double* y, proshade_double* z )
{
    //================================================ Convert Euler angles to indices
   *x                                                 = ( eulerBeta  * static_cast<proshade_double> ( band ) * 2.0 ) / M_PI;
   *y                                                 = ( eulerGamma * static_cast<proshade_double> ( band )       ) / M_PI;
   *z                                                 = ( eulerAlpha * static_cast<proshade_double> ( band )       ) / M_PI;
    
    //================================================ Deal with singularities
    if ( eulerBeta > ( M_PI - 0.05 ) )
    {
        //============================================ Rotation is 180 deg
       *z                                             = ( ( eulerAlpha - eulerGamma ) * static_cast<proshade_double> ( band ) ) / M_PI;
       *y                                             = 0;
    }
    
    //================================================ Keep value within boundaries, but do not repeat over them!
    if ( *x >= ( 2 * band ) ) { *x = ( 2 * band ) - 1; }
    if ( *y >= ( 2 * band ) ) { *y = ( 2 * band ) - 1; }
    if ( *z >= ( 2 * band ) ) { *z = ( 2 * band ) - 1; }
    
    //================================================ Done
    return ;
    
}

/*! \brief Function to find the rotation matrix from Euler angles (ZYZ convention).
 
    \param[in] eulerAlpha The Euler alpha angle value.
    \param[in] eulerBeta The Euler beta angle value.
    \param[in] eulerGamma The Euler gamma angle value.
    \param[in] matrix A pointer to array of 9 values to which the results of the function will be saved.
 */
void ProSHADE_internal_maths::getRotationMatrixFromEulerZYZAngles ( proshade_double eulerAlpha, proshade_double eulerBeta, proshade_double eulerGamma, proshade_double* matrix )
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

/*! \brief Function to find the rotation matrix from Euler angles (ZYZ convention).
 
    \param[in] eulerAlpha The Euler alpha angle value.
    \param[in] eulerBeta The Euler beta angle value.
    \param[in] eulerGamma The Euler gamma angle value.
    \param[in] matrix A pointer to array of 9 values to which the results of the function will be saved.
 */
void ProSHADE_internal_maths::getRotationMatrixFromEulerZYZAngles ( proshade_single eulerAlpha, proshade_single eulerBeta, proshade_single eulerGamma, proshade_single* matrix )
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
    \param[in] verbose Should the warnings be printed? -1 if not.
 */
 void ProSHADE_internal_maths::getAxisAngleFromRotationMatrix ( proshade_double* rotMat, proshade_double* x, proshade_double* y, proshade_double* z, proshade_double* ang, proshade_signed verbose )
{
    //================================================ Initialise
    proshade_double angleTolerance                    = 0.01;
    proshade_double closeToZero                       = 0.0000001;
    
    //================================================ Find the angle
   *ang                                               = std::acos ( ( std::max ( -1.0, std::min ( 3.0, rotMat[0] + rotMat[4] + rotMat[8] ) ) - 1.0 ) / 2.0 );
    
    //================================================ Any singularity?
    if ( std::abs ( std::sin ( *ang ) ) < angleTolerance )
    {
        //============================================ Initialise local variables
        char jobLeftEigs                              = 'N';                                   // This tells LAPACK not to compute left eigenvectors.
        char jobRightEigs                             = 'V';                                   // This tells LAPACK to compute right eigenvectors.
        int dim                                       = 3;                                     // The order of the matrix
        double* eigValReal                            = new double[dim];                       // Real parts of the eigenvalues
        double* eigValImag                            = new double[dim];                       // Imaginary parts of the eigenvalues
        double* leftEigVectors                        = new double[dim*dim*2];                 // Left eigenvectors containing matrix
        double* rightEigVectors                       = new double[dim*dim*2];                 // Right eigenvectors containing matrix
        double* work                                  = new double[10*4*dim];                  // Workspace. 4 * dim should be minimum, but more is better
        int workSize                                  = 10*4*dim;                              // Saving the work array size for passing.
        int returnValue                               = 0;                                     // This will tell if operation succeeded
        
        //============================================ Check memory allocation
        ProSHADE_internal_misc::checkMemoryAllocation ( eigValReal,      __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( eigValImag,      __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( leftEigVectors,  __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( rightEigVectors, __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( work,            __FILE__, __LINE__, __func__ );
        
        //============================================ Load input data into array in column-major order
        double* matrixToDecompose                     = new double[dim*dim];
        ProSHADE_internal_misc::checkMemoryAllocation ( matrixToDecompose, __FILE__, __LINE__, __func__ );
        for ( int rowIt = 0; rowIt < dim; rowIt++ )
        {
            for ( int colIt = 0; colIt < dim; colIt++ )
            {
                matrixToDecompose[(colIt*dim)+rowIt]  = static_cast< double > ( rotMat[(rowIt*dim)+colIt] );
            }
        }
        
        //============================================ Run LAPACK ZGESDD
        dgeev_                                        ( &jobLeftEigs, &jobRightEigs, &dim, matrixToDecompose, &dim, eigValReal, eigValImag, leftEigVectors, &dim,
                                                        rightEigVectors, &dim, work, &workSize, &returnValue );
        
        //============================================ Check for errors
        if ( returnValue != 0 )
        {
            //======================================== Report error and return zero values
            ProSHADE_internal_messages::printWarningMessage ( verbose, "!!! ProSHADE WARNING !!! Eigenval/Eigenvector algorithm did not converge.", "WS00069" );
           *x                                         = 0.0;
           *y                                         = 0.0;
           *z                                         = 0.0;
            
            //======================================== Release memory
            delete[] eigValReal;
            delete[] eigValImag;
            delete[] leftEigVectors;
            delete[] rightEigVectors;
            delete[] work;
            delete[] matrixToDecompose;
            
            //======================================== Done
            return ;
        }
        
        //============================================ If values are close to zero, just set them to zero
        for ( int i = 0; i < 9; i++ ) { if ( std::abs(rightEigVectors[i]) < closeToZero ) { rightEigVectors[i] = 0.0; } }
        for ( int i = 0; i < 3; i++ ) { if ( std::abs(eigValReal[i]) < closeToZero ) { eigValReal[i] = 0.0; } if ( std::abs(eigValImag[i]) < closeToZero ) { eigValImag[i] = 0.0; } }
        
        //============================================ Find which eigenvalue is close to 1 and has imaginary part close to zero
        proshade_signed eigIt                         = -1;
        for ( size_t it = 0; it < 3; it++ )
        {
            if ( ( eigValReal[it] > ( 1.0 - closeToZero ) ) && ( eigValReal[it] < ( 1.0 + closeToZero ) ) )
            {
                if ( ( eigValImag[it] > ( 0.0 - closeToZero ) ) && ( eigValImag[it] < ( 0.0 + closeToZero ) ) )
                {
                    eigIt                             = static_cast< proshade_signed > ( it );
                    break;
                }
            }
        }
        
        //============================================ Any axis found?
        if ( eigIt == -1 )
        {
            ProSHADE_internal_messages::printWarningMessage ( verbose, "!!! ProSHADE WARNING !!! Failed to find eigenvalue with value 1 for this rotation matrix. Is this a rotation matrix?", "WS00072" );
           *x                                         = 0.0;
           *y                                         = 0.0;
           *z                                         = 0.0;
        }
        else
        {
            //======================================== Parse LAPACK eigenvectors matrix
            int colIt;
            for( int rowIt = 0; rowIt < dim; rowIt++ )
            {
                colIt                                 = 0;
                while( colIt < dim )
                {
                    if( std::abs ( eigValImag[colIt] ) < closeToZero )
                    {
                        if ( colIt == eigIt ) { if ( rowIt == 0 ) { *x = rightEigVectors[rowIt+colIt*dim]; } if ( rowIt == 1 ) { *y = rightEigVectors[rowIt+colIt*dim]; } if ( rowIt == 2 ) { *z = rightEigVectors[rowIt+colIt*dim]; } }
                        colIt++;
                    }
                    else
                    {
// In order to access the
//                        std::cout << " ( " << rightEigVectors[rowIt+colIt*dim] << " + " <<  rightEigVectors[rowIt+(colIt+1)*dim] << " i )";
// other eigenvectors, use this:
//                        std::cout << " ( " << rightEigVectors[rowIt+colIt*dim] << " + " << -rightEigVectors[rowIt+(colIt+1)*dim] << " i )";
                        colIt                        += 2;
                    }
                }
            }
        }
        
        //============================================ Normalise axis length
        proshade_double normFactor                    = std::sqrt ( pow ( *x, 2.0 ) + pow ( *y, 2.0 ) + pow ( *z, 2.0 ) );
       *x                                            /= normFactor;
       *y                                            /= normFactor;
       *z                                            /= normFactor;
        

        //============================================ Free memory
        delete[] eigValReal;
        delete[] eigValImag;
        delete[] leftEigVectors;
        delete[] rightEigVectors;
        delete[] work;
        delete[] matrixToDecompose;
    }
    else
    {
       //============================================= Axis
      *x                                              = rotMat[7] - rotMat[5];
      *y                                              = rotMat[2] - rotMat[6];
      *z                                              = rotMat[3] - rotMat[1];
       
       proshade_double normFactor                     = std::sqrt ( pow ( *x, 2.0 ) + pow ( *y, 2.0 ) + pow ( *z, 2.0 ) );
      *x                                             /= normFactor;
      *y                                             /= normFactor;
      *z                                             /= normFactor;
       
       //============================================= Make sure largest axis is positive
       const FloatingPoint< proshade_double > lhs1 ( std::max ( std::abs ( *x ), std::max ( std::abs ( *y ), std::abs ( *z ) ) ) );
       const FloatingPoint< proshade_double > rhs1 ( std::abs ( *x ) );
       const FloatingPoint< proshade_double > rhs2 ( std::abs ( *y ) );
       const FloatingPoint< proshade_double > rhs3 ( std::abs ( *z ) );
       if ( ( ( lhs1.AlmostEquals ( rhs1 ) ) && ( *x < 0.0 ) ) ||
            ( ( lhs1.AlmostEquals ( rhs2 ) ) && ( *y < 0.0 ) ) ||
            ( ( lhs1.AlmostEquals ( rhs3 ) ) && ( *z < 0.0 ) ) )
       {
           *x                                        *= -1.0;
           *y                                        *= -1.0;
           *z                                        *= -1.0;
           *ang                                      *= -1.0;
       }
    }
    
    //================================================ Standardise angle to range 0 to 2pi
    if ( *ang < 0.0 ) { *ang = ( 2.0 * M_PI ) + *ang; }

    //================================================ Done
    return ;
    
}

/*! \brief This function converts rotation matrix to the axis-angle representation.
 
    This function takes a rotation matrix as a pointer to a vector of doubles and converts it to the Angle-Axis
    representation, which is the main rotation representation used in ProSHADE. This function deals with both
    the North and South pole singularity of the rotation matrices.
 
    \param[in] rotMat Rotation matrix as a  pointer to a vector of doubles.
    \param[in] x Pointer to which the x-axis value of the axis vector will be saved.
    \param[in] y Pointer to which the y-axis value of the axis vector will be saved.
    \param[in] z Pointer to which the z-axis value of the axis vector will be saved.
    \param[in] ang Pointer to which the angle value will be saved.
    \param[in] verbose Should the warnings be printed? -1 if not.
 */
 void ProSHADE_internal_maths::getAxisAngleFromRotationMatrix ( std::vector< proshade_double >* rotMat, proshade_double* x, proshade_double* y, proshade_double* z, proshade_double* ang, proshade_signed verbose )
{
    //================================================ Initialise
    proshade_double angleTolerance                    = 0.01;
    proshade_double closeToZero                       = 0.0000001;
    
    //================================================ Find the angle
   *ang                                               = std::acos ( ( std::max ( -1.0, std::min ( 3.0, rotMat->at(0) + rotMat->at(4) + rotMat->at(8) ) ) - 1.0 ) / 2.0 );
    
    //================================================ Any singularity?
    if ( std::abs ( std::sin ( *ang ) ) < angleTolerance )
    {
        //============================================ Initialise local variables
        char jobLeftEigs                              = 'N';                                   // This tells LAPACK not to compute left eigenvectors.
        char jobRightEigs                             = 'V';                                   // This tells LAPACK to compute right eigenvectors.
        int dim                                       = 3;                                     // The order of the matrix
        double* eigValReal                            = new double[dim];                       // Real parts of the eigenvalues
        double* eigValImag                            = new double[dim];                       // Imaginary parts of the eigenvalues
        double* leftEigVectors                        = new double[dim*dim*2];                 // Left eigenvectors containing matrix
        double* rightEigVectors                       = new double[dim*dim*2];                 // Right eigenvectors containing matrix
        double* work                                  = new double[10*4*dim];                  // Workspace. 4 * dim should be minimum, but more is better
        int workSize                                  = 10*4*dim;                              // Saving the work array size for passing.
        int returnValue                               = 0;                                     // This will tell if operation succeeded
        
        //============================================ Check memory allocation
        ProSHADE_internal_misc::checkMemoryAllocation ( eigValReal,      __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( eigValImag,      __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( leftEigVectors,  __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( rightEigVectors, __FILE__, __LINE__, __func__ );
        ProSHADE_internal_misc::checkMemoryAllocation ( work,            __FILE__, __LINE__, __func__ );
        
        //============================================ Load input data into array in column-major order
        double* matrixToDecompose                     = new double[dim*dim];
        ProSHADE_internal_misc::checkMemoryAllocation ( matrixToDecompose, __FILE__, __LINE__, __func__ );
        for ( int rowIt = 0; rowIt < dim; rowIt++ )
        {
            for ( int colIt = 0; colIt < dim; colIt++ )
            {
                matrixToDecompose[(colIt*dim)+rowIt]  = static_cast< double > ( rotMat->at( static_cast< size_t > ( ( rowIt * dim ) + colIt ) ) );
            }
        }
        
        //============================================ Run LAPACK ZGESDD
        dgeev_                                        ( &jobLeftEigs, &jobRightEigs, &dim, matrixToDecompose, &dim, eigValReal, eigValImag, leftEigVectors, &dim,
                                                        rightEigVectors, &dim, work, &workSize, &returnValue );
        
        //============================================ Check for errors
        if ( returnValue != 0 )
        {
            //======================================== Report error and return zero values
            ProSHADE_internal_messages::printWarningMessage ( verbose, "!!! ProSHADE WARNING !!! Eigenval/Eigenvector algorithm did not converge.", "WS00069" );
           *x                                         = 0.0;
           *y                                         = 0.0;
           *z                                         = 0.0;
            
            //======================================== Release memory
            delete[] eigValReal;
            delete[] eigValImag;
            delete[] leftEigVectors;
            delete[] rightEigVectors;
            delete[] work;
            delete[] matrixToDecompose;
            
            //======================================== Done
            return ;
        }
        
        //============================================ If values are close to zero, just set them to zero
        for ( int i = 0; i < 9; i++ ) { if ( std::abs(rightEigVectors[i]) < closeToZero ) { rightEigVectors[i] = 0.0; } }
        for ( int i = 0; i < 3; i++ ) { if ( std::abs(eigValReal[i]) < closeToZero ) { eigValReal[i] = 0.0; } if ( std::abs(eigValImag[i]) < closeToZero ) { eigValImag[i] = 0.0; } }
        
        //============================================ Find which eigenvalue is close to 1 and has imaginary part close to zero
        proshade_signed eigIt                         = -1;
        for ( size_t it = 0; it < 3; it++ )
        {
            if ( ( eigValReal[it] > ( 1.0 - closeToZero ) ) && ( eigValReal[it] < ( 1.0 + closeToZero ) ) )
            {
                if ( ( eigValImag[it] > ( 0.0 - closeToZero ) ) && ( eigValImag[it] < ( 0.0 + closeToZero ) ) )
                {
                    eigIt                             = static_cast< proshade_signed > ( it );
                    break;
                }
            }
        }
        
        //============================================ Any axis found?
        if ( eigIt == -1 )
        {
            ProSHADE_internal_messages::printWarningMessage ( verbose, "!!! ProSHADE WARNING !!! Failed to find eigenvalue with value 1 for this rotation matrix. Is this a rotation matrix?", "WS00072" );
           *x                                         = 0.0;
           *y                                         = 0.0;
           *z                                         = 0.0;
        }
        else
        {
            //======================================== Parse LAPACK eigenvectors matrix
            int colIt;
            for( int rowIt = 0; rowIt < dim; rowIt++ )
            {
                colIt                                 = 0;
                while( colIt < dim )
                {
                    if( std::abs ( eigValImag[colIt] ) < closeToZero )
                    {
                        if ( colIt == eigIt ) { if ( rowIt == 0 ) { *x = rightEigVectors[rowIt+colIt*dim]; } if ( rowIt == 1 ) { *y = rightEigVectors[rowIt+colIt*dim]; } if ( rowIt == 2 ) { *z = rightEigVectors[rowIt+colIt*dim]; } }
                        colIt++;
                    }
                    else
                    {
// In order to access the
//                        std::cout << " ( " << rightEigVectors[rowIt+colIt*dim] << " + " <<  rightEigVectors[rowIt+(colIt+1)*dim] << " i )";
// other eigenvectors, use this:
//                        std::cout << " ( " << rightEigVectors[rowIt+colIt*dim] << " + " << -rightEigVectors[rowIt+(colIt+1)*dim] << " i )";
                        colIt                        += 2;
                    }
                }
            }
        }
        
        //============================================ Normalise axis length
        proshade_double normFactor                    = std::sqrt ( pow ( *x, 2.0 ) + pow ( *y, 2.0 ) + pow ( *z, 2.0 ) );
       *x                                            /= normFactor;
       *y                                            /= normFactor;
       *z                                            /= normFactor;
        
        //============================================= Make sure largest axis is positive
        const FloatingPoint< proshade_double > lhs1 ( std::max ( std::abs ( *x ), std::max ( std::abs ( *y ), std::abs ( *z ) ) ) );
        const FloatingPoint< proshade_double > rhs1 ( std::abs ( *x ) );
        const FloatingPoint< proshade_double > rhs2 ( std::abs ( *y ) );
        const FloatingPoint< proshade_double > rhs3 ( std::abs ( *z ) );
        if ( ( ( lhs1.AlmostEquals ( rhs1 ) ) && ( *x < 0.0 ) ) ||
             ( ( lhs1.AlmostEquals ( rhs2 ) ) && ( *y < 0.0 ) ) ||
             ( ( lhs1.AlmostEquals ( rhs3 ) ) && ( *z < 0.0 ) ) )
        {
            *x                                       *= -1.0;
            *y                                       *= -1.0;
            *z                                       *= -1.0;
            *ang                                     *= -1.0;
        }

        //================================================ Free memory
        delete[] eigValReal;
        delete[] eigValImag;
        delete[] leftEigVectors;
        delete[] rightEigVectors;
        delete[] work;
        delete[] matrixToDecompose;
    }
    else
    {
       //============================================= Axis
      *x                                              = rotMat->at(7) - rotMat->at(5);
      *y                                              = rotMat->at(2) - rotMat->at(6);
      *z                                              = rotMat->at(3) - rotMat->at(1);
       
       proshade_double normFactor                     = std::sqrt ( pow ( *x, 2.0 ) + pow ( *y, 2.0 ) + pow ( *z, 2.0 ) );
      *x                                             /= normFactor;
      *y                                             /= normFactor;
      *z                                             /= normFactor;
       
       //============================================= Make sure largest axis is positive
       const FloatingPoint< proshade_double > lhs1 ( std::max ( std::abs ( *x ), std::max ( std::abs ( *y ), std::abs ( *z ) ) ) );
       const FloatingPoint< proshade_double > rhs1 ( std::abs ( *x ) );
       const FloatingPoint< proshade_double > rhs2 ( std::abs ( *y ) );
       const FloatingPoint< proshade_double > rhs3 ( std::abs ( *z ) );
       if ( ( ( lhs1.AlmostEquals ( rhs1 ) ) && ( *x < 0.0 ) ) ||
            ( ( lhs1.AlmostEquals ( rhs2 ) ) && ( *y < 0.0 ) ) ||
            ( ( lhs1.AlmostEquals ( rhs3 ) ) && ( *z < 0.0 ) ) )
       {
           *x                                        *= -1.0;
           *y                                        *= -1.0;
           *z                                        *= -1.0;
           *ang                                      *= -1.0;
       }
    }
    
    //================================================ Standardise angle to range 0 to 2pi
    if ( *ang < 0.0 ) { *ang = ( 2.0 * M_PI ) + *ang; }

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

/*! \brief This function converts the axis-angle representation to the rotation matrix representation.
 
    \param[in] rotMat Rotation matrix as an array of 9 values will be saved to this pointer, must already be allocated.
    \param[in] x The x-axis value of the axis vector.
    \param[in] y The y-axis value of the axis vector.
    \param[in] z The z-axis value of the axis vector.
    \param[in] angThe angle value.
 */
void ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( proshade_single* rotMat, proshade_double x, proshade_double y, proshade_double z, proshade_double ang )
{
    //================================================ If angle is 0 or infinity (anything divided by 0), return identity matrix
    if ( ( ang == 0.0 ) || ( std::isinf ( ang ) ) )
    {
        //============================================ Create identity
        for ( size_t i = 0; i < 9; i++ ) { rotMat[i] = 0.0f; }
        rotMat[0]                                     = 1.0f;
        rotMat[4]                                     = 1.0f;
        rotMat[8]                                     = 1.0f;
        
        //============================================ Done
        return ;
    }
    
    //================================================ Compute the matrix
    proshade_single cAng                              = cos ( static_cast< proshade_single > ( ang ) );
    proshade_single sAng                              = sin ( static_cast< proshade_single > ( ang ) );
    proshade_single tAng                              = 1.0f - cAng;
            
    rotMat[0]                                         = cAng + static_cast< proshade_single > ( x ) * static_cast< proshade_single > ( x ) * tAng;
    rotMat[4]                                         = cAng + static_cast< proshade_single > ( y ) * static_cast< proshade_single > ( y ) * tAng;
    rotMat[8]                                         = cAng + static_cast< proshade_single > ( z ) * static_cast< proshade_single > ( z ) * tAng;
            
    proshade_single tmp1                              = static_cast< proshade_single > ( x ) * static_cast< proshade_single > ( y ) * tAng;
    proshade_single tmp2                              = static_cast< proshade_single > ( z ) * sAng;
    rotMat[3]                                         = tmp1 + tmp2;
    rotMat[1]                                         = tmp1 - tmp2;
            
    tmp1                                              = static_cast< proshade_single > ( x ) * static_cast< proshade_single > ( z ) * tAng;
    tmp2                                              = static_cast< proshade_single > ( y ) * sAng;
    rotMat[6]                                         = tmp1 - tmp2;
    rotMat[2]                                         = tmp1 + tmp2;
            
    tmp1                                              = static_cast< proshade_single > ( y ) * static_cast< proshade_single > ( z ) * tAng;
    tmp2                                              = static_cast< proshade_single > ( x ) * sAng;
    rotMat[7]                                         = tmp1 + tmp2;
    rotMat[5]                                         = tmp1 - tmp2;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function converts rotation matrix to the Euler ZYZ angles representation.
 
    \param[in] rotMat Rotation matrix as an array of 9 values.
    \param[in] eA Pointer to which the Euler angle alpha value will be saved.
    \param[in] eB Pointer to which the Euler angle beta value will be saved.
    \param[in] eG Pointer to which the Euler angle gamma value will be saved.
 */
void ProSHADE_internal_maths::getEulerZYZFromRotMatrix ( proshade_double* rotMat, proshade_double* eA, proshade_double* eB, proshade_double* eG )
{
    //================================================ Convert to Eulers
    if ( std::abs( rotMat[8] ) < 0.99999 )
    {
        //============================================ This case occurs when there is no singularity in the rotation matrix (i.e. it does not have 0 or 180 degrees angle)
       *eA                                            = std::atan2 ( rotMat[7],  rotMat[6] );
       *eB                                            = std::acos  ( rotMat[8] );
       *eG                                            = std::atan2 ( rotMat[5], -rotMat[2] );
    }
    else
    {
        //============================================ This case occurs when there is either 0 or 180 degrees rotation angle in the rotation matrix and therefore when beta is zero.
        if ( rotMat[8] >= 0.99999 )
        {
            //======================================== In this case, beta = 0 and alpha and gamma are only defined in terms of their sum. So we arbitrarily set gamma to 0 and solve alpha.
           *eA                                        = std::atan2 ( rotMat[3], rotMat[0] );
           *eB                                        = 0.0;
           *eG                                        = 0.0;
        }
        if ( rotMat[8] <= -0.99999 )
        {
            //======================================== In this case, beta = PI and alpha and gamma are only defined in terms of their difference. So we arbitrarily set gamma to 0 and solve alpha.
           *eA                                        = std::atan2 ( rotMat[3], rotMat[0] );
           *eB                                        = M_PI;
           *eG                                        = 0.0;
        }
    }
    
    //================================================ Get the angles to proper range
    if ( *eA < 0.0 ) { *eA                            = 2.0 * M_PI + *eA; }
    if ( *eB < 0.0 ) { *eB                            =       M_PI + *eB; }
    if ( *eG < 0.0 ) { *eG                            = 2.0 * M_PI + *eG; }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function converts angle-axis representation to the Euler ZYZ angles representation.
 
    This function does the angle-axis to Euler ZYZ conversion and if a problem around the Z axis arises, it deal with it.
 
    \param[in] axX Angle-axis representation axis x element.
    \param[in] axY Angle-axis representation axis y element.
    \param[in] axZ Angle-axis representation axis z element.
    \param[in] axAng Angle-axis representation angle.
    \param[in] eA Pointer to which the Euler angle alpha value will be saved.
    \param[in] eB Pointer to which the Euler angle beta value will be saved.
    \param[in] eG Pointer to which the Euler angle gamma value will be saved.
 */
void ProSHADE_internal_maths::getEulerZYZFromAngleAxis ( proshade_double axX, proshade_double axY, proshade_double axZ, proshade_double axAng, proshade_double* eA, proshade_double* eB, proshade_double* eG )
{
    //================================================ If angle is 0 or infinity (anything divided by 0), return no rotation
    if ( ( axAng == 0.0 ) || ( std::isinf ( axAng ) ) )
    {
        //============================================ Return 0 ; 0 ; 0 for no angle
       *eA                                            = 0.0;
       *eB                                            = 0.0;
       *eG                                            = 0.0;
        
        //============================================ Done
        return ;
    }
    
    //================================================ Compute required rotation matrix elements
    proshade_double cAng                              = std::cos ( axAng );
    proshade_double sAng                              = std::sin ( axAng );
    proshade_double tAng                              = 1.0 - cAng;
    
    proshade_double element22                         = cAng + axZ * axZ * tAng;
            
    proshade_double tmp1                              = axX * axZ * tAng;
    proshade_double tmp2                              = axY * sAng;
    proshade_double element20                         = tmp1 - tmp2;
    proshade_double element02                         = tmp1 + tmp2;
            
    tmp1                                              = axY * axZ * tAng;
    tmp2                                              = axX * sAng;
    proshade_double element21                         = tmp1 + tmp2;
    proshade_double element12                         = tmp1 - tmp2;
    
    //================================================ Convert to Eulers
    if ( std::abs( element22 ) <= 0.99999 )
    {
        //============================================ This case occurs when there is no singularity in the rotation matrix (i.e. it does not have 0 or 180 degrees angle)
       *eA                                            = std::atan2 ( element21,  element20 );
       *eB                                            = std::acos  ( element22 );
       *eG                                            = std::atan2 ( element12, -element02 );
    }
    else
    {
        //============================================ Compute some extra rotation matrix elements
        tmp1                                          = axX * axY * tAng;
        tmp2                                          = axZ * sAng;
        proshade_double element10                     = tmp1 + tmp2;
        proshade_double element00                     = cAng + axX * axX * tAng;
        
        //============================================ This case occurs when there is either 0 or 180 degrees rotation angle in the rotation matrix and therefore when beta is zero.
        if ( element22 >= 0.99999 )
        {
            //======================================== In this case, beta = 0 and alpha and gamma are only defined in terms of their sum. So we arbitrarily set gamma to 0 and solve alpha.
           *eA                                        = std::atan2 ( element10, element00 );
           *eB                                        = 0.0;
           *eG                                        = 0.0;
        }
        if ( element22 <= -0.99999 )
        {
            //======================================== In this case, beta = PI and alpha and gamma are only defined in terms of their difference. So we arbitrarily set gamma to 0 and solve alpha.
           *eA                                        = std::atan2 ( element10, element00 );
           *eB                                        = M_PI;
           *eG                                        = 0.0;
        }
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
    //================================================ Set res to 0.0s
    for ( proshade_unsign iter = 0; iter < 9; iter++ ) { res[iter] = 0.0; }
    
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
    for ( proshade_double posDiv = 3; posDiv <= sqrt ( static_cast< proshade_double > ( number ) ); posDiv += 2.0 )
    {
        // If posDiv is a divisor of the number, save the result
        while ( number % static_cast< proshade_signed > ( posDiv ) == 0 )
        {
            ProSHADE_internal_misc::addToSignedVector ( &ret, static_cast< proshade_signed > ( posDiv ) );
            number                                    = number / static_cast< proshade_signed > ( posDiv );
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
    \param[out] X The dot product of the two input vectors.
 */
proshade_double ProSHADE_internal_maths::computeDotProduct ( proshade_double* x1, proshade_double* y1, proshade_double* z1, proshade_double* x2, proshade_double* y2, proshade_double* z2 )
{
    //================================================ Compute and return
    return                                            ( (*x1 * *x2) + (*y1 * *y2) + (*z1 * *z2) );
}

/*! \brief Simple 3D vector dot product computation.
 
    \param[in] x1 The x-axis element of the first vector.
    \param[in] y1 The y-axis element of the first vector.
    \param[in] z1 The z-axis element of the first vector.
    \param[in] x2 The x-axis element of the second vector.
    \param[in] y2 The y-axis element of the second vector.
    \param[in] z2 The z-axis element of the second vector.
    \param[out] X The dot product of the two input vectors.
 */
proshade_double ProSHADE_internal_maths::computeDotProduct ( proshade_double x1, proshade_double y1, proshade_double z1, proshade_double x2, proshade_double y2, proshade_double z2 )
{
    //================================================ Compute and return
    return                                            ( (x1 * x2) + (y1 * y2) + (z1 * z2) );
}

/*! \brief Simple 3D vector cross product computation.
 
    \param[in] x1 The x-axis element of the first vector.
    \param[in] y1 The y-axis element of the first vector.
    \param[in] z1 The z-axis element of the first vector.
    \param[in] x2 The x-axis element of the second vector.
    \param[in] y2 The y-axis element of the second vector.
    \param[in] z2 The z-axis element of the second vector.
    \param[out] crossProd The vector representing the cross product of the two input vectors.
 */
proshade_double* ProSHADE_internal_maths::computeCrossProduct ( proshade_double* x1, proshade_double* y1, proshade_double* z1, proshade_double* x2, proshade_double* y2, proshade_double* z2 )
{
    //================================================ Allocate memory
    proshade_double* crossProd                        = new proshade_double[3];
    ProSHADE_internal_misc::checkMemoryAllocation     ( crossProd, __FILE__, __LINE__, __func__ );
    
    //================================================ Compute
    crossProd[0]                                      = ( (*y1) * (*z2) ) - ( (*z1) * (*y2) );
    crossProd[1]                                      = ( (*z1) * (*x2) ) - ( (*x1) * (*z2) );
    crossProd[2]                                      = ( (*x1) * (*y2) ) - ( (*y1) * (*x2) );
    
    //================================================ Done
    return                                            ( crossProd );
    
}

/*! \brief Simple 3D vector cross product computation.
 
    \param[in] x1 The x-axis element of the first vector.
    \param[in] y1 The y-axis element of the first vector.
    \param[in] z1 The z-axis element of the first vector.
    \param[in] x2 The x-axis element of the second vector.
    \param[in] y2 The y-axis element of the second vector.
    \param[in] z2 The z-axis element of the second vector.
    \param[out] crossProd The vector representing the cross product of the two input vectors.
 */
proshade_double* ProSHADE_internal_maths::computeCrossProduct ( proshade_double x1, proshade_double y1, proshade_double z1, proshade_double x2, proshade_double y2, proshade_double z2 )
{
    //================================================ Allocate memory
    proshade_double* crossProd                        = new proshade_double[3];
    ProSHADE_internal_misc::checkMemoryAllocation     ( crossProd, __FILE__, __LINE__, __func__ );
    
    //================================================ Compute
    crossProd[0]                                      = ( y1 * z2 ) - ( z1 * y2 );
    crossProd[1]                                      = ( z1 * x2 ) - ( x1 * z2 );
    crossProd[2]                                      = ( x1 * y2 ) - ( y1 * x2 );
    
    //================================================ Done
    return                                            ( crossProd );
    
}

/*! \brief Function for computing a 3x3 matrix multiplication.
 
    \param[in] mat1 The matrix to multiply mat2.
    \param[in] mat2 The matrix to be multiplied by mat1.
    \param[out] ret The matrix resulting from matrix multiplication of mat1 and mat2 in this order.
 */
proshade_double* ProSHADE_internal_maths::compute3x3MatrixMultiplication ( proshade_double* mat1, proshade_double* mat2 )
{
    //================================================ Allocate memory
    proshade_double* ret                              = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( ret, __FILE__, __LINE__, __func__ );
    
    //================================================ Multiply
    ret[0]                                            = ( mat1[0] * mat2[0] ) + ( mat1[1] * mat2[3] ) + ( mat1[2] * mat2[6] );
    ret[1]                                            = ( mat1[0] * mat2[1] ) + ( mat1[1] * mat2[4] ) + ( mat1[2] * mat2[7] );
    ret[2]                                            = ( mat1[0] * mat2[2] ) + ( mat1[1] * mat2[5] ) + ( mat1[2] * mat2[8] );
    ret[3]                                            = ( mat1[3] * mat2[0] ) + ( mat1[4] * mat2[3] ) + ( mat1[5] * mat2[6] );
    ret[4]                                            = ( mat1[3] * mat2[1] ) + ( mat1[4] * mat2[4] ) + ( mat1[5] * mat2[7] );
    ret[5]                                            = ( mat1[3] * mat2[2] ) + ( mat1[4] * mat2[5] ) + ( mat1[5] * mat2[8] );
    ret[6]                                            = ( mat1[6] * mat2[0] ) + ( mat1[7] * mat2[3] ) + ( mat1[8] * mat2[6] );
    ret[7]                                            = ( mat1[6] * mat2[1] ) + ( mat1[7] * mat2[4] ) + ( mat1[8] * mat2[7] );
    ret[8]                                            = ( mat1[6] * mat2[2] ) + ( mat1[7] * mat2[5] ) + ( mat1[8] * mat2[8] );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief Function for computing a 3x3 matrix to 3x1 vector multiplication.
 
    \param[in] mat The matrix to multiply the vector with..
    \param[in] x The x-axis element of the vector which is to be multiplied by the matrix.
    \param[in] y The x-axis element of the vector which is to be multiplied by the matrix.
    \param[in] z The x-axis element of the vector which is to be multiplied by the matrix.
    \param[out] ret The vector resulting from matrix multiplication of mat and the vector in this order.
 */
proshade_double* ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( proshade_double* mat, proshade_double x, proshade_double y, proshade_double z )
{
    //================================================ Allocate memory
    proshade_double* ret                              = new proshade_double[3];
    ProSHADE_internal_misc::checkMemoryAllocation     ( ret, __FILE__, __LINE__, __func__ );
    
    //================================================ Compute the multiplication
    ret[0]                                            = ( x * mat[0] ) + ( y * mat[1] ) + ( z * mat[2] );
    ret[1]                                            = ( x * mat[3] ) + ( y * mat[4] ) + ( z * mat[5] );
    ret[2]                                            = ( x * mat[6] ) + ( y * mat[7] ) + ( z * mat[8] );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief Function for computing a 3x3 matrix to 3x1 vector multiplication.
 
    \param[in] mat The matrix to multiply the vector with..
    \param[in] x The x-axis element of the vector which is to be multiplied by the matrix.
    \param[in] y The x-axis element of the vector which is to be multiplied by the matrix.
    \param[in] z The x-axis element of the vector which is to be multiplied by the matrix.
    \param[out] ret The vector resulting from matrix multiplication of mat and the vector in this order.
 */
proshade_single* ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( proshade_single* mat, proshade_single x, proshade_single y, proshade_single z )
{
    //================================================ Allocate memory
    proshade_single* ret                              = new proshade_single[3];
    ProSHADE_internal_misc::checkMemoryAllocation     ( ret, __FILE__, __LINE__, __func__ );
    
    //================================================ Compute the multiplication
    ret[0]                                            = ( x * mat[0] ) + ( y * mat[1] ) + ( z * mat[2] );
    ret[1]                                            = ( x * mat[3] ) + ( y * mat[4] ) + ( z * mat[5] );
    ret[2]                                            = ( x * mat[6] ) + ( y * mat[7] ) + ( z * mat[8] );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief Function for computing a 3x3 matrix inverse.
 
    \param[in] mat The matrix to be inverted.
    \param[out] inverse The inverse of matrix mat.
 */
proshade_double* ProSHADE_internal_maths::compute3x3MatrixInverse ( proshade_double* mat )
{
    //================================================ Allocate memory
    proshade_double* inverse                          = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( inverse, __FILE__, __LINE__, __func__ );
    
    //================================================ Compute determinant
    proshade_double matDet                            = ( mat[0] * mat[4] * mat[8] ) +
                                                        ( mat[1] * mat[5] * mat[6] ) +
                                                        ( mat[2] * mat[3] * mat[7] ) -
                                                        ( mat[0] * mat[5] * mat[7] ) -
                                                        ( mat[1] * mat[3] * mat[8] ) -
                                                        ( mat[2] * mat[4] * mat[6] );
    
    //================================================ Compute inverse matrix
    inverse[0]                                        = ( mat[4] * mat[8] - mat[5] * mat[7] ) / matDet;
    inverse[1]                                        = ( mat[2] * mat[7] - mat[1] * mat[8] ) / matDet;
    inverse[2]                                        = ( mat[1] * mat[5] - mat[2] * mat[4] ) / matDet;
    inverse[3]                                        = ( mat[5] * mat[6] - mat[3] * mat[8] ) / matDet;
    inverse[4]                                        = ( mat[0] * mat[8] - mat[2] * mat[6] ) / matDet;
    inverse[5]                                        = ( mat[2] * mat[3] - mat[0] * mat[5] ) / matDet;
    inverse[6]                                        = ( mat[3] * mat[7] - mat[4] * mat[6] ) / matDet;
    inverse[7]                                        = ( mat[1] * mat[6] - mat[0] * mat[7] ) / matDet;
    inverse[8]                                        = ( mat[0] * mat[4] - mat[1] * mat[3] ) / matDet;
    
    //================================================ Done
    return                                            ( inverse );
}

/*! \brief Transposes 3x3 matrix in place.
 
    \param[in] mat The matrix to be transposed.
 */
void ProSHADE_internal_maths::transpose3x3MatrixInPlace ( proshade_double* mat )
{
    //================================================ Initialise variables
    proshade_double tmp;
    
    //================================================ Transpose the non-diagonal values
    tmp                                               = mat[1];
    mat[1]                                            = mat[3];
    mat[3]                                            = tmp;
    
    tmp                                               = mat[2];
    mat[2]                                            = mat[6];
    mat[6]                                            = tmp;
    
    tmp                                               = mat[5];
    mat[5]                                            = mat[7];
    mat[7]                                            = tmp;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for building a 3x3 matrix from diagonal (and assuming zero padding).
 
    \param[in] diag Array of 3 values that should form the diagonal.
    \param[out] ret The matrix having zeroes everywhere except for the diagonal, where the supplied values will be.
 */
proshade_double* ProSHADE_internal_maths::build3x3MatrixFromDiag ( proshade_double* diag )
{
    //================================================ Allocate memory
    proshade_double* ret                              = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( ret, __FILE__, __LINE__, __func__ );
    
    //================================================ Build
    ret[0]                                            = diag[0];
    ret[1]                                            = 0.0;
    ret[2]                                            = 0.0;
    ret[3]                                            = 0.0;
    ret[4]                                            = diag[1];
    ret[5]                                            = 0.0;
    ret[6]                                            = 0.0;
    ret[7]                                            = 0.0;
    ret[8]                                            = diag[2];
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief Function for building a 3x3 rotation matrix from the x, y and z rotations in degrees.
 
    \param[in] xRot The counter-clockwise rotation about the x axis.
    \param[in] yRot The counter-clockwise rotation about the y axis.
    \param[in] zRot The counter-clockwise rotation about the z axis.
    \param[out] ret The rotation matrix in a vector form.
 */
proshade_double* ProSHADE_internal_maths::build3x3MatrixFromXYZRotations ( proshade_double xRot, proshade_double yRot, proshade_double zRot )
{
    //================================================ Allocate memory
    proshade_double* ret                              = new proshade_double[9];
    proshade_double* XRM                              = new proshade_double[9];
    proshade_double* YRM                              = new proshade_double[9];
    proshade_double* ZRM                              = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( ret, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( XRM, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( YRM, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( ZRM, __FILE__, __LINE__, __func__ );
    
    //================================================ Convert to radians
    proshade_double xRad                              = xRot * (  M_PI / 180.0 );
    proshade_double yRad                              = yRot * (  M_PI / 180.0 );
    proshade_double zRad                              = zRot * (  M_PI / 180.0 );
    
    //================================================ Build the X, Y and Z rotation matrices
    XRM[0] =  1.0;               XRM[1] =  0.0;               XRM[2] =  0.0;
    XRM[3] =  0.0;               XRM[4] =  std::cos ( xRad ); XRM[5] = -std::sin ( xRad );
    XRM[6] =  0.0;               XRM[7] =  std::sin ( xRad ); XRM[8] =  std::cos ( xRad );
    
    YRM[0] =  std::cos ( yRad ); YRM[1] =  0.0;               YRM[2] =  std::sin ( yRad );
    YRM[3] =  0.0;               YRM[4] =  1.0;               YRM[5] =  0.0;
    YRM[6] = -std::sin ( yRad ); YRM[7] =  0.0;               YRM[8] =  std::cos ( yRad );
    
    ZRM[0] =  std::cos ( zRad ); ZRM[1] = -std::sin ( zRad ); ZRM[2] =  0.0;
    ZRM[3] =  std::sin ( zRad ); ZRM[4] =  std::cos ( zRad ); ZRM[5] =  0.0;
    ZRM[6] =  0.0;               ZRM[7] =  0.0;               ZRM[8] =  1.0;
    
    //================================================ Multiply in XYZ order
    proshade_double* tmpMat                           = compute3x3MatrixMultiplication ( XRM,    YRM );
    ret                                               = compute3x3MatrixMultiplication ( tmpMat, ZRM );
    
    //================================================ Release memory
    delete[] tmpMat;
    delete[] XRM;
    delete[] YRM;
    delete[] ZRM;
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief Computation of rotation matrix rotating one vector onto the other.
 
    This function starts by normalising both input vectors to have magnitude of 1.0. Then, it computes the cosine and sine of the angle between the two
    vectors (using the magnitude of the cross product and the dot product); these are then sufficient to build the rotation matrix for rotation in plane on
    which both of the vectors lie.
 
    It then proceeds to compute the change of basis matrix and its inverse, which are in turn sufficient to to compute the rotation matrix in the original
    basis. This rotation matrix is then returned.
 
    \param[in] x1 The x-axis element of the first vector.
    \param[in] y1 The y-axis element of the first vector.
    \param[in] z1 The z-axis element of the first vector.
    \param[in] x2 The x-axis element of the second vector.
    \param[in] y2 The y-axis element of the second vector.
    \param[in] z2 The z-axis element of the second vector.
    \param[out] rotMat Rotation matrix optimally rotating x1 ; y1 ; z1 to match x2 ; y2 ; z2.
 */
proshade_double* ProSHADE_internal_maths::findRotMatMatchingVectors ( proshade_double x1, proshade_double y1, proshade_double z1, proshade_double x2, proshade_double y2, proshade_double z2 )
{
    //================================================ Allocate required memory
    proshade_double* inPlaneRotation                  = new proshade_double[9];
    proshade_double* basisChangeMat                   = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( inPlaneRotation, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( basisChangeMat, __FILE__, __LINE__, __func__ );
    
    //================================================ Normalise inputs
    proshade_double normF                             = std::sqrt( std::pow( x1, 2.0 ) + std::pow ( y1, 2.0 ) + std::pow ( z1, 2.0 ) );
    x1 /= normF; y1 /= normF; z1 /= normF;
    
    normF                                             = std::sqrt( std::pow( x2, 2.0 ) + std::pow ( y2, 2.0 ) + std::pow ( z2, 2.0 ) );
    x2 /= normF; y2 /= normF; z2 /= normF;
    
    //================================================ Compute cross product's magnitude
    proshade_double* crossProd                        = ProSHADE_internal_maths::computeCrossProduct( &x1, &y1, &z1, &x2, &y2, &z2 );
    proshade_double  crossProdMag                     = std::sqrt( std::pow( crossProd[0], 2.0 ) + std::pow ( crossProd[1], 2.0 ) + std::pow ( crossProd[2], 2.0 ) );
    delete[] crossProd;
    
    //================================================ Compute dot product
    proshade_double dotProd                           = ProSHADE_internal_maths::computeDotProduct ( &x1, &y1, &z1, &x2, &y2, &z2 );
    
    //================================================ Construct the in-plane rotation matrix
    inPlaneRotation[0] = dotProd;        inPlaneRotation[1] = -crossProdMag;   inPlaneRotation[2] = 0.0;
    inPlaneRotation[3] = crossProdMag;   inPlaneRotation[4] = dotProd;         inPlaneRotation[5] = 0.0;
    inPlaneRotation[6] = 0.0;            inPlaneRotation[7] = 0.0;             inPlaneRotation[8] = 1.0;
    
    //================================================ Construct change of basis matrix
    crossProd                                         = ProSHADE_internal_maths::computeCrossProduct( &x2, &y2, &z2, &x1, &y1, &z1 );
    normF                                             = std::sqrt ( std::pow ( x2 - ( dotProd * x1 ), 2.0 ) + std::pow ( y2 - ( dotProd * y1 ), 2.0 ) + std::pow ( z2 - ( dotProd * z1 ), 2.0 ) );
    
    basisChangeMat[0] = x1; basisChangeMat[1] = ( x2 - ( dotProd * x1 ) ) / normF; basisChangeMat[2] = crossProd[0];
    basisChangeMat[3] = y1; basisChangeMat[4] = ( y2 - ( dotProd * y1 ) ) / normF; basisChangeMat[5] = crossProd[1];
    basisChangeMat[6] = z1; basisChangeMat[7] = ( z2 - ( dotProd * z1 ) ) / normF; basisChangeMat[8] = crossProd[2];
    
    //================================================ Invert the change of basis matrix
    proshade_double* basisChangeMatInverse            = ProSHADE_internal_maths::compute3x3MatrixInverse ( basisChangeMat );
    
    //================================================ Multiply inverse of change of basis matrix with the in plane rotation matrix, then multiply the result with the inverse
    proshade_double* tmpMat                           = ProSHADE_internal_maths::compute3x3MatrixMultiplication ( basisChangeMat, inPlaneRotation );
    proshade_double* rotMat                           = ProSHADE_internal_maths::compute3x3MatrixMultiplication ( tmpMat, basisChangeMatInverse );
    
    //================================================ Release memory
    delete[] crossProd;
    delete[] inPlaneRotation;
    delete[] basisChangeMat;
    delete[] basisChangeMatInverse;
    delete[] tmpMat;
    
    //================================================ Done
    return                                            ( rotMat );
    
}

/*! \brief This function computes the Moore-Penrose pseudo-inverse of equation I - input matrix.
 
    This function starts by setting and allocating all variables required for singular value decomposition using LAPACK. It also computes the I - input matrix
    result and saves it in column-major format (for Fortran). Next, the singular values are computed and checked for positivity, making use of the fact that
    at least one singular value has to be zero. Finally, the inverse is computed and the resulting matrix is returned.
 
    \warning This function assumes that the input matrix is a rotation matrix. If this assumption does not hold, this function will fail to produce correct results!
 
    \param[in] rMat The rotation matrix to be inverted.
    \param[in] verbose How loud the function should be.
    \param[out] pseudoInverseMat The Moore-Penrose pseudo-inverse of the identity matrix - input matrix.
 */
proshade_double* ProSHADE_internal_maths::compute3x3MoorePenrosePseudoInverseOfIMinusMat ( std::vector < proshade_double >* rMat, proshade_signed verbose )
{
    //================================================ Initialise local variables and allocate the memory
    int dim                                           = 3;                                     // Number of dimensions
    char job                                          = 'A';                                   // Save computation of parts of U and V matrices, they are not needed here
    double *singularValues                            = new double[dim];                       // The array of singular values
    double *rotMatU                                   = new double [dim*dim];                  // The U matrix space
    double *rotMatV                                   = new double [dim*dim];                  // The V^T matrix space
    double *work                                      = new double [static_cast< proshade_unsign >( ( 3 * dim ) + pow( dim, 2 ) * dim)]; // Workspace, minimum required is 4*dim^2 + 7*dim, using more for performance
    int workDim                                       = static_cast< int > ( 2 * ( ( 4 * dim * dim ) + ( 7 * dim ) ) ); // Formalism stating just that
    double* rwork                                     = new double[static_cast<proshade_unsign>((5 * dim) + 5 * pow(dim,2))]; // Required by LAPACK
    int* iwork                                        = new int[(8 * dim)];                    // Required by LAPACK
    int returnValue                                   = 0;                                     // This will tell if operation succeeded
    double *matrixToDecompose                         = new double[dim*dim];
    
    //================================================ Check the memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( singularValues,    __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMatU,           __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMatV,           __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( work,              __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rwork,             __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( iwork,             __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( matrixToDecompose, __FILE__, __LINE__, __func__ );

    //================================================ Load input data into array in column-major order
    for ( int rowIt = 0; rowIt < dim; rowIt++ )
    {
        for ( int colIt = 0; colIt < dim; colIt++ )
        {
            if ( rowIt == colIt ) { matrixToDecompose[(colIt*dim)+rowIt] = 1.0 - rMat->at( static_cast< size_t > ( ( rowIt * dim ) + colIt ) ); }
            else                  { matrixToDecompose[(colIt*dim)+rowIt] = 0.0 - rMat->at( static_cast< size_t > ( ( rowIt * dim ) + colIt ) ); }
        }
    }

    //================================================ Run LAPACK ZGESDD
    dgesdd_                                           ( &job, &dim, &dim, matrixToDecompose, &dim, singularValues, rotMatU, &dim, rotMatV, &dim,
                                                        work, &workDim, rwork, iwork, &returnValue );
    
    //================================================ Check the convergence
    if ( returnValue != 0 )
    {
        ProSHADE_internal_messages::printWarningMessage ( verbose, "!!! ProSHADE WARNING !!! SVD algorithm did not converge.", "WS00069" );
    }

    //================================================ Determine positivity
    bool anyPositive = false;
    std::vector< bool > positivityTest;
    for ( proshade_unsign it = 0; it < static_cast< proshade_unsign > ( dim ); it++ )
    {
        positivityTest.push_back                      ( singularValues[it] > 0.001 );
        if ( positivityTest.at(it) )                  { anyPositive = true; }
    }

    //================================================ Compute according to positivity
    proshade_double* pseudoInverseMat;
    if ( anyPositive )
    {
        //============================================ Set all non-positive singular values and appropriate matrix rows/columns to zero
        if ( !positivityTest.at(0) )
        {
            singularValues[0]                         = 0.0;
            rotMatU[0]                                = 0.0;
            rotMatU[1]                                = 0.0;
            rotMatU[2]                                = 0.0;
            rotMatV[0]                                = 0.0;
            rotMatV[3]                                = 0.0;
            rotMatV[6]                                = 0.0;
        }
        else { singularValues[0] = 1.0 / singularValues[0]; }

        if ( !positivityTest.at(1) )
        {
            singularValues[1]                         = 0.0;
            rotMatU[3]                                = 0.0;
            rotMatU[4]                                = 0.0;
            rotMatU[5]                                = 0.0;
            rotMatV[1]                                = 0.0;
            rotMatV[4]                                = 0.0;
            rotMatV[7]                                = 0.0;
        }
        else { singularValues[1] = 1.0 / singularValues[1]; }

        if ( !positivityTest.at(2) )
        {
            singularValues[2]                         = 0.0;
            rotMatU[6]                                = 0.0;
            rotMatU[7]                                = 0.0;
            rotMatU[8]                                = 0.0;
            rotMatV[2]                                = 0.0;
            rotMatV[5]                                = 0.0;
            rotMatV[8]                                = 0.0;
        }
        else { singularValues[2] = 1.0 / singularValues[2]; }

        //============================================ All positive values formula
        proshade_double* diagMat                      = ProSHADE_internal_maths::build3x3MatrixFromDiag ( singularValues );
        proshade_double* hlpMat                       = ProSHADE_internal_maths::compute3x3MatrixMultiplication ( diagMat, rotMatU );
        pseudoInverseMat                              = ProSHADE_internal_maths::compute3x3MatrixMultiplication ( rotMatV, hlpMat );
        
        //============================================ Release memory
        delete[] diagMat;
        delete[] hlpMat;
    }
    else
    {
        //============================================ No axis in matrix
        pseudoInverseMat                              = new proshade_double[9];
        ProSHADE_internal_misc::checkMemoryAllocation ( pseudoInverseMat, __FILE__, __LINE__, __func__ );

        for ( size_t mIt = 0; mIt < 9; mIt++ ) { pseudoInverseMat[mIt] = 0.0; }
    }

    //================================================ Free memory
    delete[] work;
    delete[] rwork;
    delete[] iwork;
    delete[] matrixToDecompose;
    delete[] singularValues;
    delete[] rotMatU;
    delete[] rotMatV;
    
    //================================================ Done
    return                                            ( pseudoInverseMat );
    
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
    const FloatingPoint< proshade_double > lhs1 ( std::max ( std::abs ( solX ), std::max( std::abs ( solY ), std::abs ( solZ ) ) ) );
    const FloatingPoint< proshade_double > rhs1 ( std::abs ( solX ) );
    const FloatingPoint< proshade_double > rhs2 ( std::abs ( solY ) );
    const FloatingPoint< proshade_double > rhs3 ( std::abs ( solZ ) );
    if ( ( ( lhs1.AlmostEquals ( rhs1 ) ) && ( solX < 0.0 ) ) ||
         ( ( lhs1.AlmostEquals ( rhs2 ) ) && ( solY < 0.0 ) ) ||
         ( ( lhs1.AlmostEquals ( rhs3 ) ) && ( solZ < 0.0 ) ) ) { solX *= -1.0; solY *= -1.0; solZ *= -1.0; }
    
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
    const FloatingPoint< proshade_double > lhs1 ( std::max ( std::abs ( solX ), std::max( std::abs ( solY ), std::abs ( solZ ) ) ) );
    const FloatingPoint< proshade_double > rhs1 ( std::abs ( solX ) );
    const FloatingPoint< proshade_double > rhs2 ( std::abs ( solY ) );
    const FloatingPoint< proshade_double > rhs3 ( std::abs ( solZ ) );
    if ( ( ( lhs1.AlmostEquals ( rhs1 ) ) && ( solX < 0.0 ) ) ||
         ( ( lhs1.AlmostEquals ( rhs2 ) ) && ( solY < 0.0 ) ) ||
         ( ( lhs1.AlmostEquals ( rhs3 ) ) && ( solZ < 0.0 ) ) ) { solX *= -1.0; solY *= -1.0; solZ *= -1.0; }
    
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
    if ( tolerance > std::abs ( trace ) ) { ret = true; }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function compares two vectors using cosine distance and decides if they are similar using tolerance.
 
    This function computes the distance between two vectors, specifically by computing the cosine distance ( ( dot( A, B ) ) / ( mag(A) x mag(B) ) ). This measure will be
    1.0 if the two vectors are identically oriented, 0.0 if they are perpendicular and -1.0 if they have opposite direction. Given that opposite direction must be regarded as
    same direction with opposite angles for symmetry axes detection purposes, this function uses the absolute value of this measure and checks if the two supplied vectors
    (supplied element by element) have cosine distance within 1.0 - tolerance, returning true if they do and false otherwise.
 
    \param[in] a1 The first element of the first vector.
    \param[in] a2 The second element of the first vector.
    \param[in] a3 The third element of the first vector.
    \param[in] b1 The first element of the second vector.
    \param[in] b2 The second element of the second vector.
    \param[in] b3 The third element of the second vector.
    \param[in] tolerance The allowed difference of the distance measure from the 1.0 for the vectors to still be considered similar.
    \param[out] res Boolean decision if the two vectors are similar or not.
 */
bool ProSHADE_internal_maths::vectorOrientationSimilarity ( proshade_double a1, proshade_double a2, proshade_double a3, proshade_double b1, proshade_double b2, proshade_double b3, proshade_double tolerance )
{
    //================================================ Initialise variables
    bool ret                                          = false;
    
    //================================================ Cosine distance
    proshade_double cosDist                           = ( ( a1 * b1 ) + ( a2 * b2 ) + ( a3 * b3 ) ) /
                                                        ( sqrt( pow( a1, 2.0 ) + pow( a2, 2.0 ) + pow( a3, 2.0 ) ) *
                                                          sqrt( pow( b1, 2.0 ) + pow( b2, 2.0 ) + pow( b3, 2.0 ) ) );
    
    //================================================ Compare the absolute value of distance to 1.0 - tolerance
    if ( std::abs( cosDist ) > ( 1.0 - tolerance ) ) { ret = true; }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function compares two vectors using cosine distance and decides if they are similar using tolerance.
 
    This function computes the distance between two vectors, specifically by computing the cosine distance ( ( dot( A, B ) ) / ( mag(A) x mag(B) ) ). This measure will be
    1.0 if the two vectors are identically oriented, 0.0 if they are perpendicular and -1.0 if they have opposite direction. Given that opposite direction must be regarded as
    different vector for peak detection purposes (spheres with different angles are covered separately), this function does not use the absolute value of this measure and
    checks if the two supplied vectors (supplied element by element) have cosine distance within 1.0 - tolerance, returning true if they do and false otherwise.
 
    \param[in] a1 The first element of the first vector.
    \param[in] a2 The second element of the first vector.
    \param[in] a3 The third element of the first vector.
    \param[in] b1 The first element of the second vector.
    \param[in] b2 The second element of the second vector.
    \param[in] b3 The third element of the second vector.
    \param[in] tolerance The allowed difference of the distance measure from the 1.0 for the vectors to still be considered similar.
    \param[out] res Boolean decision if the two vectors are similar or not.
 */
bool ProSHADE_internal_maths::vectorOrientationSimilaritySameDirection ( proshade_double a1, proshade_double a2, proshade_double a3, proshade_double b1, proshade_double b2, proshade_double b3, proshade_double tolerance )
{
    //================================================ Initialise variables
    bool ret                                          = false;
    
    //================================================ Cosine distance
    proshade_double cosDist                           = ( ( a1 * b1 ) + ( a2 * b2 ) + ( a3 * b3 ) ) /
                                                        ( sqrt( pow( a1, 2.0 ) + pow( a2, 2.0 ) + pow( a3, 2.0 ) ) *
                                                          sqrt( pow( b1, 2.0 ) + pow( b2, 2.0 ) + pow( b3, 2.0 ) ) );
    
    //================================================ Compare the absolute value of distance to 1.0 - tolerance
    if ( cosDist > ( 1.0 - tolerance ) ) { ret = true; }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function provides axis optimisation given starting lattitude and longitude indices.
 
    This function takes the initial lattitude and longitude indices as well as the current best sum over all appropriate spheres and the list of the spheres and proceeds to
    use bi-cubic interpolation and a sort of gradient ascend algorithm to search the space around the given indices for interpolated values, which would have higher
    sum of the rotation function values than the initial position. If any improvement is found, it will over-write the input variables.
 
    \param[in] bestLattitude Proshade double pointer to variable containing the best lattitude index value and to which the optimised result will be saved into.
    \param[in] bestLongitude Proshade double pointer to variable containing the best longitude index value and to which the optimised result will be saved into.
    \param[in] bestSum Proshade double pointer to variable containing the best position rotation function values sum and to which the optimised result will be saved into.
    \param[in] sphereList A vector containing the list of spheres which form the set for this symmetry.
    \param[in] step The size of the step.
 */
void ProSHADE_internal_maths::optimiseAxisBiCubicInterpolation ( proshade_double* bestLattitude, proshade_double* bestLongitude, proshade_double* bestSum, std::vector<proshade_unsign>* sphereList, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*>* sphereMappedRotFun, proshade_double step )
{
    //================================================ Initialise variables
    proshade_double lonM, lonP, latM, latP, movSum;
    std::vector<proshade_double> latVals              ( 3 );
    std::vector<proshade_double> lonVals              ( 3 );
    proshade_double learningRate                      = 0.1;
    proshade_double prevVal                           = *bestSum;
    proshade_double valChange                         = 999.9;
    proshade_double origBestLat                       = std::round ( *bestLattitude );
    proshade_double origBestLon                       = std::round ( *bestLongitude );
    proshade_double tmpVal;
    
    //================================================ Initialise interpolators in all directions around the point of interest
    std::vector<ProSHADE_internal_maths::BicubicInterpolator*> interpolsMinusMinus;
    std::vector<ProSHADE_internal_maths::BicubicInterpolator*> interpolsMinusPlus;
    std::vector<ProSHADE_internal_maths::BicubicInterpolator*> interpolsPlusMinus;
    std::vector<ProSHADE_internal_maths::BicubicInterpolator*> interpolsPlusPlus;
    prepareBiCubicInterpolatorsMinusMinus             ( std::round ( *bestLattitude ), std::round ( *bestLongitude ), sphereList, &interpolsMinusMinus, sphereMappedRotFun );
    prepareBiCubicInterpolatorsMinusPlus              ( std::round ( *bestLattitude ), std::round ( *bestLongitude ), sphereList, &interpolsMinusPlus,  sphereMappedRotFun );
    prepareBiCubicInterpolatorsPlusMinus              ( std::round ( *bestLattitude ), std::round ( *bestLongitude ), sphereList, &interpolsPlusMinus,  sphereMappedRotFun );
    prepareBiCubicInterpolatorsPlusPlus               ( std::round ( *bestLattitude ), std::round ( *bestLongitude ), sphereList, &interpolsPlusPlus,   sphereMappedRotFun );
    
    //================================================ Start the pseudo gradient ascent (while there is some change)
    while ( valChange > 0.0001 )
    {
        //============================================ Find the surrounding points to the currently best position
        lonM                                          = *bestLongitude - step;
        lonP                                          = *bestLongitude + step;
        latM                                          = *bestLattitude - step;
        latP                                          = *bestLattitude + step;
        
        //============================================ Deal with optimising outside of prepared range - recursion
        const FloatingPoint< proshade_double > lhs1 ( *bestLattitude ), rhs1 ( origBestLat - 1.0 );
        const FloatingPoint< proshade_double > lhs2 ( *bestLattitude ), rhs2 ( origBestLat + 1.0 );
        const FloatingPoint< proshade_double > lhs3 ( *bestLongitude ), rhs3 ( origBestLon - 1.0 );
        const FloatingPoint< proshade_double > lhs4 ( *bestLongitude ), rhs4 ( origBestLon + 1.0 );
        if ( latM < ( origBestLat - 1.0 ) ) { tmpVal = *bestLattitude; *bestLattitude = origBestLat - 1.0; optimiseAxisBiCubicInterpolation ( bestLattitude, bestLongitude, bestSum, sphereList, sphereMappedRotFun, step ); if ( lhs1.AlmostEquals ( rhs1 ) ) { *bestLattitude = tmpVal; } break; }
        if ( latP > ( origBestLat + 1.0 ) ) { tmpVal = *bestLattitude; *bestLattitude = origBestLat + 1.0; optimiseAxisBiCubicInterpolation ( bestLattitude, bestLongitude, bestSum, sphereList, sphereMappedRotFun, step ); if ( lhs2.AlmostEquals ( rhs2 ) ) { *bestLattitude = tmpVal; } break; }
        if ( lonM < ( origBestLon - 1.0 ) ) { tmpVal = *bestLongitude; *bestLongitude = origBestLon - 1.0; optimiseAxisBiCubicInterpolation ( bestLattitude, bestLongitude, bestSum, sphereList, sphereMappedRotFun, step ); if ( lhs3.AlmostEquals ( rhs3 ) ) { *bestLongitude = tmpVal; } break; }
        if ( lonP > ( origBestLon + 1.0 ) ) { tmpVal = *bestLongitude; *bestLongitude = origBestLon + 1.0; optimiseAxisBiCubicInterpolation ( bestLattitude, bestLongitude, bestSum, sphereList, sphereMappedRotFun, step ); if ( lhs4.AlmostEquals ( rhs4 ) ) { *bestLongitude = tmpVal; } break; }

        //============================================ Prepare vectors of tested positions
        latVals.at(0) = latM; latVals.at(1) = *bestLattitude; latVals.at(2) = latP;
        lonVals.at(0) = lonM; lonVals.at(1) = *bestLongitude; lonVals.at(2) = lonP;
        
        //============================================ Find the best change
        for ( proshade_unsign laIt = 0; laIt < static_cast<proshade_unsign> ( latVals.size() ); laIt++ )
        {
            for ( proshade_unsign loIt = 0; loIt < static_cast<proshade_unsign> ( lonVals.size() ); loIt++ )
            {
                //==================================== For this combination of lat and lon, find sum over spheres
                movSum                                = 1.0;
                for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( sphereList->size() ); iter++ )
                {
                    //================================ Interpolate using correct interpolators
                    if ( ( latVals.at(laIt) <= origBestLat ) && ( lonVals.at(loIt) <= origBestLon ) ) { movSum += interpolsMinusMinus.at(iter)->getValue ( latVals.at(laIt), lonVals.at(loIt) ); }
                    if ( ( latVals.at(laIt) <= origBestLat ) && ( lonVals.at(loIt) >  origBestLon ) ) { movSum += interpolsMinusPlus.at(iter)->getValue  ( latVals.at(laIt), lonVals.at(loIt) ); }
                    if ( ( latVals.at(laIt) >  origBestLat ) && ( lonVals.at(loIt) <= origBestLon ) ) { movSum += interpolsPlusMinus.at(iter)->getValue  ( latVals.at(laIt), lonVals.at(loIt) ); }
                    if ( ( latVals.at(laIt) >  origBestLat ) && ( lonVals.at(loIt) >  origBestLon ) ) { movSum += interpolsPlusPlus.at(iter)->getValue   ( latVals.at(laIt), lonVals.at(loIt) ); }
                }
                
                //==================================== If position has improved, save it
                if ( *bestSum < movSum )
                {
                   *bestSum                           = movSum;
                   *bestLongitude                     = lonVals.at(loIt);
                   *bestLattitude                     = latVals.at(laIt);
                }
            }
        }
        
        //============================================ Prepare for next iteration
        valChange                                     = std::floor ( 100000.0 * ( *bestSum - prevVal ) ) / 100000.0;
        prevVal                                       = std::floor ( 100000.0 * ( *bestSum           ) ) / 100000.0;
        step                                          = std::max ( ( valChange / step ) * learningRate, 0.01 );
        if ( learningRate >= 0.02 ) { learningRate -= 0.01; }
        
    }
    
    //================================================ Release interpolators memory
    for ( proshade_unsign intIt = 0; intIt < static_cast<proshade_unsign> ( interpolsMinusMinus.size() ); intIt++ ) { delete interpolsMinusMinus.at(intIt); }
    for ( proshade_unsign intIt = 0; intIt < static_cast<proshade_unsign> ( interpolsMinusPlus.size()  ); intIt++ ) { delete interpolsMinusPlus.at(intIt);  }
    for ( proshade_unsign intIt = 0; intIt < static_cast<proshade_unsign> ( interpolsPlusMinus.size()  ); intIt++ ) { delete interpolsPlusMinus.at(intIt);  }
    for ( proshade_unsign intIt = 0; intIt < static_cast<proshade_unsign> ( interpolsPlusPlus.size()   ); intIt++ ) { delete interpolsPlusPlus.at(intIt);   }
    
    //================================================ Done
    return ;
}

/*! \brief This function prepares the interpolation objects for the bi-cubic interpolation.
 
    This function takes the position around which the interpolation is to be done and proceeds to create the interpolator
    objects for bi-cubic interpolation in the -- direction (i.e. when both interpolated values will be lower than the best
    lattitude and longitude) using the correct spheres in the correct ranges.
 
    \param[in] bestLattitude The lattitude index value around which interpolation is to be prepared.
    \param[in] bestLongitude The longitude index value around which interpolation is to be prepared.
    \param[in] sphereList A vector containing the list of spheres which form the set for this symmetry.
    \param[in] interpols A pointer to a vector of ProSHADE interpolator objects to which the interpolators will be saved into.
 */
void ProSHADE_internal_maths::prepareBiCubicInterpolatorsMinusMinus ( proshade_double bestLattitude, proshade_double bestLongitude, std::vector<proshade_unsign>* sphereList, std::vector<ProSHADE_internal_maths::BicubicInterpolator*>* interpols, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*>* sphereMappedRotFun )
{
    //================================================ Initialise local variables
    proshade_signed latHlp, lonHlp;
    proshade_signed angDim                            = static_cast< proshade_signed > ( sphereMappedRotFun->at(0)->getAngularDim() );
    
    //================================================ Prepare the interpolator objects for interpolation around the position
    for ( proshade_unsign sphereIt = 0; sphereIt < static_cast<proshade_unsign> ( sphereList->size() ); sphereIt++ )
    {
        //============================================ Allocate memory for the value grid on which the interpolation is to be done (along first dimension)
        proshade_double** interpGrid                  = new proshade_double*[4];
        ProSHADE_internal_misc::checkMemoryAllocation ( interpGrid, __FILE__, __LINE__, __func__ );

        //============================================ Allocate memory for the value grid on which the interpolation is to be done (along second dimension)
        for ( proshade_unsign iter = 0; iter < 4; iter++ )
        {
            interpGrid[iter]                          = new proshade_double[4];
            ProSHADE_internal_misc::checkMemoryAllocation ( interpGrid[iter], __FILE__, __LINE__, __func__ );
        }

        //============================================ Fill in the value grid on which the interpolation is to be done
        for ( proshade_signed latIt = 0; latIt < 4; latIt++ )
        {
            for ( proshade_signed lonIt = 0; lonIt < 4; lonIt++ )
            {
                latHlp = static_cast< proshade_signed > ( bestLattitude - 2.0 + static_cast< proshade_double > ( latIt ) ); if ( latHlp < 0 ) { latHlp += angDim; } if ( latHlp >= angDim ) { latHlp -= angDim; }
                lonHlp = static_cast< proshade_signed > ( bestLongitude - 2.0 + static_cast< proshade_double > ( lonIt ) ); if ( lonHlp < 0 ) { lonHlp += angDim; } if ( lonHlp >= angDim ) { lonHlp -= angDim; }
                interpGrid[latIt][lonIt]              = sphereMappedRotFun->at(sphereList->at(sphereIt))->getSphereLatLonPosition ( static_cast< proshade_unsign > ( latHlp ), static_cast< proshade_unsign > ( lonHlp ) );
            }
        }

        //============================================ Create the interpolators
        ProSHADE_internal_maths::BicubicInterpolator* biCubInterp = new ProSHADE_internal_maths::BicubicInterpolator ( interpGrid, bestLattitude - 1.0, bestLongitude - 1.0 );
        interpols->emplace_back                       ( biCubInterp );
        
        //============================================ Release memory
        for ( proshade_unsign iter = 0; iter < 4; iter++ ) { delete[] interpGrid[iter]; }
        delete[] interpGrid;
    }
    
    //================================================ Done
    return ;
}

/*! \brief This function prepares the interpolation objects for the bi-cubic interpolation.
 
    This function takes the position around which the interpolation is to be done and proceeds to create the interpolator
    objects for bi-cubic interpolation in the -+ direction (i.e. when interpolated lattitude is lower than best lattitude, but
    interpolated longitude is higher than best longitude) using the correct spheres in the correct ranges.
 
    \param[in] bestLattitude The lattitude index value around which interpolation is to be prepared.
    \param[in] bestLongitude The longitude index value around which interpolation is to be prepared.
    \param[in] sphereList A vector containing the list of spheres which form the set for this symmetry.
    \param[in] interpols A pointer to a vector of ProSHADE interpolator objects to which the interpolators will be saved into.
 */
void ProSHADE_internal_maths::prepareBiCubicInterpolatorsMinusPlus ( proshade_double bestLattitude, proshade_double bestLongitude, std::vector<proshade_unsign>* sphereList, std::vector<ProSHADE_internal_maths::BicubicInterpolator*>* interpols, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*>* sphereMappedRotFun )
{
    //================================================ Initialise local variables
    proshade_signed latHlp, lonHlp;
    proshade_signed angDim                            = static_cast< proshade_signed > ( sphereMappedRotFun->at(0)->getAngularDim() );
    
    //================================================ Prepare the interpolator objects for interpolation around the position
    for ( proshade_unsign sphereIt = 0; sphereIt < static_cast<proshade_unsign> ( sphereList->size() ); sphereIt++ )
    {
        //============================================ Allocate memory for the value grid on which the interpolation is to be done (along first dimension)
        proshade_double** interpGrid                  = new proshade_double*[4];
        ProSHADE_internal_misc::checkMemoryAllocation ( interpGrid, __FILE__, __LINE__, __func__ );

        //============================================ Allocate memory for the value grid on which the interpolation is to be done (along second dimension)
        for ( proshade_unsign iter = 0; iter < 4; iter++ )
        {
            interpGrid[iter]                          = new proshade_double[4];
            ProSHADE_internal_misc::checkMemoryAllocation ( interpGrid[iter], __FILE__, __LINE__, __func__ );
        }

        //============================================ Fill in the value grid on which the interpolation is to be done
        for ( proshade_unsign latIt = 0; latIt < 4; latIt++ )
        {
            for ( proshade_unsign lonIt = 0; lonIt < 4; lonIt++ )
            {
                latHlp = static_cast< proshade_signed > ( bestLattitude - 2 + static_cast< proshade_double > ( latIt ) ); if ( latHlp < 0 ) { latHlp += angDim; } if ( latHlp >= angDim ) { latHlp -= angDim; }
                lonHlp = static_cast< proshade_signed > ( bestLongitude - 1 + static_cast< proshade_double > ( lonIt ) ); if ( lonHlp < 0 ) { lonHlp += angDim; } if ( lonHlp >= angDim ) { lonHlp -= angDim; }
                interpGrid[latIt][lonIt]              = sphereMappedRotFun->at(sphereList->at(sphereIt))->getSphereLatLonPosition ( static_cast< proshade_unsign > ( latHlp ) , static_cast< proshade_unsign > ( lonHlp ) );
            }
        }

        //============================================ Create the interpolators
        ProSHADE_internal_maths::BicubicInterpolator* biCubInterp = new ProSHADE_internal_maths::BicubicInterpolator ( interpGrid, bestLattitude - 1.0, bestLongitude );
        interpols->emplace_back                       ( biCubInterp );
        
        //============================================ Release memory
        for ( proshade_unsign iter = 0; iter < 4; iter++ ) { delete[] interpGrid[iter]; }
        delete[] interpGrid;
    }
    
    //================================================ Done
    return ;
}

/*! \brief This function prepares the interpolation objects for the bi-cubic interpolation.
 
    This function takes the position around which the interpolation is to be done and proceeds to create the interpolator
    objects for bi-cubic interpolation in the +- direction (i.e. when interpolated lattitude is higher than best lattitude, but
    interpolated longitude is lower than best longitude) using the correct spheres in the correct ranges.
 
    \param[in] bestLattitude The lattitude index value around which interpolation is to be prepared.
    \param[in] bestLongitude The longitude index value around which interpolation is to be prepared.
    \param[in] sphereList A vector containing the list of spheres which form the set for this symmetry.
    \param[in] interpols A pointer to a vector of ProSHADE interpolator objects to which the interpolators will be saved into.
 */
void ProSHADE_internal_maths::prepareBiCubicInterpolatorsPlusMinus ( proshade_double bestLattitude, proshade_double bestLongitude, std::vector<proshade_unsign>* sphereList, std::vector<ProSHADE_internal_maths::BicubicInterpolator*>* interpols, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*>* sphereMappedRotFun )
{
    //================================================ Initialise local variables
    proshade_signed latHlp, lonHlp;
    proshade_signed angDim                            = static_cast< proshade_signed > ( sphereMappedRotFun->at(0)->getAngularDim() );
    
    //================================================ Prepare the interpolator objects for interpolation around the position
    for ( proshade_unsign sphereIt = 0; sphereIt < static_cast<proshade_unsign> ( sphereList->size() ); sphereIt++ )
    {
        //============================================ Allocate memory for the value grid on which the interpolation is to be done (along first dimension)
        proshade_double** interpGrid                  = new proshade_double*[4];
        ProSHADE_internal_misc::checkMemoryAllocation ( interpGrid, __FILE__, __LINE__, __func__ );

        //============================================ Allocate memory for the value grid on which the interpolation is to be done (along second dimension)
        for ( proshade_unsign iter = 0; iter < 4; iter++ )
        {
            interpGrid[iter]                          = new proshade_double[4];
            ProSHADE_internal_misc::checkMemoryAllocation ( interpGrid[iter], __FILE__, __LINE__, __func__ );
        }

        //============================================ Fill in the value grid on which the interpolation is to be done
        for ( proshade_unsign latIt = 0; latIt < 4; latIt++ )
        {
            for ( proshade_unsign lonIt = 0; lonIt < 4; lonIt++ )
            {
                latHlp = static_cast< proshade_signed > ( bestLattitude - 1 + static_cast< proshade_double > ( latIt ) ); if ( latHlp < 0 ) { latHlp += angDim; } if ( latHlp >= angDim ) { latHlp -= angDim; }
                lonHlp = static_cast< proshade_signed > ( bestLongitude - 2 + static_cast< proshade_double > ( lonIt ) ); if ( lonHlp < 0 ) { lonHlp += angDim; } if ( lonHlp >= angDim ) { lonHlp -= angDim; }
                interpGrid[latIt][lonIt]              = sphereMappedRotFun->at(sphereList->at(sphereIt))->getSphereLatLonPosition ( static_cast< proshade_unsign > ( latHlp ), static_cast< proshade_unsign > ( lonHlp ) );
            }
        }

        //============================================ Create the interpolators
        ProSHADE_internal_maths::BicubicInterpolator* biCubInterp = new ProSHADE_internal_maths::BicubicInterpolator ( interpGrid, bestLattitude, bestLongitude - 1.0 );
        interpols->emplace_back                       ( biCubInterp );
        
        //============================================ Release memory
        for ( proshade_unsign iter = 0; iter < 4; iter++ ) { delete[] interpGrid[iter]; }
        delete[] interpGrid;
    }
    
    //================================================ Done
    return ;
}

/*! \brief This function prepares the interpolation objects for the bi-cubic interpolation.
 
    This function takes the position around which the interpolation is to be done and proceeds to create the interpolator
    objects for bi-cubic interpolation in the ++ direction (i.e. when both interpolated values will be larger than the best
    lattitude and longitude) using the correct spheres in the correct ranges.
 
    \param[in] bestLattitude The lattitude index value around which interpolation is to be prepared.
    \param[in] bestLongitude The longitude index value around which interpolation is to be prepared.
    \param[in] sphereList A vector containing the list of spheres which form the set for this symmetry.
    \param[in] interpols A pointer to a vector of ProSHADE interpolator objects to which the interpolators will be saved into.
 */
void ProSHADE_internal_maths::prepareBiCubicInterpolatorsPlusPlus ( proshade_double bestLattitude, proshade_double bestLongitude, std::vector<proshade_unsign>* sphereList, std::vector<ProSHADE_internal_maths::BicubicInterpolator*>* interpols, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*>* sphereMappedRotFun )
{
    //================================================ Initialise local variables
    proshade_signed latHlp, lonHlp;
    proshade_signed angDim                            = static_cast< proshade_signed > ( sphereMappedRotFun->at(0)->getAngularDim() );
    
    //================================================ Prepare the interpolator objects for interpolation around the position
    for ( proshade_unsign sphereIt = 0; sphereIt < static_cast<proshade_unsign> ( sphereList->size() ); sphereIt++ )
    {
        //============================================ Allocate memory for the value grid on which the interpolation is to be done (along first dimension)
        proshade_double** interpGrid                  = new proshade_double*[4];
        ProSHADE_internal_misc::checkMemoryAllocation ( interpGrid, __FILE__, __LINE__, __func__ );

        //============================================ Allocate memory for the value grid on which the interpolation is to be done (along second dimension)
        for ( proshade_unsign iter = 0; iter < 4; iter++ )
        {
            interpGrid[iter]                          = new proshade_double[4];
            ProSHADE_internal_misc::checkMemoryAllocation ( interpGrid[iter], __FILE__, __LINE__, __func__ );
        }

        //============================================ Fill in the value grid on which the interpolation is to be done
        for ( proshade_unsign latIt = 0; latIt < 4; latIt++ )
        {
            for ( proshade_unsign lonIt = 0; lonIt < 4; lonIt++ )
            {
                latHlp = static_cast< proshade_signed > ( bestLattitude - 1 + static_cast< proshade_double > ( latIt ) ); if ( latHlp < 0 ) { latHlp += angDim; } if ( latHlp >= angDim ) { latHlp -= angDim; }
                lonHlp = static_cast< proshade_signed > ( bestLongitude - 1 + static_cast< proshade_double > ( lonIt ) ); if ( lonHlp < 0 ) { lonHlp += angDim; } if ( lonHlp >= angDim ) { lonHlp -= angDim; }
                interpGrid[latIt][lonIt]              = sphereMappedRotFun->at(sphereList->at(sphereIt))->getSphereLatLonPosition ( static_cast< proshade_unsign > ( latHlp ), static_cast< proshade_unsign > ( lonHlp ) );
            }
        }

        //============================================ Create the interpolators
        ProSHADE_internal_maths::BicubicInterpolator* biCubInterp = new ProSHADE_internal_maths::BicubicInterpolator ( interpGrid, bestLattitude, bestLongitude );
        interpols->emplace_back                       ( biCubInterp );
        
        //============================================ Release memory
        for ( proshade_unsign iter = 0; iter < 4; iter++ ) { delete[] interpGrid[iter]; }
        delete[] interpGrid;
    }
    
    //================================================ Done
    return ;
}

/*! \brief This function checks if new axis is unique, or already detected.
 
    This function compares the supplied axis against all members of the axes vector. If the axis has the same fold and very similar
    axis vector (i.e. all three elements are within tolerance), then the function returns false. If no such match is found, true is returned.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] axis The axis to be checked against CSymList to see if it not already present.
    \param[in] tolerance The allowed error on each dimension of the axis.
    \param[in] improve If a similar axis is found and if this already existing axis has lower peak height, should the CSymList be updated with the higher peak height axis?
    \param[out] ret Boolean specifying whether a similar axis was found or not.
 */
bool ProSHADE_internal_maths::isAxisUnique ( std::vector< proshade_double* >* CSymList, proshade_double* axis, proshade_double tolerance, bool improve )
{
    //================================================ Initialise variables
    bool ret                                          = true;
    proshade_unsign whichImprove                      = 0;
    
    //================================================ For each already detected member
    for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( CSymList->size() ); grIt++ )
    {
        //============================================ Is fold the same?
        const FloatingPoint< proshade_double > lhs ( CSymList->at(grIt)[0] ), rhs ( axis[0] );
        if ( lhs.AlmostEquals ( rhs ) )
        {
            if ( ProSHADE_internal_maths::vectorOrientationSimilarity ( CSymList->at(grIt)[1], CSymList->at(grIt)[2], CSymList->at(grIt)[3], axis[1], axis[2], axis[3], tolerance ) )
            {
                ret                                   = false;
                whichImprove                          = grIt;
                break;
            }
        }
    }
    
    //================================================ Improve, if required
    if ( improve && !ret )
    {
        if ( axis[5] > CSymList->at(whichImprove)[5] )
        {
            CSymList->at(whichImprove)[1]             = axis[1];
            CSymList->at(whichImprove)[2]             = axis[2];
            CSymList->at(whichImprove)[3]             = axis[3];
            CSymList->at(whichImprove)[4]             = axis[4];
            CSymList->at(whichImprove)[5]             = axis[5];
        }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function checks if new axis is unique, or already detected.
 
    This function compares the supplied axis against all members of the axes vector. If the axis has the same fold and very similar
    axis vector (i.e. all three elements are within tolerance), then the function returns false. If no such match is found, true is returned.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] X The axis x-element to be checked against CSymList to see if it not already present.
    \param[in] Y The axis x-element to be checked against CSymList to see if it not already present.
    \param[in] Z The axis x-element to be checked against CSymList to see if it not already present.
    \param[in] tolerance The allowed error on each dimension of the axis.
    \param[out] ret Boolean specifying whether a similar axis was found or not.
 */
bool ProSHADE_internal_maths::isAxisUnique ( std::vector< proshade_double* >* CSymList, proshade_double X, proshade_double Y, proshade_double Z, proshade_double fold, proshade_double tolerance )
{
    //================================================ Initialise variables
    bool ret                                          = true;
    
    //================================================ For each already detected member
    for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( CSymList->size() ); grIt++ )
    {
        const FloatingPoint< proshade_double > lhs ( fold ), rhs ( CSymList->at(grIt)[0] );
        if ( lhs.AlmostEquals ( rhs ) )
        {
            if ( ProSHADE_internal_maths::vectorOrientationSimilarity ( CSymList->at(grIt)[1], CSymList->at(grIt)[2], CSymList->at(grIt)[3], X, Y, Z, tolerance ) )
            {
                ret                                       = false;
                break;
            }
        }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function checks if new axis is unique, or already detected and returns the position of match or -1.
 
    This function is a variation on the isAxisUnique function, returning the index of the match or -1.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] axis The axis to be checked against CSymList to see if it not already present.
    \param[in] tolerance The allowed error on each dimension of the axis.
    \param[in] improve If a similar axis is found and if this already existing axis has lower peak height, should the CSymList be updated with the higher peak height axis?
    \param[out] ret Index of the match or -1.
 */
proshade_signed ProSHADE_internal_maths::whichAxisUnique ( std::vector< proshade_double* >* CSymList, proshade_double* axis, proshade_double tolerance )
{
    //================================================ Initialise variables
    proshade_signed ret                               = -1;
    
    //================================================ For each already detected member
    for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( CSymList->size() ); grIt++ )
    {
        //============================================ Is fold the same?
        const FloatingPoint< proshade_double > lhs ( CSymList->at(grIt)[0] ), rhs ( axis[0] );
        if ( lhs.AlmostEquals ( rhs ) )
        {
            if ( ProSHADE_internal_maths::vectorOrientationSimilarity ( CSymList->at(grIt)[1], CSymList->at(grIt)[2], CSymList->at(grIt)[3], axis[1], axis[2], axis[3], tolerance ) )
            {
                ret                                   = static_cast< proshade_signed > ( grIt );
                break;
            }
        }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function checks if new axis is unique, or already detected and returns the position of match or -1.
 
    This function is a variation on the isAxisUnique function, returning the index of the match or -1.
 
    \param[in] CSymList A vector containing the already detected Cyclic symmetries.
    \param[in] X The axis x-element to be checked against CSymList to see if it not already present.
    \param[in] Y The axis x-element to be checked against CSymList to see if it not already present.
    \param[in] Z The axis x-element to be checked against CSymList to see if it not already present.
    \param[in] tolerance The allowed error on each dimension of the axis.
    \param[out] ret Index of the match or -1.
 */
proshade_signed ProSHADE_internal_maths::whichAxisUnique ( std::vector< proshade_double* >* CSymList, proshade_double X, proshade_double Y, proshade_double Z, proshade_double fold, proshade_double tolerance )
{
    //================================================ Initialise variables
    proshade_signed ret                               = -1;
    
    //================================================ For each already detected member
    for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( CSymList->size() ); grIt++ )
    {
        const FloatingPoint< proshade_double > lhs ( fold ), rhs ( CSymList->at(grIt)[0] );
        if ( lhs.AlmostEquals ( rhs ) )
        {
            if ( ProSHADE_internal_maths::vectorOrientationSimilarity ( CSymList->at(grIt)[1], CSymList->at(grIt)[2], CSymList->at(grIt)[3], X, Y, Z, tolerance ) )
            {
                ret                                   = static_cast< proshade_signed > ( grIt );
                break;
            }
        }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function finds all prime numbers up to the supplied limit.
 
    This function uses the sieve of Eratosthenes algorithm to find all prime numbers from 2 to the supplied limit. This is not
    the fastest algorithm and it may become slow when the limit is high, but it is fine for small numbers and given that we
    will use it for symmetry folds, which should not got much over 20, this should be more than fast enough.
 
    \param[in] upTo The limit to which prime numbers should be sought.
 */
std::vector< proshade_unsign > ProSHADE_internal_maths::findAllPrimes ( proshade_unsign upTo )
{
    //================================================ Initialise variables
    std::vector< proshade_unsign > ret;
    std::vector< std::pair< proshade_unsign, bool > > sieveOfEratosthenesArray;
    
    //================================================ Sanity check
    if ( upTo < 2 ) { return ( ret ); }
    
    //================================================ Initialise the sieve array up to the required number
    for ( proshade_unsign iter = 2; iter <= upTo; iter++ ) { sieveOfEratosthenesArray.emplace_back ( std::pair< proshade_unsign, bool > ( iter, true ) ); }
    
    //================================================ For each entry in the array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( sieveOfEratosthenesArray.size() ); iter++ )
    {
        //============================================ If this entry is still true
        if ( sieveOfEratosthenesArray.at(iter).second )
        {
            //======================================== Set all entries with the position x * [this entry value] to false
            for ( proshade_unsign it = iter + sieveOfEratosthenesArray.at(iter).first; it < static_cast<proshade_unsign> ( sieveOfEratosthenesArray.size() ); it += sieveOfEratosthenesArray.at(iter).first )
            {
                sieveOfEratosthenesArray.at(it).second = false;
            }
        }
    }
    
    //================================================ Copy passing results to return vector
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( sieveOfEratosthenesArray.size() ); iter++ )
    {
        if ( sieveOfEratosthenesArray.at(iter).second ) { ProSHADE_internal_misc::addToUnsignVector ( &ret, sieveOfEratosthenesArray.at(iter).first ); }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function computes a Gaussian (normal) distribution value given distance from mean and sigma.
 
    This function simply returns the height of a normal distribution with a given sigma for a value specific distance from the mean.
 
    \param[in] val The distance from the mean for which the Gaussian height should be computed.
    \param[in] sigma The standard deviation of the Gaussian for which the computation is done.
    \param[out] height The height of the Gaussian distribution as desctibed by the sigma.
 */
proshade_double ProSHADE_internal_maths::computeGaussian ( proshade_double val, proshade_double sigma )
{
    //================================================ Compute cumulative probability from Z-score
    proshade_double zScore                            = ( val / sigma );
    proshade_double cumulativeProbability             = 0.5 * std::erfc ( zScore * M_SQRT1_2 );
    
    //================================================ Symmetrise
    if ( cumulativeProbability > 0.5 ) { cumulativeProbability = 1.0 - cumulativeProbability; }
    
    //================================================ Done
    return                                            ( cumulativeProbability );
    
}

/*! \brief This function takes a 1D vector and computes smoothened version based on the parameters.
 
    This function firstly computes the Gaussian weights for each position in a window size accordingly to the parameters amd then
    proceeds to compute the weighted sum for each position created by sliding this window along the data. This results in smoothening
    of the data in accordance with the parameters.
 
    \param[in] step The size of the step on scale from 0.0 to 1.0 including boarders.
    \param[in] windowSize The size of the averaged over window. It is assumed to be odd.
    \param[in] sigma The standard deviation of the Gaussian to be used for smoothening.
    \param[in] data The data to be smoothened.
    \param[out] smoothened A vector of smoothened values for the input data with length hist.size() - (windowSize - 1).
 */
std::vector < proshade_double > ProSHADE_internal_maths::smoothen1D ( proshade_double step, proshade_signed windowSize, proshade_double sigma, std::vector< proshade_double > data )
{
    //================================================ Initialise local variables
    proshade_signed windowHalf                        = ( windowSize - 1 ) / 2;
    proshade_signed totSize                           = static_cast< proshade_signed > ( ( 1.0 / step ) + 1 );
    std::vector< proshade_double > smoothened         ( static_cast< size_t > ( totSize - ( windowSize - 1 ) ), 0.0 );
    std::vector< proshade_double > winWeights         ( static_cast< size_t > ( windowSize ), 0.0 );
    
    //================================================ Prepare window weights
    for ( proshade_double winIt = 0.0; winIt < static_cast< proshade_double > ( windowSize ); winIt += 1.0 ) { winWeights.at( static_cast< proshade_unsign > ( winIt ) ) = ProSHADE_internal_maths::computeGaussian ( ( winIt - static_cast< proshade_double > ( windowHalf ) ) * step, sigma ); }
    
    //================================================ Compute smoothened data
    for ( proshade_unsign it = 0; it < static_cast< proshade_unsign > ( smoothened.size() ); it++ )
    {
        //============================================ Compute window weighted average
        for ( proshade_signed winIt = 0; winIt < windowSize; winIt++ )
        {
            smoothened.at(it)                        += winWeights.at( static_cast< size_t > (  winIt ) ) * data.at( static_cast< size_t > ( static_cast< proshade_signed > ( it ) + winIt ) );
        }
    }
    
    //================================================ Done
    return                                            ( smoothened );
    
}

/*! \brief This function computes the resolution of a particular reflection.
 
    \param[in] h The index of the reflection in reciprocal space along the x-axis.
    \param[in] k The index of the reflection in reciprocal space along the y-axis.
    \param[in] l The index of the reflection in reciprocal space along the z-axis.
    \param[in] xDim The dimension of the cell along the x-axis in Angstroms.
    \param[in] yDim The dimension of the cell along the y-axis in Angstroms.
    \param[in] zDim The dimension of the cell along the z-axis in Angstroms.
    \param[out] ret The resolution of the particular reflection.
 */
proshade_single ProSHADE_internal_maths::getResolutionOfReflection ( proshade_single h, proshade_single k, proshade_single l, proshade_single xDim, proshade_single yDim, proshade_single zDim )
{
    //================================================ Compute volume and proportions
    proshade_single vol                               = ( xDim * yDim * zDim );
    proshade_single sa                                = ( yDim * zDim ) / vol;
    proshade_single sb                                = ( xDim * zDim ) / vol;
    proshade_single sc                                = ( xDim * yDim ) / vol;
    
    //================================================ Compute distance
    proshade_single s2                                = ( std::pow ( h * sa, 2.0f ) +
                                                          std::pow ( k * sb, 2.0f ) +
                                                          std::pow ( l * sc, 2.0f ) ) / 4.0f;
    
    //================================================ Deal with F000
    if ( s2 == 0.0f ) { s2 = 0.0000000001f; }
    
    //================================================ Done
    return                                            ( 1.0f / ( 2.0f * std::sqrt ( s2 ) ) );
    
}

/*! \brief This function does binning of the reciprocal space reflections.
 
    This funcion uses the knowledge of the cell dimensions to firstly decide which dimension is the limitting one in terms
    of FSC binning and then it proceeds to compute the bins, their positions in terms of distance from F000. Finally, it uses
    this information to compute a "mask" map, where each position has the value corresponding to the appropriate bin
    index, thus allowing fast binning of any map with the same dimensions as supplied to this function.
 
    \param[in] xInds The number of indices along the x-axis.
    \param[in] yInds The number of indices along the y-axis.
    \param[in] zInds The number of indices along the z-axis.
    \param[in] noBin Variable to which the number of binds found will be saved into.
    \param[in] binIndexing A pointer to which the map of bin belonging for each reflection will be saved into.
 */
void ProSHADE_internal_maths::binReciprocalSpaceReflections ( proshade_unsign xInds, proshade_unsign yInds, proshade_unsign zInds, proshade_signed* noBin, proshade_signed*& binIndexing )
{
    //================================================ Allocate output bin indexing memory and set to -100
    binIndexing                                       = new proshade_signed [xInds * yInds * zInds];
    ProSHADE_internal_misc::checkMemoryAllocation     ( binIndexing, __FILE__, __LINE__, __func__ );
    for ( size_t iter = 0; iter < static_cast< size_t > ( xInds * yInds * zInds ); iter++ ) { binIndexing[iter] = -100; }
    proshade_single xIndsF = static_cast< proshade_single > ( xInds ), yIndsF = static_cast< proshade_single > ( yInds ), zIndsF = static_cast< proshade_single > ( zInds );
    
    //================================================ Allocate local memory
    proshade_single *mins                             = new proshade_single[3];
    proshade_single *maxs                             = new proshade_single[3];
    proshade_single *resMins                          = new proshade_single[3];
    proshade_signed *resMinLoc                        = new proshade_signed[3];
    proshade_single *steps                            = new proshade_single[3];
    
    //================================================ Check local memory
    ProSHADE_internal_misc::checkMemoryAllocation     ( mins,      __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( maxs,      __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( resMins,   __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( resMinLoc, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( steps,     __FILE__, __LINE__, __func__ );
    
    //================================================ Initialise local variables
    proshade_single resol                             = 0.0f;
    proshade_signed reciX, reciY, reciZ, arrPos = 0, minLoc = -1;
   *noBin                                             = 0;
    
    //================================================ Determine reciprocal space indexing
    mins[0]                                           = std::floor ( xIndsF / -2.0f );
    mins[1]                                           = std::floor ( yIndsF / -2.0f );
    mins[2]                                           = std::floor ( zIndsF / -2.0f );
        
    maxs[0]                                           = -mins[0];
    maxs[1]                                           = -mins[1];
    maxs[2]                                           = -mins[2];
    
    if ( xInds % 2 == 0 ) { mins[0] += 1.0f; }
    if ( yInds % 2 == 0 ) { mins[1] += 1.0f; }
    if ( zInds % 2 == 0 ) { mins[2] += 1.0f; }
    
    //================================================ Get minimum resolution based on dims for each dimension
    resMins[0]                                        = ProSHADE_internal_maths::getResolutionOfReflection ( maxs[0], 0.0f, 0.0f, xIndsF, yIndsF, zIndsF );
    resMins[1]                                        = ProSHADE_internal_maths::getResolutionOfReflection ( 0.0f, maxs[1], 0.0f, xIndsF, yIndsF, zIndsF );
    resMins[2]                                        = ProSHADE_internal_maths::getResolutionOfReflection ( 0.0f, 0.0f, maxs[2], xIndsF, yIndsF, zIndsF );

    //================================================ Decide which dimension to work with (the one with the lowest resolution)
    resMinLoc[0] = 0; resMinLoc[1] = 0; resMinLoc[2] = 0;
    const FloatingPoint< proshade_single > lhs1 ( resMins[0] ), lhs2 ( resMins[1] ), lhs3 ( resMins[2] ), rhs1 ( std::min( resMins[0], std::min( resMins[1], resMins[2] ) ) );
    if ( lhs1.AlmostEquals ( rhs1 ) ) { resMinLoc[0] = 1; minLoc = 0; }
    if ( lhs2.AlmostEquals ( rhs1 ) ) { resMinLoc[1] = 1; minLoc = 1; }
    if ( lhs3.AlmostEquals ( rhs1 ) ) { resMinLoc[2] = 1; minLoc = 2; }
    
    //================================================ Find the bins and corresponding cut-offs
    std::vector< proshade_single > resArray           ( static_cast< size_t > ( maxs[minLoc] - 1 ), 0.0f );
    std::vector< proshade_single > binArray           ( static_cast< size_t > ( maxs[minLoc] - 1 ), 0.0f );
    for ( proshade_signed dimIt = 0; dimIt < static_cast< proshade_signed > ( maxs[minLoc] - 1 ); dimIt++ )
    {
        //============================================ Prepare steps
        steps[0]                                      = ( static_cast< proshade_single > ( dimIt ) + 2.5f ) * static_cast< proshade_single > ( resMinLoc[0] );
        steps[1]                                      = ( static_cast< proshade_single > ( dimIt ) + 2.5f ) * static_cast< proshade_single > ( resMinLoc[1] );
        steps[2]                                      = ( static_cast< proshade_single > ( dimIt ) + 2.5f ) * static_cast< proshade_single > ( resMinLoc[2] );
        
        //============================================ Find resolution
        resol                                         = ProSHADE_internal_maths::getResolutionOfReflection ( steps[0], steps[1], steps[2], xIndsF, yIndsF, zIndsF );
        
        //============================================ Assign to arrays
        resArray.at( static_cast< size_t > ( dimIt ) ) = resol;
        binArray.at( static_cast< size_t > ( dimIt ) ) = static_cast< proshade_single > ( dimIt ) + 2.5f;
       *noBin                                          = dimIt + 1;
    }
    
    //================================================ Assign reflections to bins
    for ( proshade_signed xIt = 0; xIt < static_cast< proshade_signed > ( xInds ); xIt++ )
    {
        for ( proshade_signed yIt = 0; yIt < static_cast< proshade_signed > ( yInds ); yIt++ )
        {
            for ( proshade_signed zIt = 0; zIt < static_cast< proshade_signed > ( zInds / 2 ) + 1; zIt++ )
            {
                //==================================== Deal with reciprocal indices ordering
                reciX = xIt; if ( reciX > static_cast< proshade_signed > ( maxs[0] ) ) { reciX -= static_cast< proshade_signed > ( xInds ); }
                reciY = yIt; if ( reciY > static_cast< proshade_signed > ( maxs[1] ) ) { reciY -= static_cast< proshade_signed > ( yInds ); }
                reciZ = zIt; if ( reciZ > static_cast< proshade_signed > ( maxs[2] ) ) { reciZ -= static_cast< proshade_signed > ( zInds ); }
                
                //==================================== For each bin, check if this reflection belongs to it
                for ( proshade_signed binIt = 0; binIt < (*noBin); binIt++ )
                {
                    //================================ Check by comparing distances
                    if ( std::sqrt ( std::pow ( static_cast< proshade_single > ( reciX ), 2.0f ) +
                                     std::pow ( static_cast< proshade_single > ( reciY ), 2.0f ) +
                                     std::pow ( static_cast< proshade_single > ( reciZ ), 2.0f ) ) <= binArray.at( static_cast< size_t > ( binIt ) ) )
                    {
                        //============================ This is the bin for this reflection. Assign it.
                        arrPos                        = zIt + static_cast< proshade_signed > ( zInds ) * ( yIt + static_cast< proshade_signed > ( yInds ) * xIt );
                        binIndexing[ static_cast< size_t > ( arrPos ) ] = binIt;
                        
                        //============================ If one of the uneven ends, do not use Friedel's Law
                        if ( reciX == static_cast< proshade_signed > ( mins[0] ) || -reciX == static_cast< proshade_signed > ( mins[0] ) ) { break; }
                        if ( reciY == static_cast< proshade_signed > ( mins[1] ) || -reciY == static_cast< proshade_signed > ( mins[1] ) ) { break; }
                        if ( reciZ == static_cast< proshade_signed > ( mins[2] ) || -reciZ == static_cast< proshade_signed > ( mins[2] ) ) { break; }
                        
                        //============================ Use Friedel's Law to find the second index (this is why we can use zDim / 2)
                        reciX *= -1; if ( reciX < 0 ) { reciX += static_cast< proshade_signed > ( xInds ); }
                        reciY *= -1; if ( reciY < 0 ) { reciY += static_cast< proshade_signed > ( yInds ); }
                        reciZ *= -1; if ( reciZ < 0 ) { reciZ += static_cast< proshade_signed > ( zInds ); }

                        //============================ Apply Friedel's Law
                        arrPos                        = reciZ + static_cast< proshade_signed > ( zInds ) * ( reciY + static_cast< proshade_signed > ( yInds ) * reciX );
                        binIndexing[ static_cast< size_t > ( arrPos ) ] = binIt;
                        
                        //============================ Done, exit bins loop
                        break;
                    }
                }
            }
        }
    }
    
    //================================================ Release memory
    delete[] mins;
    delete[] maxs;
    delete[] resMins;
    delete[] resMinLoc;
    delete[] steps;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the FSC.
 
    This funcion computes the Fourier Shell Correlation weighted average from two identically sized arrays of Fourier coefficients. It requires these to
    have been pre-computed as well as the number of bins and bin mapping to be pre-computed (using the binReciprocalSpaceReflections
    function). Given all these inputs, this function simply computes all the required sums for each bin, processes them and outputs the weighted
    average FSC over all bins.
 
    \param[in] fCoeffs1 The Fourier coefficients of the first map.
    \param[in] fCoeffs2 The Fourier coefficients of the second map.
    \param[in] xInds The number of indices along the x-axis.
    \param[in] yInds The number of indices along the y-axis.
    \param[in] zInds The number of indices along the z-axis.
    \param[in] noBin Number of bins.
    \param[in] binIndexing The map of bin belonging for each reflection.
    \param[in] binData Array of arrays for holding temporary results of the FSC computation. It needs to have been already allocated and have dimensions of noBins x 12. This array is modified by the function in case the caller would like access to these.
    \param[in] binCounts Array of counts for each bin. It needs to be pre-allocated and have dimension of noBins. This array is modified by the function in case the caller would like access to these.
    \param[in] fscByBin This array will hold FSC values for each bin. This is useful in further computations, but could be internal for FSC only computation.
    \param[out] fsc The Fourier Shell Correlation between the two supplied Fourier coefficient maps.
 */
proshade_double ProSHADE_internal_maths::computeFSC ( fftw_complex *fCoeffs1, fftw_complex *fCoeffs2, proshade_unsign xInds, proshade_unsign yInds, proshade_unsign zInds, proshade_signed noBins, proshade_signed* binIndexing, proshade_double**& binData, proshade_signed*& binCounts, proshade_double*& fscByBin )
{
    //================================================ Initialise local variables
    proshade_double realOrig, realRot, imagOrig, imagRot, fsc = 0.0;;
    proshade_signed indx, arrPos;
    std::vector< proshade_double > covarByBin         ( static_cast< size_t > ( noBins ), 0.0 );
    
    //================================================ Clean FSC computation memory
    for ( size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ ) { for ( size_t valIt = 0; valIt < 12; valIt++ ) { binData[binIt][valIt] = 0.0; } }
    for ( size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ ) { binCounts[binIt] = 0; }
    
    //================================================ Compute bin sums
    for ( proshade_signed xIt = 0; xIt < static_cast< proshade_signed > ( xInds ); xIt++ )
    {
        for ( proshade_signed yIt = 0; yIt < static_cast< proshade_signed > ( yInds ); yIt++ )
        {
            for ( proshade_signed zIt = 0; zIt < static_cast< proshade_signed > ( zInds ); zIt++ )
            {
                //==================================== Find array position
                arrPos                                = zIt + static_cast< proshade_signed > ( zInds ) * ( yIt + static_cast< proshade_signed > ( yInds ) * xIt );
                
                //==================================== If no bin is associated, skip this reflection
                indx                                  = binIndexing[ static_cast< size_t > ( arrPos ) ];
                if ( ( indx < 0 ) || ( indx > noBins ) ) { continue; }
                
                //==================================== Calculate the sums
                realOrig                              = fCoeffs1[arrPos][0];
                imagOrig                              = fCoeffs1[arrPos][1];
                realRot                               = fCoeffs2[arrPos][0];
                imagRot                               = fCoeffs2[arrPos][1];
                    
                binData[indx][0]                     += realOrig;
                binData[indx][1]                     += imagOrig;
                binData[indx][2]                     += realRot;
                binData[indx][3]                     += imagRot;
                binData[indx][4]                     += realOrig * realRot;
                binData[indx][5]                     += imagOrig * imagRot;
                binData[indx][6]                     += std::pow ( realOrig, 2.0 );
                binData[indx][7]                     += std::pow ( imagOrig, 2.0 );
                binData[indx][8]                     += std::pow ( realRot,  2.0 );
                binData[indx][9]                     += std::pow ( imagRot,  2.0 );
                
                //==================================== Update bin counts
                binCounts[indx]                      += 1;
            }
        }
    }
    
    //================================================ Compute covariance by bin
    for ( size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ )
    {
        covarByBin.at(binIt)                          = ( ( binData[binIt][4] + binData[binIt][5] ) / static_cast< proshade_double > ( binCounts[binIt] )  -
                                                        ( ( binData[binIt][0]                       / static_cast< proshade_double > ( binCounts[binIt] )   *
                                                            binData[binIt][2]                       / static_cast< proshade_double > ( binCounts[binIt] ) ) +
                                                          ( binData[binIt][1]                       / static_cast< proshade_double > ( binCounts[binIt] )   *
                                                            binData[binIt][3]                       / static_cast< proshade_double > ( binCounts[binIt] ) ) ) );
    }
    
    //================================================ Get FSC by bin
    for ( size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ )
    {
        binData[binIt][10]                            = ( binData[binIt][6] + binData[binIt][7] ) / static_cast< proshade_double > ( binCounts[binIt] ) -
                                                        ( std::pow ( binData[binIt][0] / static_cast< proshade_double > ( binCounts[binIt] ), 2.0 ) +
                                                          std::pow ( binData[binIt][1] / static_cast< proshade_double > ( binCounts[binIt] ), 2.0 ) );
        binData[binIt][11]                            = ( binData[binIt][8] + binData[binIt][9] ) / static_cast< proshade_double > ( binCounts[binIt] ) -
                                                        ( std::pow ( binData[binIt][2] / static_cast< proshade_double > ( binCounts[binIt] ), 2.0 ) +
                                                          std::pow ( binData[binIt][3] / static_cast< proshade_double > ( binCounts[binIt] ), 2.0 ) );
        fscByBin[binIt]                               = covarByBin.at(binIt) / ( std::sqrt ( binData[binIt][10] ) * std::sqrt ( binData[binIt][11] ) );
    }
    
    //================================================ Get average FSC over all bins
    proshade_double binSizeSum                        = 0.0;
    for ( size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ )
    {
        fsc                                          += fscByBin[binIt] * static_cast< proshade_double > ( binCounts[binIt] );
        binSizeSum                                   += static_cast< proshade_double > ( binCounts[binIt] );
    }
    fsc                                              /= static_cast< proshade_double > ( binSizeSum );
    
    //================================================ Done
    return                                            ( fsc );
    
}

/*! \brief This function computes the weights for each reflection using its bin belonging.
 
    This function computes the weights for tralsation optimisation - the bin FSC for each reflection according to its bin belonging for weights1 and
    the square of this value for weights2. Note that weights1 and weights2 are allocated within the function, but since they are the output, the caller
    is expected to delete them after usage.
 
    \param[in] weights1 A pointer to which the FSC by bin weights will be saved into.
    \param[in] weights1 A pointer to which the squared FSC by bin weights will be saved into.
    \param[in] binIndexing The map of bin belonging for each reflection.
    \param[in] fscByBin This array holds FSC values for each bin and is produced by the computeFSC() function.
    \param[in] noBin Number of bins.
    \param[in] xDim The size of the x-dimension of the map in indices.
    \param[in] yDim The size of the y-dimension of the map in indices.
    \param[in] zDim The size of the z-dimension of the map in indices.
 */
void ProSHADE_internal_maths::computeFSCWeightByBin ( proshade_double*& weights1, proshade_double*& weights2, proshade_signed* binIndexing, proshade_double* fscByBin, proshade_signed noBins, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim )
{
    //================================================ Initialise local variables
    proshade_signed indx, arrPos, reciX, reciY, reciZ;
    
    //================================================ Allocate memmory
    weights1                                          = new proshade_double[xDim * yDim * zDim];
    weights2                                          = new proshade_double[xDim * yDim * zDim];
    proshade_single *mins                             = new proshade_single[3];
    proshade_single *maxs                             = new proshade_single[3];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( weights1, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( weights2, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( mins,     __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( maxs,     __FILE__, __LINE__, __func__ );
    
    //================================================ Assign values to the memory
    for ( size_t iter = 0; iter < static_cast< size_t > ( xDim * yDim * zDim ); iter++ ) { weights1[iter] = -100.0; weights2[iter] = -100.0; }
    
    //================================================ Determine reciprocal space indexing ranges
    mins[0]                                           = std::floor ( static_cast< proshade_single > ( xDim ) / -2.0f );
    mins[1]                                           = std::floor ( static_cast< proshade_single > ( yDim ) / -2.0f );
    mins[2]                                           = std::floor ( static_cast< proshade_single > ( zDim ) / -2.0f );
        
    maxs[0]                                           = -mins[0];
    maxs[1]                                           = -mins[1];
    maxs[2]                                           = -mins[2];
    
    if ( xDim % 2 == 0 )                              { mins[0] += 1.0f; }
    if ( yDim % 2 == 0 )                              { mins[1] += 1.0f; }
    if ( zDim % 2 == 0 )                              { mins[2] += 1.0f; }
    
    //================================================ For each reflection
    for ( proshade_signed xIt = 0; xIt < xDim; xIt++ )
    {
        for ( proshade_signed yIt = 0; yIt < yDim; yIt++ )
        {
            for ( proshade_signed zIt = 0; zIt < ( ( zDim / 2 ) + 1 ); zIt++ )
            {
                //==================================== Deal with reciprocal indices ordering
                reciX = xIt; if ( reciX > static_cast< proshade_signed > ( maxs[0] ) ) { reciX -= static_cast< proshade_signed > ( xDim ); }
                reciY = yIt; if ( reciY > static_cast< proshade_signed > ( maxs[1] ) ) { reciY -= static_cast< proshade_signed > ( yDim ); }
                reciZ = zIt; if ( reciZ > static_cast< proshade_signed > ( maxs[2] ) ) { reciZ -= static_cast< proshade_signed > ( zDim ); }
                
                //==================================== Find array position and bin index
                arrPos                                = zIt + zDim * ( yIt + yDim * xIt );
                indx                                  = binIndexing[ static_cast< size_t > ( arrPos ) ];

                //==================================== Ignore unassigned bin reflections
                if ( ( indx < 0 ) || ( indx > noBins ) ) { continue; }

                //==================================== Set weights using the bin FSC
                weights1[arrPos]                      = fscByBin[indx];
                weights2[arrPos]                      = std::pow ( fscByBin[indx], 2.0 );
                
                //==================================== If one of the uneven ends, do not use Friedel's Law
                if ( reciX == static_cast< proshade_signed > ( mins[0] ) || -reciX == static_cast< proshade_signed > ( mins[0] ) ) { break; }
                if ( reciY == static_cast< proshade_signed > ( mins[1] ) || -reciY == static_cast< proshade_signed > ( mins[1] ) ) { break; }
                if ( reciZ == static_cast< proshade_signed > ( mins[2] ) || -reciZ == static_cast< proshade_signed > ( mins[2] ) ) { break; }
                
                //==================================== Use Friedel's Law to find the second index (this is why we can use zDim / 2)
                reciX *= -1; if ( reciX < 0 ) { reciX += xDim; }
                reciY *= -1; if ( reciY < 0 ) { reciY += yDim; }
                reciZ *= -1; if ( reciZ < 0 ) { reciZ += zDim; }

                //==================================== Apply Friedel's Law
                arrPos                                = reciZ + zDim * ( reciY + yDim * reciX );
                weights1[arrPos]                      = fscByBin[indx];
                weights2[arrPos]                      = std::pow ( fscByBin[indx], 2.0 );
            }
        }
    }
    
    //================================================ Release memory
    delete[] mins;
    delete[] maxs;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the real part of the sum of all coefficients except where the weight is less than -2.
 
    \param[in] fCoeffs The Fourier coefficients to be moved.
    \param[in] weights The weights to be applied to the shift.
    \param[in] xDim The size of the x-dimension of the map in indices.
    \param[in] yDim The size of the y-dimension of the map in indices.
    \param[in] zDim The size of the z-dimension of the map in indices.
    \param[out] sum The F value.
 */
proshade_double ProSHADE_internal_maths::computeTheFValue ( proshade_complex* fCoeffs, proshade_double* weights, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim )
{
    //================================================ Initialise local variables
    proshade_double sum                               = 0.0;
    proshade_signed arrPos;
    
    //================================================ For each reflection
    for ( proshade_signed xIt = 0; xIt < xDim; xIt++ )
    {
        for ( proshade_signed yIt = 0; yIt < yDim; yIt++ )
        {
            for ( proshade_signed zIt = 0; zIt < ( ( zDim / 2 ) + 1 ); zIt++ )
            {
                //==================================== Find array position and bin index
                arrPos                                = zIt + zDim * ( yIt + yDim * xIt );
                
                //==================================== Sum real parts if weight exists
                if ( weights[arrPos] > -2.0 ) { sum += fCoeffs[arrPos][0]; }
            }
        }
    }
    
    //================================================ Done
    return                                            ( sum );
    
}

/*! \brief This function computes the first and second derivatives of the translation function at coefficient [0,0,0].
 
    \param[in] fCoeffs The Fourier coefficients to be used.
    \param[in] weights1 The weights to be applied to the computation of the first derivatives (given by function computeFSCWeightByBin() ).
    \param[in] weights2 The weights to be applied to the computation of the second derivatives (given by function computeFSCWeightByBin() ).
    \param[in] xDim The size of the x-dimension of the map in indices.
    \param[in] yDim The size of the y-dimension of the map in indices.
    \param[in] zDim The size of the z-dimension of the map in indices.
    \param[in] firstDers Pointer to array where the first derivatives (array of 3) will be stored. This function will allocate the memory, but the caller will have to delete it.
    \param[in] secondDers Pointer to array where the second derivatives (array of 9) will be stored. This function will allocate the memory, but the caller will have to delete it.
 */
void ProSHADE_internal_maths::computeTrFunDerivatives ( proshade_complex* fCoeffs, proshade_double* weights1, proshade_double* weights2, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim, proshade_double*& firstDers, proshade_double*& secondDers )
{
    //================================================ Allocate memmory
    firstDers                                         = new proshade_double[3];
    secondDers                                        = new proshade_double[9];
    proshade_single *mins                             = new proshade_single[3];
    proshade_single *maxs                             = new proshade_single[3];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( firstDers,  __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( secondDers, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( mins,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( maxs,       __FILE__, __LINE__, __func__ );
    
    //================================================ Initialise variables
    std::complex< proshade_double > piConstFirst      ( 2.0 * M_PI, 1.0 );
    proshade_double piConstSecond                     = std::pow ( 2.0 * M_PI, 2.0 );
    for ( size_t iter = 0; iter < 3; iter++ ) { firstDers[iter]  = 0.0; }
    for ( size_t iter = 0; iter < 9; iter++ ) { secondDers[iter] = 0.0; }
    proshade_signed reciX, reciY, reciZ, arrPos;
    
    //================================================ Determine reciprocal space indexing ranges
    mins[0]                                           = std::floor ( static_cast< proshade_single > ( xDim ) / -2.0f );
    mins[1]                                           = std::floor ( static_cast< proshade_single > ( yDim ) / -2.0f );
    mins[2]                                           = std::floor ( static_cast< proshade_single > ( zDim ) / -2.0f );
        
    maxs[0]                                           = -mins[0];
    maxs[1]                                           = -mins[1];
    maxs[2]                                           = -mins[2];
    
    if ( xDim % 2 == 0 )                              { mins[0] += 1.0f; }
    if ( yDim % 2 == 0 )                              { mins[1] += 1.0f; }
    if ( zDim % 2 == 0 )                              { mins[2] += 1.0f; }
    
    //================================================ For each reflection
    for ( proshade_signed xIt = 0; xIt < xDim; xIt++ )
    {
        for ( proshade_signed yIt = 0; yIt < yDim; yIt++ )
        {
            for ( proshade_signed zIt = 0; zIt < ( ( zDim / 2 ) + 1 ); zIt++ )
            {
                //==================================== Deal with reciprocal indices ordering
                reciX = xIt; if ( reciX > static_cast< proshade_signed > ( maxs[0] ) ) { reciX -= static_cast< proshade_signed > ( xDim ); }
                reciY = yIt; if ( reciY > static_cast< proshade_signed > ( maxs[1] ) ) { reciY -= static_cast< proshade_signed > ( yDim ); }
                reciZ = zIt; if ( reciZ > static_cast< proshade_signed > ( maxs[2] ) ) { reciZ -= static_cast< proshade_signed > ( zDim ); }
                
                //==================================== Find array position and bin index
                arrPos                                = zIt + zDim * ( yIt + yDim * xIt );
                
                //==================================== Ignore if outside of weights
                if ( weights1[arrPos] < -2.0 ) { continue; }
                
                //==================================== Add to the first derivatives sum
                firstDers[0] += ( weights1[arrPos] * fCoeffs[arrPos][0] * std::conj( fCoeffs[arrPos][1] * piConstFirst * static_cast< proshade_double > ( reciX ) ) ).real();
                firstDers[1] += ( weights1[arrPos] * fCoeffs[arrPos][0] * std::conj( fCoeffs[arrPos][1] * piConstFirst * static_cast< proshade_double > ( reciY ) ) ).real();
                firstDers[2] += ( weights1[arrPos] * fCoeffs[arrPos][0] * std::conj( fCoeffs[arrPos][1] * piConstFirst * static_cast< proshade_double > ( reciZ ) ) ).real();
                
                //==================================== Add to the second derivatives sum
                secondDers[0]                        += weights2[arrPos] * reciX * reciX;
                secondDers[1]                        += weights2[arrPos] * reciX * reciY;
                secondDers[2]                        += weights2[arrPos] * reciX * reciZ;
                secondDers[4]                        += weights2[arrPos] * reciY * reciY;
                secondDers[5]                        += weights2[arrPos] * reciY * reciZ;
                secondDers[8]                        += weights2[arrPos] * reciZ * reciZ;
            }
        }
    }
    
    //================================================ Complete second darivatives matrix
    secondDers[3]                                     = secondDers[1];
    secondDers[6]                                     = secondDers[2];
    secondDers[7]                                     = secondDers[5];
    for ( size_t iter = 0; iter < 9; iter++ ) { secondDers[iter] *= -piConstSecond; }
    
    //================================================ Release memory
    delete[] mins;
    delete[] maxs;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the step sizes for translation function optimisation from the first and second derivatives.
 
    \param[in] firstDers Pointer to array where the first derivatives are stored (as computed by computeTrFunDerivatives() ).
    \param[in] secondDers Pointer to array where the second derivatives are stored (as computed by computeTrFunDerivatives() ).
    \param[out] stepArr An array holding the step sizes along the three dimensions. It is allocated here, but the caller is required to delete the pointer.
 */
proshade_double* ProSHADE_internal_maths::computeTrFunStep ( proshade_double* firstDers, proshade_double* secondDers )
{
    //================================================ Change format of second derivatives and add I matrix (the inversion function will subtract it)
    std::vector < proshade_double > tmpMap            ( 9, 0.0 );
    for ( size_t iter = 0; iter < 9; iter++ ) { tmpMap.at(iter) = secondDers[iter]; }
    tmpMap.at(0) += 1.0; tmpMap.at(4) += 1.0; tmpMap.at(8) += 1.0;
    
    //================================================ Compute matrix inversion for the second derivatives matrix
    proshade_double* secondInv                        = compute3x3MoorePenrosePseudoInverseOfIMinusMat ( &tmpMap, -1 );
    
    //================================================ Compute dot product between the inverted matrix and the first derivatives vector
    proshade_double* stepArr                          = compute3x3MatrixVectorMultiplication ( secondInv, -firstDers[0], -firstDers[1], -firstDers[2] );
    
    //================================================ Release memory
    delete[] secondInv;
    
    //================================================ Done
    return                                            ( stepArr );
    
}

/*! \brief This function simply finds all the peaks in a 1D data array.
 
    Simple function for detecting all points in the supplied array which have higher value than all surrounding points along the single dimension of the array.
 
    \param[in] data The input array containning (pressumably smoothened) data.
    \param[out] peaks A vector containing all the peak indices in the  input array.
 */
std::vector< proshade_signed > ProSHADE_internal_maths::findPeaks1D ( std::vector< proshade_double > data )
{
    //================================================ Initialise local variables
    std::vector< proshade_signed > ret;
    
    //================================================ Peak is simply any position with both neighbours having lower position (with special care for borders)
    for ( proshade_signed index = 0; index < static_cast< proshade_signed > ( data.size() ); index++ )
    {
        //============================================ Starting border?
        if ( index == 0 )
        {
            if ( data.size() > 1 ) { if ( data.at(0) > data.at(1) ) { ProSHADE_internal_misc::addToSignedVector ( &ret, index ); } }
            continue;
        }
        
        //============================================ End border?
        if ( index == static_cast< proshade_signed > ( data.size() - 1 ) )
        {
            if ( data.at( static_cast< size_t > ( index ) ) > data.at( static_cast< size_t > ( index ) - 1 ) ) { ProSHADE_internal_misc::addToSignedVector ( &ret, index ); }
            continue;
        }
        
        //============================================ Is this a peak?
        if ( ( data.at( static_cast< size_t > ( index ) ) > data.at( static_cast< size_t > ( index ) - 1 ) ) &&
             ( data.at( static_cast< size_t > ( index ) ) > data.at( static_cast< size_t > ( index ) + 1 ) ) ) { ProSHADE_internal_misc::addToSignedVector ( &ret, index ); }
        
        //============================================ Deal with equally sized values
        if ( index < static_cast< proshade_signed > ( data.size() - 2 ) ) { if ( data.at( static_cast< size_t > ( index ) ) >= data.at( static_cast< size_t > ( index ) - 1 ) ) { if ( data.at( static_cast< size_t > ( index ) ) >= data.at( static_cast< size_t > ( index ) + 1 ) ) { if ( data.at( static_cast< size_t > ( index ) ) > data.at( static_cast< size_t > ( index ) + 2 ) ) { ProSHADE_internal_misc::addToSignedVector ( &ret, index ); } } } }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}


/*! \brief This function finds a subgroup of axes with distinctly higher correlation value.
 
    This function starts by getting a vector of all peak heights detected in all C symmetries detected in the structure. It then proceeds to convert
    this vector into a histogram, which it then smoothens using Gaussian convolution according to the input parameters. Finally, it searches for
    peaks in the smoothened histogram function and reports the minimal value that average peak height must be in order for the axis to be
    considered part of this top group.
 
    \param[in] CSym A vector of pointers to double arrays, each array being a single Cyclic symmetry entry.
    \param[in] step The granulosity of the interval <0,1> using which the search should be done.
    \param[in] sigma The variance of the Gaussian used to smoothen the peak height histogram.
    \param[in] windowSize The width of the window over which smoothening is done.
    \param[in] maxLim The maximum value that can be reached - this is to step a single extremely high peak overshadowing very good peaks.
    \param[out] threshold The minimum peak height that an axis needs to have to be considered a member of the distinct top group.
 */
proshade_double ProSHADE_internal_maths::findTopGroupSmooth ( std::vector< proshade_double* >* CSym, size_t peakPos, proshade_double step, proshade_double sigma, proshade_signed windowSize, proshade_double maxLim )
{
    //================================================ Initialise local variables
    proshade_double threshold                         = 0.0;
    proshade_signed totSize                           = static_cast< proshade_signed > ( ( 1.0 / step ) + 1 );
    std::vector< std::pair < proshade_double, proshade_unsign > > vals;
    std::vector< proshade_double > hist               ( static_cast< unsigned long int > ( totSize ), 0.0 );
    proshade_unsign histPos                           = 0;
    
    //================================================ Make sure window size is odd
    if ( windowSize % 2 == 0 )                        { windowSize += 1; }
    
    //================================================ Get vector of pairs of peak heights and indices in CSym array
    for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( CSym->size() ); symIt++ ) { if ( !std::isinf( CSym->at(symIt)[peakPos] ) ) { vals.emplace_back ( std::pair < proshade_double, proshade_unsign > ( CSym->at(symIt)[peakPos], symIt ) ); } }
    
    //================================================ Bump all top peaks together - we do not want single high peak overshadowing many good peaks
    for ( proshade_unsign vIt = 0; vIt < static_cast< proshade_unsign > ( vals.size() ); vIt++ ) { if ( vals.at(vIt).first > maxLim ) { vals.at(vIt).first = maxLim; } }
    
    //================================================ Convert all found heights to histogram from 0.0 to 1.0 by step
    for ( proshade_double it = 0.0; it <= 1.0; it = it + step )
    {
        for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( vals.size() ); symIt++ )
        {
            //======================================== Is this height in the range?
            if ( ( vals.at(symIt).first > it ) && ( vals.at(symIt).first <= ( it + step ) ) ) { hist.at(histPos) += 1.0; }
        }
        
        //============================================ Update counter and continue
        histPos                                      += 1;
    }
    
    //================================================ Smoothen the distribution
    std::vector< proshade_double > smoothened         = ProSHADE_internal_maths::smoothen1D ( step, windowSize, sigma, hist );
    
    //================================================ Find peaks in smoothened data
    std::vector< proshade_signed > peaks              = ProSHADE_internal_maths::findPeaks1D ( smoothened );
    
    //================================================ Take best peaks surroundings and produce a new set of "high" axes
    proshade_signed bestHistPos;
    if ( peaks.size() > 0 ) { bestHistPos = peaks.at(peaks.size()-1) + ( ( windowSize - 1 ) / 2 ); }
    else                    { bestHistPos = 0.0; }
    
    threshold                                         = ( static_cast< proshade_double > ( bestHistPos ) * step );
    
    //================================================ Done
    return                                            ( threshold );
    
}

/*! \brief This function finds a subgroup of axes with distinctly higher correlation value.
 
    This function starts by getting a vector of all peak heights detected in all C symmetries detected in the structure. It then proceeds to convert
    this vector into a histogram, which it then smoothens using Gaussian convolution according to the input parameters. Finally, it searches for
    peaks in the smoothened histogram function and reports the minimal value that average peak height must be in order for the axis to be
    considered part of this top group.
 
    \param[in] CSym A vector of vectors of doubles, each array being a single Cyclic symmetry entry.
    \param[in] step The granulosity of the interval <0,1> using which the search should be done.
    \param[in] sigma The variance of the Gaussian used to smoothen the peak height histogram.
    \param[in] windowSize The width of the window over which smoothening is done.
    \param[in] maxLim The maximum value that can be reached - this is to step a single extremely high peak overshadowing very good peaks.
    \param[out] threshold The minimum peak height that an axis needs to have to be considered a member of the distinct top group.
 */
proshade_double ProSHADE_internal_maths::findTopGroupSmooth ( std::vector< std::vector< proshade_double > >* CSym, size_t peakPos, proshade_double step, proshade_double sigma, proshade_signed windowSize, proshade_double maxLim )
{
    //================================================ Initialise local variables
    proshade_double threshold                         = 0.0;
    proshade_signed totSize                           = static_cast< proshade_signed > ( ( 1.0 / step ) + 1 );
    std::vector< std::pair < proshade_double, proshade_unsign > > vals;
    std::vector< proshade_double > hist               ( static_cast< unsigned long int > ( totSize ), 0.0 );
    proshade_unsign histPos                           = 0;
    
    //================================================ Make sure window size is odd
    if ( windowSize % 2 == 0 )                        { windowSize += 1; }
    
    //================================================ Get vector of pairs of peak heights and indices in CSym array
    for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( CSym->size() ); symIt++ ) { if ( !std::isinf( CSym->at(symIt).at(peakPos) ) ) { vals.emplace_back ( std::pair < proshade_double, proshade_unsign > ( CSym->at(symIt).at(peakPos), symIt ) ); } }
    
    //================================================ Bump all top peaks together - we do not want single high peak overshadowing many good peaks
    for ( proshade_unsign vIt = 0; vIt < static_cast< proshade_unsign > ( vals.size() ); vIt++ ) { if ( vals.at(vIt).first > maxLim ) { vals.at(vIt).first = maxLim; } }
    
    
    //================================================ Convert all found heights to histogram from 0.0 to 1.0 by step
    for ( proshade_double it = 0.0; it <= 1.0; it = it + step )
    {
        for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( vals.size() ); symIt++ )
        {
            //======================================== Is this height in the range?
            if ( ( vals.at(symIt).first > it ) && ( vals.at(symIt).first <= ( it + step ) ) ) { hist.at(histPos) += 1.0; }
        }
        
        //============================================ Update counter and continue
        histPos                                      += 1;
    }
    
    //================================================ Smoothen the distribution
    std::vector< proshade_double > smoothened         = ProSHADE_internal_maths::smoothen1D ( step, windowSize, sigma, hist );
    
    //================================================ Find peaks in smoothened data
    std::vector< proshade_signed > peaks              = ProSHADE_internal_maths::findPeaks1D ( smoothened );
    
    //================================================ Take best peaks surroundings and produce a new set of "high" axes
    proshade_signed bestHistPos;
    if ( peaks.size() > 0 ) { bestHistPos = peaks.at(peaks.size()-1) + ( ( windowSize - 1 ) / 2 ); }
    else                    { bestHistPos = 0.0; }
    
    threshold                                         = ( static_cast< proshade_double > ( bestHistPos ) * step );
    
    //================================================ Done
    return                                            ( threshold );
    
}

/*! \brief This function combines Fourier coefficients of two structures in a way, so that inverse Fourier of the combination will be the translation function.
 
    \param[in] tmpOut1 Array holding the static structure Fourier outputs.
    \param[in] tmpOut2 Array holding the moving structure Fourier outputs.
    \param[in] resOut Array to hold the combined Fourier coefficients of both structures.
    \param[in] xD The dimension of the X axis of the structures (assumes both structures have the same sizes and sampling).
    \param[in] yD The dimension of the Y axis of the structures (assumes both structures have the same sizes and sampling).
    \param[in] zD The dimension of the Z axis of the structures (assumes both structures have the same sizes and sampling).
 */
void ProSHADE_internal_maths::combineFourierForTranslation ( fftw_complex* tmpOut1, fftw_complex* tmpOut2, fftw_complex*& resOut, proshade_unsign xD, proshade_unsign yD, proshade_unsign zD )
{
    //================================================ Initialise local variables
    double normFactor                                 = static_cast<double> ( xD * yD * zD );
    proshade_signed arrPos;
    
    //================================================ Combine the coefficients
    for ( proshade_signed xIt = 0; xIt < static_cast< proshade_signed > ( xD ); xIt++ )
    {
        for ( proshade_signed yIt = 0; yIt < static_cast< proshade_signed > ( yD ); yIt++ )
        {
            for ( proshade_signed zIt = 0; zIt < static_cast< proshade_signed > ( zD ); zIt++ )
            {
                //==================================== Find indices
                arrPos                                = zIt   + static_cast< proshade_signed > ( zD ) * ( yIt   + static_cast< proshade_signed > ( yD ) * xIt );
                
                //==================================== Combine
                ProSHADE_internal_maths::complexMultiplicationConjug ( &tmpOut1[arrPos][0],
                                                                       &tmpOut1[arrPos][1],
                                                                       &tmpOut2[arrPos][0],
                                                                       &tmpOut2[arrPos][1],
                                                                       &resOut[arrPos][0],
                                                                       &resOut[arrPos][1] );
                
                //==================================== Save
                resOut[arrPos][0]                    /= normFactor;
                resOut[arrPos][1]                    /= normFactor;
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function simply finds the highest value in fftw_complex map and returns its position and value.
 
    \param[in] resIn Array holding the translation function values.
    \param[in] xD The dimension of the X axis of the structures (assumes both structures have the same sizes and sampling).
    \param[in] yD The dimension of the Y axis of the structures (assumes both structures have the same sizes and sampling).
    \param[in] zD The dimension of the Z axis of the structures (assumes both structures have the same sizes and sampling).
    \param[in] trsX Variable to which the X axis translation function peak position will be saved to.
    \param[in] trsY Variable to which the Y axis translation function peak position will be saved to.
    \param[in] trsZ Variable to which the Z axis translation function peak position will be saved to.
    \param[in] mapPeak Variable to which the height of the translation function peak will be saved to.
 */
void ProSHADE_internal_maths::findHighestValueInMap ( fftw_complex* resIn, proshade_unsign xD, proshade_unsign yD, proshade_unsign zD, proshade_double* trsX, proshade_double* trsY, proshade_double* trsZ, proshade_double* mapPeak )
{
    //================================================ Initialise variables
    proshade_signed arrPos;
   *mapPeak                                           = 0.0;
    
    //================================================ Search the map
    for ( proshade_signed uIt = 0; uIt < static_cast<proshade_signed> ( xD ); uIt++ )
    {
        for ( proshade_signed vIt = 0; vIt < static_cast<proshade_signed> ( yD ); vIt++ )
        {
            for ( proshade_signed wIt = 0; wIt < static_cast<proshade_signed> ( zD ); wIt++ )
            {
                arrPos                                = wIt + static_cast< proshade_signed > ( zD ) * ( vIt + static_cast< proshade_signed > ( yD ) * uIt );
                if ( resIn[arrPos][0] > *mapPeak )
                {
                   *mapPeak                           = resIn[arrPos][0];
                   *trsX                              = static_cast< proshade_double > ( uIt );
                   *trsY                              = static_cast< proshade_double > ( vIt );
                   *trsZ                              = static_cast< proshade_double > ( wIt );
                }
            }
        }
    }
    
    //================================================ Done
    return ;
    
}
