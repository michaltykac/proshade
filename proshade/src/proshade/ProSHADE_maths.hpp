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
    \version   0.7.5.4
    \date      MAR 2021
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

//==================================================== ProSHADE_internal_spheres Namespace
// ... NOTE: The full namespace is defined in ProSHADE_spheres.h ; however, as forward declaration
// ... would not work here, the ProSHADE_rotFun_sphere class needs to be defined here, making the
// ... organisation a bit messier.
namespace ProSHADE_internal_spheres
{
/*! \class ProSHADE_rotFun_sphere
    \brief This class contains all inputed data for the rotation function angle-axis converted spheres.
 
    This class codes the object that contains all the information about a single concentric sphere in the angle-axis space
    obtained by conversion of the rotation function. It also contains all the functionality required
    to process such data.
 */
    class ProSHADE_rotFun_sphere
    {
    private:
        proshade_double radius;
        proshade_double radiusMax;
        proshade_double radiusMin;
        proshade_unsign angularDim;
        proshade_double representedAngle;
        proshade_unsign sphereNumber;
        
        proshade_double* axesValues;
        std::vector< std::pair< proshade_unsign,proshade_unsign > > peaks;
    public:
        ProSHADE_rotFun_sphere                        ( proshade_double rad, proshade_double radRange, proshade_unsign dim, proshade_double repAng, proshade_unsign sphNo );
       ~ProSHADE_rotFun_sphere                        ( void );
        
    public:
        proshade_double getRadius                     ( void );
        proshade_double getMaxRadius                  ( void );
        proshade_double getMinRadius                  ( void );
        proshade_unsign getAngularDim                 ( void );
        proshade_double getRepresentedAngle           ( void );
        proshade_unsign getSphereNumber               ( void );
        std::vector<std::pair<proshade_unsign,proshade_unsign>> getPeaks ( void );
        proshade_double getSphereLatLonPosition       ( proshade_unsign lattitude, proshade_unsign longitude );
        proshade_double getSphereLatLonLinearInterpolationPos ( proshade_double lattitude, proshade_double longitude );
        
        
    public:
        void interpolateSphereValues                  ( proshade_complex* rotFun );
        void findAllPeaks                             ( proshade_signed noSmNeighbours, std::vector< proshade_double >* allHeights );
        void removeSmallPeaks                         ( proshade_double peakThres );
    };
}

//==================================================== ProSHADE_internal_io Namespace
/*! \namespace ProSHADE_internal_maths
    \brief This namespace contains the internal functions for common mathematical operations.
 
    The ProSHADE_internal_maths namespace contains a set of common mathematical operations used in many places by ProSHADE. These
    typically include complex number operations and angle conversions.
 */
namespace ProSHADE_internal_maths
{
/*! \class ProSHADE_internal_maths_bicubicInterpolator
    \brief This class takes 4 by 4 array of input point values and allows bicubic interpolation of any value between them.
 
    This class provides a constructor, which takes a 4 by 4 matrix of values. These values signify a function f(x,y) values at equidistantly
    placed indices on a unit square. Once the constructor is finished, the function getValue() an be used to obtain the interpolated value
    of any position on this unit square.
 */
    class BicubicInterpolator
    {
    private:
        proshade_double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
        proshade_double xStartIndex, xRange, yStartIndex, yRange;

    public:
/*!     \brief This is the constructor for the BicubicInterpolator class.

        This constructor takes the 4 by 4 array of values and pre-computes all the interpolators required for fast computation of any
        positions value.

        \param[in] areaToInterpolate Pointer to pointer of 4 by 4 equidistantly spaced grid of values.
 */
        BicubicInterpolator ( proshade_double** areaToInterpolate, proshade_double xStart, proshade_double yStart )
        {
            //======================================== Save the original non-unit square positions
            this->xStartIndex                         = xStart;
            this->yStartIndex                         = yStart;
            
            //======================================== Prepare variables for converting from original indices to unit square
            this->xRange                              = 1.0;
            this->yRange                              = 1.0;
            
            //======================================== Precompute interpolators
            this->a00                                 = areaToInterpolate[1][1];
            this->a01                                 = - ( 0.5 * areaToInterpolate[1][0] ) +
                                                          ( 0.5 * areaToInterpolate[1][2] );
            this->a02                                 =   areaToInterpolate[1][0] -
                                                          ( 2.5 * areaToInterpolate[1][1] ) +
                                                          ( 2.0 * areaToInterpolate[1][2] ) -
                                                          ( 0.5 * areaToInterpolate[1][3] );
            this->a03                                 = - ( 0.5 * areaToInterpolate[1][0] ) +
                                                          ( 1.5 * areaToInterpolate[1][1] ) -
                                                          ( 1.5 * areaToInterpolate[1][2] ) +
                                                          ( 0.5 * areaToInterpolate[1][3] );
            this->a10                                 = - ( 0.5 * areaToInterpolate[0][1] ) +
                                                          ( 0.5 * areaToInterpolate[2][1] );
            this->a11                                 =   ( 0.25 * areaToInterpolate[0][0] ) -
                                                          ( 0.25 * areaToInterpolate[0][2] ) -
                                                          ( 0.25 * areaToInterpolate[2][0] ) +
                                                          ( 0.25 * areaToInterpolate[2][2] );
            this->a12                                 = - ( 0.5 * areaToInterpolate[0][0] ) +
                                                          ( 1.25 * areaToInterpolate[0][1] ) -
                                                          areaToInterpolate[0][2] +
                                                          ( 0.25 * areaToInterpolate[0][3] ) +
                                                          ( 0.5 * areaToInterpolate[2][0] ) -
                                                          ( 1.25 * areaToInterpolate[2][1] ) +
                                                          areaToInterpolate[2][2] -
                                                          ( 0.25 * areaToInterpolate[2][3] );
            this->a13                                 =   ( 0.25 * areaToInterpolate[0][0] ) -
                                                          ( 0.75 * areaToInterpolate[0][1] ) +
                                                          ( 0.75 * areaToInterpolate[0][2] ) -
                                                          ( 0.25 * areaToInterpolate[0][3] ) -
                                                          ( 0.25 * areaToInterpolate[2][0] ) +
                                                          ( 0.75 * areaToInterpolate[2][1] ) -
                                                          ( 0.75 * areaToInterpolate[2][2] ) +
                                                          ( 0.25 * areaToInterpolate[2][3] );
            this->a20                                 =   areaToInterpolate[0][1] -
                                                          ( 2.5 * areaToInterpolate[1][1] ) +
                                                          ( 2.0 * areaToInterpolate[2][1] ) -
                                                          ( 0.5 * areaToInterpolate[3][1] );
            this->a21                                 = - ( 0.5 * areaToInterpolate[0][0] ) +
                                                          ( 0.5 * areaToInterpolate[0][2] ) +
                                                          ( 1.25 * areaToInterpolate[1][0] ) -
                                                          ( 1.25 * areaToInterpolate[1][2] ) -
                                                          areaToInterpolate[2][0] + areaToInterpolate[2][2] +
                                                          ( 0.25 * areaToInterpolate[3][0] ) -
                                                          ( 0.25 * areaToInterpolate[3][2] );
            this->a22                                 =   areaToInterpolate[0][0] -
                                                          ( 2.5 * areaToInterpolate[0][1] ) +
                                                          ( 2.0 * areaToInterpolate[0][2] ) -
                                                          ( 0.5 * areaToInterpolate[0][3] ) -
                                                          ( 2.5 * areaToInterpolate[1][0] ) +
                                                          ( 6.25 * areaToInterpolate[1][1] ) -
                                                          ( 5.0 * areaToInterpolate[1][2] ) +
                                                          ( 1.25 * areaToInterpolate[1][3] ) +
                                                          ( 2.0 * areaToInterpolate[2][0] ) -
                                                          ( 5.0 * areaToInterpolate[2][1] ) +
                                                          ( 4.0 * areaToInterpolate[2][2] ) -
                                                          areaToInterpolate[2][3] -
                                                          ( 0.5 * areaToInterpolate[3][0] ) +
                                                          ( 1.25 * areaToInterpolate[3][1] ) -
                                                          areaToInterpolate[3][2] +
                                                          ( 0.25 * areaToInterpolate[3][3] );
            this->a23                                 = - ( 0.5 * areaToInterpolate[0][0] ) +
                                                          ( 1.5 * areaToInterpolate[0][1] ) -
                                                          ( 1.5 * areaToInterpolate[0][2] ) +
                                                          ( 0.5 * areaToInterpolate[0][3] ) +
                                                          ( 1.25 * areaToInterpolate[1][0] ) -
                                                          ( 3.75 * areaToInterpolate[1][1] ) +
                                                          ( 3.75 * areaToInterpolate[1][2] ) -
                                                          ( 1.25 * areaToInterpolate[1][3] ) -
                                                          areaToInterpolate[2][0] +
                                                          ( 3.0 * areaToInterpolate[2][1] ) -
                                                          ( 3.0 * areaToInterpolate[2][2] ) +
                                                          areaToInterpolate[2][3] +
                                                          ( 0.25 * areaToInterpolate[3][0] ) -
                                                          ( 0.75 * areaToInterpolate[3][1] ) +
                                                          ( 0.75 * areaToInterpolate[3][2] ) -
                                                          ( 0.25 * areaToInterpolate[3][3] );
            this->a30                                 = - ( 0.5 * areaToInterpolate[0][1] ) +
                                                          ( 1.5 * areaToInterpolate[1][1] ) -
                                                          ( 1.5 * areaToInterpolate[2][1] ) +
                                                          ( 0.5*areaToInterpolate[3][1] );
            this->a31                                 =   ( 0.25 * areaToInterpolate[0][0] ) -
                                                          ( 0.25 * areaToInterpolate[0][2] ) -
                                                          ( 0.75 * areaToInterpolate[1][0] ) +
                                                          ( 0.75 * areaToInterpolate[1][2] ) +
                                                          ( 0.75 * areaToInterpolate[2][0] ) -
                                                          ( 0.75 * areaToInterpolate[2][2] ) -
                                                          ( 0.25 * areaToInterpolate[3][0] ) +
                                                          ( 0.25 * areaToInterpolate[3][2] );
            this->a32                                 = - ( 0.5 * areaToInterpolate[0][0] ) +
                                                          ( 1.25 * areaToInterpolate[0][1] ) -
                                                          areaToInterpolate[0][2] +
                                                          ( 0.25 * areaToInterpolate[0][3] ) +
                                                          ( 1.5 * areaToInterpolate[1][0] ) -
                                                          ( 3.75 * areaToInterpolate[1][1] ) +
                                                          ( 3.0 * areaToInterpolate[1][2] ) -
                                                          ( 0.75 * areaToInterpolate[1][3] ) -
                                                          ( 1.5 * areaToInterpolate[2][0] ) +
                                                          ( 3.75 * areaToInterpolate[2][1] ) -
                                                          ( 3.0 * areaToInterpolate[2][2] ) +
                                                          ( 0.75 * areaToInterpolate[2][3] ) +
                                                          ( 0.5 * areaToInterpolate[3][0] ) -
                                                          ( 1.25 * areaToInterpolate[3][1] ) +
                                                          areaToInterpolate[3][2] -
                                                          ( 0.25 * areaToInterpolate[3][3] );
            this->a33                                 =   ( 0.25 * areaToInterpolate[0][0] ) -
                                                          ( 0.75 * areaToInterpolate[0][1] ) +
                                                          ( 0.75 * areaToInterpolate[0][2] ) -
                                                          ( 0.25 * areaToInterpolate[0][3] ) -
                                                          ( 0.75 * areaToInterpolate[1][0] ) +
                                                          ( 2.25 * areaToInterpolate[1][1] ) -
                                                          ( 2.25 * areaToInterpolate[1][2] ) +
                                                          ( 0.75 * areaToInterpolate[1][3] ) +
                                                          ( 0.75 * areaToInterpolate[2][0] ) -
                                                          ( 2.25 * areaToInterpolate[2][1] ) +
                                                          ( 2.25 * areaToInterpolate[2][2] ) -
                                                          ( 0.75 * areaToInterpolate[2][3] ) -
                                                          ( 0.25 * areaToInterpolate[3][0] ) +
                                                          ( 0.75 * areaToInterpolate[3][1] ) -
                                                          ( 0.75 * areaToInterpolate[3][2] ) +
                                                          ( 0.25 * areaToInterpolate[3][3] );
        }
        
/*!     \brief This is the destructor for the BicubicInterpolator class.

        This destructor does nothing.
 */
       ~BicubicInterpolator ( void ) { ; }

/*!     \brief This function allows accessing the interpolated value for a given position x and y.

        This constructor takes the 4 by 4 array of values and pre-computes all the interpolators required for fast computation of any
        positions value.

        \param[in] x The x-axis position on the unit square for which the interpolated value should be computed.
        \param[in] y The y-axis position on the unit square for which the interpolated value should be computed.
        \param[out] res The interpolated value for the position x and y.
 */
        proshade_double getValue ( proshade_double x, proshade_double y )
        {
            //======================================== Sanity check
            if ( ( ( x < this->xStartIndex ) || ( x > ( this->xStartIndex + this->xRange ) ) ) ||
                 ( ( y < this->yStartIndex ) || ( y > ( this->yStartIndex + this->yRange ) ) ) )
            {
                if ( ( x < this->xStartIndex ) || ( x > ( this->xStartIndex + this->xRange ) ) ) { std::cout << "PROBLEM WITH LAT" << std::endl; }
                if ( ( y < this->yStartIndex ) || ( y > ( this->yStartIndex + this->yRange ) ) ) { std::cout << "PROBLEM WITH LON" << std::endl; }
                
                throw ProSHADE_exception ( "Requested bicubic interpolation outside of pre-computed\n                    : square.", "ES00064", __FILE__, __LINE__, __func__, "The supplied x or y value(s) is outside of the range of\n                    : the bi-cubic interpolator's pre-computed square. Please\n                    : make sure the start values were correctly supplied when\n                    : the constructor was called or create a new interpolator\n                    : for these values." );
            }
            
            //======================================== Convert x and y to unit square
            proshade_double unitSquareX               = ( x - this->xStartIndex ) / this->xRange;
            proshade_double unitSquareY               = ( y - this->yStartIndex ) / this->yRange;
            
            //======================================== Precompute powers
            proshade_double x2                        = std::pow ( unitSquareX, 2.0 );
            proshade_double x3                        = std::pow ( unitSquareX, 3.0 );
            proshade_double y2                        = std::pow ( unitSquareY, 2.0 );
            proshade_double y3                        = std::pow ( unitSquareY, 3.0 );

            //======================================== Done
            return                                    ( ( this->a00 + this->a01 * unitSquareY + this->a02 * y2 + this->a03 * y3 ) +
                                                        ( this->a10 + this->a11 * unitSquareY + this->a12 * y2 + this->a13 * y3 ) * unitSquareX +
                                                        ( this->a20 + this->a21 * unitSquareY + this->a22 * y2 + this->a23 * y3 ) * x2 +
                                                        ( this->a30 + this->a31 * unitSquareY + this->a32 * y2 + this->a33 * y3 ) * x3 );
        }
    };

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
    void getSOFTPositionFromEulerZXZ                  ( proshade_signed band, proshade_double eulerAlpha, proshade_double eulerBeta, proshade_double eulerGamma,
                                                        proshade_double* x, proshade_double* y, proshade_double* z );
    void getRotationMatrixFromEulerZXZAngles          ( proshade_double eulerAlpha, proshade_double eulerBeta, proshade_double eulerGamma, proshade_double* matrix );
    void getAxisAngleFromRotationMatrix               ( proshade_double* rotMat, proshade_double* x, proshade_double* y, proshade_double* z, proshade_double* ang );
    void getAxisAngleFromRotationMatrix               ( std::vector< proshade_double >* rotMat, proshade_double* x, proshade_double* y, proshade_double* z, proshade_double* ang );
    void getRotationMatrixFromAngleAxis               ( proshade_double* rotMat, proshade_double x, proshade_double y, proshade_double z, proshade_double ang );
    void getEulerZXZFromRotMatrix                     ( proshade_double* rotMat, proshade_double* eA, proshade_double* eB, proshade_double* eG );
    void getEulerZXZFromAngleAxis                     ( proshade_double axX, proshade_double axY, proshade_double axZ, proshade_double axAng, proshade_double* eA,
                                                        proshade_double* eB, proshade_double* eG, proshade_unsign angDim );
    void getEulerZXZFromAngleAxisFullSearch           ( proshade_double axX, proshade_double axY, proshade_double axZ, proshade_double axAng, proshade_double* eA,
                                                        proshade_double* eB, proshade_double* eG, proshade_signed angDim );
    void multiplyTwoSquareMatrices                    ( proshade_double* A, proshade_double* B, proshade_double* res, proshade_unsign dim = 3 );
    std::vector < proshade_signed > primeFactorsDecomp ( proshade_signed number );
    proshade_double normalDistributionValue           ( proshade_double mean, proshade_double standardDev, proshade_double value );
    proshade_double computeDotProduct                 ( proshade_double* x1, proshade_double* y1, proshade_double* z1, proshade_double* x2, proshade_double* y2,
                                                        proshade_double* z2 );
    proshade_double computeDotProduct                 ( proshade_double x1, proshade_double y1, proshade_double z1, proshade_double x2, proshade_double y2,
                                                        proshade_double z2 );
    proshade_double* computeCrossProduct              ( proshade_double* x1, proshade_double* y1, proshade_double* z1, proshade_double* x2, proshade_double* y2,
                                                        proshade_double* z2 );
    proshade_double* compute3x3MatrixMultiplication   ( proshade_double* mat1, proshade_double* mat2 );
    proshade_double* compute3x3MatrixVectorMultiplication ( proshade_double* mat, proshade_double x, proshade_double y, proshade_double z );
    proshade_double* compute3x3MatrixInverse          ( proshade_double* mat );
    void transpose3x3MatrixInPlace                    ( proshade_double* mat );
    proshade_double* findRotMatMatchingVectors        ( proshade_double x1, proshade_double y1, proshade_double z1, proshade_double x2, proshade_double y2,
                                                        proshade_double z2 );
    std::vector < proshade_double > findVectorFromTwoVAndTwoD ( proshade_double x1, proshade_double y1, proshade_double z1, proshade_double x2, proshade_double y2,
                                                                proshade_double z2, proshade_double dot1, proshade_double dot2 );
    std::vector < proshade_double > findVectorFromThreeVAndThreeD ( proshade_double x1, proshade_double y1, proshade_double z1, proshade_double x2, proshade_double y2,
                                                                    proshade_double z2, proshade_double x3, proshade_double y3, proshade_double z3, proshade_double dot1,
                                                                    proshade_double dot2, proshade_double dot3 );
    std::vector< proshade_double > multiplyGroupElementMatrices ( std::vector< proshade_double >* el1, std::vector< proshade_double >* el2 );
    bool rotationMatrixSimilarity                     ( std::vector< proshade_double >* mat1, std::vector< proshade_double >* mat2, proshade_double tolerance = 0.1 );
    bool vectorOrientationSimilarity                  ( proshade_double a1, proshade_double a2, proshade_double a3, proshade_double b1, proshade_double b2,
                                                        proshade_double b3, proshade_double tolerance = 0.1 );
    bool vectorOrientationSimilaritySameDirection     ( proshade_double a1, proshade_double a2, proshade_double a3, proshade_double b1, proshade_double b2,
                                                        proshade_double b3, proshade_double tolerance = 0.1 );
    void optimiseAxisBiCubicInterpolation             ( proshade_double* bestLattitude, proshade_double* bestLongitude, proshade_double* bestSum, std::vector<proshade_unsign>* sphereList,
                                                        std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*>* sphereMappedRotFun, proshade_double step = 0.05 );
    void prepareBiCubicInterpolatorsMinusMinus        ( proshade_double bestLattitude, proshade_double bestLongitude, std::vector<proshade_unsign>* sphereList,
                                                        std::vector<ProSHADE_internal_maths::BicubicInterpolator*>* interpols, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*>* sphereMappedRotFun );
    void prepareBiCubicInterpolatorsMinusPlus         ( proshade_double bestLattitude, proshade_double bestLongitude, std::vector<proshade_unsign>* sphereList,
                                                        std::vector<ProSHADE_internal_maths::BicubicInterpolator*>* interpols, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*>* sphereMappedRotFun );
    void prepareBiCubicInterpolatorsPlusMinus         ( proshade_double bestLattitude, proshade_double bestLongitude, std::vector<proshade_unsign>* sphereList,
                                                        std::vector<ProSHADE_internal_maths::BicubicInterpolator*>* interpols, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*>* sphereMappedRotFun );
    void prepareBiCubicInterpolatorsPlusPlus          ( proshade_double bestLattitude, proshade_double bestLongitude, std::vector<proshade_unsign>* sphereList,
                                                        std::vector<ProSHADE_internal_maths::BicubicInterpolator*>* interpols, std::vector<ProSHADE_internal_spheres::ProSHADE_rotFun_sphere*>* sphereMappedRotFun );
    bool isAxisUnique                                 ( std::vector< proshade_double* >* CSymList, proshade_double* axis, proshade_double tolerance = 0.1, bool improve = false );
    bool isAxisUnique                                 ( std::vector< proshade_double* >* CSymList, proshade_double X, proshade_double Y, proshade_double Z, proshade_double fold, proshade_double tolerance );
    std::vector< proshade_unsign > findAllPrimes      ( proshade_unsign upTo );
}

#endif
