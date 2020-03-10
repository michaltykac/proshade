/*! \file ProSHADE_wignerMatrices.hpp
 \brief ...
 
 ...
 
 This file is part of the ProSHADE library for calculating
 shape descriptors and symmetry operators of protein structures.
 This is a prototype code, which is by no means complete or fully
 tested. Its use is at your own risk only. There is no quarantee
 that the results are correct.
 
 \author    Michal Tykac
 \author    Garib N. Murshudov
 \version   0.7.2
 \date      DEC 2019
 */

//============================================ ProSHADE library code
#include "ProSHADE_peakSearch.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE_WIGNERMATRICES__
#define __PROSHADE_WIGNERMATRICES__

//============================================ ProSHADE Namespace
/*! \namespace ProSHADE_internal_wigner
 \brief This namespace contains all the functionality required to perform object rotation in SO(3) space.
 
 The ProSHADE namespace wraps around all the functions required to compute Wigner d and D matrices for any perticular rotation (as given
 by the ZXZ Euler angles) and applying these to the SO(3) coefficients of any object.
 */
namespace ProSHADE_internal_wigner
{
    void allocateWignerWorkspace              ( proshade_double*& matIn, proshade_double*& matOut, proshade_double*& sqrts, proshade_double*& workspace,
                                                proshade_double*& alphaExponentReal, proshade_double*& alphaExponentImag, proshade_double*& gammaExponentReal,
                                                proshade_double*& gammaExponentImag, proshade_double*& trigs, proshade_unsign compBand );
    void releaseWignerWorkspace               ( proshade_double*& matIn, proshade_double*& matOut, proshade_double*& sqrts, proshade_double*& workspace,
                                                proshade_double*& alphaExponentReal, proshade_double*& alphaExponentImag, proshade_double*& gammaExponentReal,
                                                proshade_double*& gammaExponentImag, proshade_double*& trigs );
    void prepareTrigsSqrtsAndExponents        ( proshade_double* sqrts, proshade_double* alphaExponentReal, proshade_double* alphaExponentImag,
                                                proshade_double* gammaExponentReal, proshade_double* gammaExponentImag, proshade_double* trigs, proshade_unsign compBand,
                                                proshade_double angAlpha, proshade_double angBeta, proshade_double angGamma );
    void computeWignerMatrices                ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* obj, proshade_double* alphaExponentReal,
                                                proshade_double* alphaExponentImag, proshade_double* gammaExponentReal, proshade_double* gammaExponentImag,
                                                proshade_double* matIn, proshade_double* matOut, proshade_double* trigs, proshade_double* sqrts, proshade_double* workspace );
    void computeWignerMatricesForRotation     ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* obj, proshade_double eulerAlpha,
                                                proshade_double eulerBeta, proshade_double eulerGamma );
}

#endif
