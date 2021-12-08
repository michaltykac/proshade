/*! \file ProSHADE_wignerMatrices.hpp
    \brief This header declares the functions required to compute the Wigner D matrices.
 
    The functions declared in this header file in the ProSHADE_internal_wigner namespace are required to compute the Wigner D matrices, which in turn allow the computation of the inverse SOFT transform as well as
    rotation of internal map representations in the spherical harmonics coefficient space.
 
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

//==================================================== ProSHADE library code
#include "ProSHADE_data.hpp"

//==================================================== Overinclusion protection
#ifndef PROSHADE_WIGNERMATRICES
#define PROSHADE_WIGNERMATRICES

//==================================================== ProSHADE Namespace
/*! \namespace ProSHADE_internal_wigner
    \brief This namespace contains all the functionality required to perform object rotation in SO(3) space.
 
    The ProSHADE namespace wraps around all the functions required to compute Wigner d and D matrices for any perticular rotation (as given
    by the ZXZ Euler angles) and applying these to the SO(3) coefficients of any object.
 */
namespace ProSHADE_internal_wigner
{
    void allocateWignerWorkspace                      ( proshade_double*& matIn, proshade_double*& matOut, proshade_double*& sqrts, proshade_double*& workspace,
                                                        proshade_double*& alphaExponentReal, proshade_double*& alphaExponentImag, proshade_double*& gammaExponentReal,
                                                        proshade_double*& gammaExponentImag, proshade_double*& trigs, proshade_unsign compBand );
    void releaseWignerWorkspace                       ( proshade_double*& matIn, proshade_double*& matOut, proshade_double*& sqrts, proshade_double*& workspace,
                                                        proshade_double*& alphaExponentReal, proshade_double*& alphaExponentImag, proshade_double*& gammaExponentReal,
                                                        proshade_double*& gammaExponentImag, proshade_double*& trigs );
    void prepareTrigsSqrtsAndExponents                ( proshade_double* sqrts, proshade_double* alphaExponentReal, proshade_double* alphaExponentImag,
                                                        proshade_double* gammaExponentReal, proshade_double* gammaExponentImag, proshade_double* trigs, proshade_unsign compBand,
                                                        proshade_double angAlpha, proshade_double angBeta, proshade_double angGamma );
    void computeWignerMatrices                        ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* obj, proshade_double* alphaExponentReal,
                                                        proshade_double* alphaExponentImag, proshade_double* gammaExponentReal, proshade_double* gammaExponentImag,
                                                        proshade_double* matIn, proshade_double* matOut, proshade_double* trigs, proshade_double* sqrts, proshade_double* workspace );
    void computeWignerMatricesForRotation             ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* obj, proshade_double eulerAlpha,
                                                        proshade_double eulerBeta, proshade_double eulerGamma );
}

#endif
