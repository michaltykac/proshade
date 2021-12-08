/*! \file ProSHADE_overlay.hpp
    \brief This header file declares the functions required for the structure overlay computation.
 
    The function grouped in the ProSHADE_internal_overlay namespace and declared here deal with structure overlay computation, including
    FT coefficients equalising structure padding and inverse SOFT computation.
 
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
#include "ProSHADE_symmetry.hpp"

//==================================================== Overinclusion protection
#ifndef PROSHADE_OVERLAY
#define PROSHADE_OVERLAY

//==================================================== ProSHADE_internal_overlay Namespace
/*! \namespace ProSHADE_internal_overlay
    \brief This namespace contains all the functions required for map overlays.
 
    The map overlay task does have some specific functions that are betters groupped together for higher
    clarity of the code and also simple modifications later. This is where these live.
 */
namespace ProSHADE_internal_overlay
{
    void getOptimalRotation                           ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* staticStructure,
                                                        ProSHADE_internal_data::ProSHADE_data* movingStructure, proshade_double* eulA, proshade_double* eulB,
                                                        proshade_double* eulG );
    void getOptimalTranslation                        ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* staticStructure,
                                                        ProSHADE_internal_data::ProSHADE_data* movingStructure, proshade_double* trsX,
                                                        proshade_double* trsY, proshade_double* trsZ, proshade_double eulA, proshade_double eulB,
                                                        proshade_double eulG );
    void computeTranslationsFromPeak                  ( ProSHADE_internal_data::ProSHADE_data* staticStructure, ProSHADE_internal_data::ProSHADE_data* movingStructure, proshade_double *trsX,
                                                        proshade_double *trsY, proshade_double *trsZ );
    void computeBeforeAfterZeroCounts                 ( proshade_unsign* addXPre, proshade_unsign* addYPre, proshade_unsign* addZPre, proshade_unsign* addXPost,
                                                        proshade_unsign* addYPost, proshade_unsign* addZPost, proshade_unsign xDim, proshade_unsign yDim,
                                                        proshade_unsign zDim, proshade_unsign xDimIndices, proshade_unsign yDimIndices, proshade_unsign zDimIndices );
    void paddMapWithZeroes                            ( proshade_double* oldMap, proshade_double*& newMap, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim,
                                                        proshade_unsign xDimIndices, proshade_unsign yDimIndices, proshade_unsign zDimIndices, proshade_unsign addXPre,
                                                        proshade_unsign addYPre, proshade_unsign addZPre );
    void allocateTranslationFunctionMemory            ( fftw_complex*& tmpIn1, fftw_complex*& tmpOut1, fftw_complex*& tmpIn2, fftw_complex*& tmpOut2, fftw_complex*& resIn,
                                                        fftw_complex*& resOut, fftw_plan& forwardFourierObj1, fftw_plan& forwardFourierObj2, fftw_plan& inverseFourierCombo,
                                                        proshade_unsign xD, proshade_unsign yD, proshade_unsign zD );
    void freeTranslationFunctionMemory                ( fftw_complex*& tmpIn1, fftw_complex*& tmpOut1, fftw_complex*& tmpIn2, fftw_complex*& tmpOut2, fftw_complex*& resOut,
                                                        fftw_plan& forwardFourierObj1, fftw_plan& forwardFourierObj2, fftw_plan& inverseFourierCombo );
    void computeAngularThreshold                      ( std::vector<proshade_double>* lonCO, std::vector<proshade_double>* latCO, proshade_unsign angRes );
    void initialiseInverseSHComputation               ( proshade_unsign shBand, double*& sigR, double*& sigI, double*& rcoeffs, double*& icoeffs, double*& weights, double*& workspace,
                                                        fftw_plan& idctPlan, fftw_plan& ifftPlan );
}

#endif
