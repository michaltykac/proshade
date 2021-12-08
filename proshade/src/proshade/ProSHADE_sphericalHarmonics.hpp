/*! \file ProSHADE_sphericalHarmonics.hpp
    \brief This header file declares the functions required to compute the spherical harmonics decomposition for a particular sphere.
 
    This header file declares the ProSHADE_internal_sphericalHarmonics namespace, which groups all the function required to compute the spherical harmonics
    decompostion on a surface of a sphere. There are also support functions for tasks like memory allocation and release and such.
 
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
#include "ProSHADE_mapManip.hpp"

//==================================================== Overinclusion protection
#ifndef PROSHADE_SPHERICAL_HARMONICS
#define PROSHADE_SPHERICAL_HARMONICS

//==================================================== ProSHADE_internal_sphericalHarmonics namespace
/*! \namespace ProSHADE_internal_sphericalHarmonics
    \brief This namespace contains the internal functions for computing spherical harmonics and their related computations.
 
    The ProSHADE_internal_sphericalHarmonics namespace contains helper functions for spherical harmonics computation for each sphere as
    created in the sphere object as well as the further processing computations. None of these functions should be used directly be
    the user.
 */
namespace ProSHADE_internal_sphericalHarmonics
{
    void allocateComputationMemory                    ( proshade_unsign band, proshade_double*& inputReal, proshade_double*& inputImag, proshade_double*& outputReal,
                                                        proshade_double*& outputImag, double*& shWeights, double*& tableSpaceHelper, fftw_complex*& workspace );
    void placeWithinWorkspacePointers                 ( fftw_complex*& workspace, proshade_unsign oDim, proshade_double*& rres, proshade_double*& ires,
                                                        proshade_double*& fltres, proshade_double*& scratchpad );
    void initialiseFFTWPlans                          ( proshade_unsign band, fftw_plan& fftPlan, fftw_plan& dctPlan, proshade_double*& inputReal,
                                                        proshade_double*& inputImag, proshade_double*& rres, proshade_double*& ires,
                                                        proshade_double*& scratchpad );
    void releaseSphericalMemory                       ( proshade_double*& inputReal, proshade_double*& inputImag, proshade_double*& outputReal,
                                                        proshade_double*& outputImag, double*& tableSpaceHelper, double**& tableSpace,
                                                        double*& shWeights, fftw_complex*& workspace, fftw_plan& fftPlan, fftw_plan& dctPlan );
    void initialiseAllMemory                          ( proshade_unsign band, proshade_double*& inputReal, proshade_double*& inputImag,
                                                        proshade_double*& outputReal,
                                                        proshade_double*& outputImag, double*& shWeights, double**& tableSpace,
                                                        double*& tableSpaceHelper, fftw_complex*& workspace, proshade_double*& rres, proshade_double*& ires,
                                                        proshade_double*& fltres, proshade_double*& scratchpad, fftw_plan& fftPlan, fftw_plan& dctPlan );
    void initialSplitDiscreteTransform                ( proshade_unsign oneDim, proshade_double*& inputReal, proshade_double*& inputImag,
                                                        proshade_double*& rres,
                                                        proshade_double*& ires, proshade_double* mappedData, fftw_plan& fftPlan, proshade_double normCoeff );
    void computeSphericalTransformCoeffs              ( proshade_unsign band, proshade_double*& rdataptr, proshade_double*& idataptr,
                                                        proshade_double*& outputReal,
                                                        proshade_double*& outputImag, proshade_double*& rres, proshade_double*& ires, proshade_double*& fltres,
                                                        proshade_double*& scratchpad, double**& tablePml, double*& shWeights, fftw_plan& dctPlan );
    void applyCondonShortleyPhase                     ( proshade_unsign band, proshade_double* outputReal, proshade_double* outputImag,
                                                        proshade_complex*& shArray );
    void computeSphericalHarmonics                    ( proshade_unsign band, proshade_double* sphereMappedData, proshade_complex*& shArray );
}

#endif 
