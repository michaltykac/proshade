/*! \file ProSHADE_sphericalHarmonics.hpp
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

//============================================ ProSHADE
#include "ProSHADE_mapManip.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE_SPHERICAL_HARMONICS__
#define __PROSHADE_SPHERICAL_HARMONICS__

//============================================ ProSHADE_internal_sphericalHarmonics namespace
/*! \namespace ProSHADE_internal_sphericalHarmonics
 \brief This namespace contains the internal functions for computing spherical harmonics and their related computations.
 
 The ProSHADE_internal_sphericalHarmonics namespace contains helper functions for spherical harmonics computation for each sphere as
 created in the sphere object as well as the further processing computations. None of these functions should be used directly be
 the user.
 */
namespace ProSHADE_internal_sphericalHarmonics
{
    void allocateComputationMemory            ( proshade_unsign band, proshade_double*& inputReal, proshade_double*& inputImag, proshade_double*& outputReal,
                                                proshade_double*& outputImag, double*& shWeights, double*& tableSpaceHelper, fftw_complex*& workspace );
    void placeWithinWorkspacePointers         ( fftw_complex*& workspace, proshade_unsign oDim, proshade_double*& rres, proshade_double*& ires,
                                                proshade_double*& fltres, proshade_double*& scratchpad );
    void initialiseFFTWPlans                  ( proshade_unsign band, fftw_plan& fftPlan, fftw_plan& dctPlan, proshade_double*& inputReal,
                                                proshade_double*& inputImag, proshade_double*& rres, proshade_double*& ires, proshade_double*& scratchpad );
    void releaseSphericalMemory               ( proshade_double*& inputReal, proshade_double*& inputImag, proshade_double*& outputReal,
                                                proshade_double*& outputImag, double*& tableSpaceHelper, double**& tableSpace,
                                                double*& shWeights, fftw_complex*& workspace, fftw_plan& fftPlan, fftw_plan& dctPlan );
    void initialiseAllMemory                  ( proshade_unsign band, proshade_double*& inputReal, proshade_double*& inputImag, proshade_double*& outputReal,
                                                proshade_double*& outputImag, double*& shWeights, double**& tableSpace,
                                                double*& tableSpaceHelper, fftw_complex*& workspace, proshade_double*& rres, proshade_double*& ires,
                                                proshade_double*& fltres, proshade_double*& scratchpad, fftw_plan& fftPlan, fftw_plan& dctPlan );
    void initialSplitDiscreteTransform        ( proshade_unsign oneDim, proshade_double*& inputReal, proshade_double*& inputImag, proshade_double*& rres,
                                                proshade_double*& ires, proshade_double* mappedData, fftw_plan& fftPlan, proshade_double normCoeff );
    void computeSphericalTransformCoeffs      ( proshade_unsign band, proshade_double*& rdataptr, proshade_double*& idataptr, proshade_double*& outputReal,
                                                proshade_double*& outputImag, proshade_double*& rres, proshade_double*& ires, proshade_double*& fltres,
                                                proshade_double*& scratchpad, double**& tablePml, double*& shWeights, fftw_plan& dctPlan );
    void applyCondonShortleyPhase             ( proshade_unsign band, proshade_double* outputReal, proshade_double* outputImag, proshade_complex*& shArray );
    void computeSphericalHarmonics            ( proshade_unsign band, proshade_double* sphereMappedData, proshade_complex*& shArray );
}

#endif 
