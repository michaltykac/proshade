/*! \file ProSHADE_overlay.hpp
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
#include "ProSHADE_symmetry.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE_OVERLAY__
#define __PROSHADE_OVERLAY__

//============================================ ProSHADE_internal_overlay Namespace
/*! \namespace ProSHADE_internal_overlay
 \brief This namespace contains all the functions required for map overlays.
 
 The map overlay task does have some specific functions that are betters groupped together for higher
 clarity of the code and also simple modifications later. This is where these live.
 */
namespace ProSHADE_internal_overlay
{
    void getOptimalRotation                   ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* staticStructure,
                                                ProSHADE_internal_data::ProSHADE_data* movingStructure, proshade_double* eulA, proshade_double* eulB,
                                                proshade_double* eulG );
    void getOptimalTranslation                ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* staticStructure,
                                                ProSHADE_internal_data::ProSHADE_data* movingStructure, proshade_double* trsX,
                                                proshade_double* trsY, proshade_double* trsZ, proshade_double eulA, proshade_double eulB,
                                                proshade_double eulG );
    void computeBeforeAfterZeroCounts         ( proshade_unsign* addXPre, proshade_unsign* addYPre, proshade_unsign* addZPre, proshade_unsign* addXPost,
                                                proshade_unsign* addYPost, proshade_unsign* addZPost, proshade_unsign xDim, proshade_unsign yDim,
                                                proshade_unsign zDim, proshade_unsign xDimIndices, proshade_unsign yDimIndices, proshade_unsign zDimIndices );
    void paddMapWithZeroes                    ( proshade_double* oldMap, proshade_double*& newMap, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim,
                                                proshade_unsign xDimIndices, proshade_unsign yDimIndices, proshade_unsign zDimIndices, proshade_unsign addXPre,
                                                proshade_unsign addYPre, proshade_unsign addZPre );
    void allocateTranslationFunctionMemory    ( fftw_complex*& tmpIn1, fftw_complex*& tmpOut1, fftw_complex*& tmpIn2, fftw_complex*& tmpOut2, fftw_complex*& resIn,
                                                fftw_complex*& resOut, fftw_plan& forwardFourierObj1, fftw_plan& forwardFourierObj2, fftw_plan& inverseFourierCombo,
                                                proshade_unsign xD, proshade_unsign yD, proshade_unsign zD );
    void combineFourierForTranslation         ( fftw_complex* tmpOut1, fftw_complex* tmpOut2, fftw_complex*& resOut, proshade_unsign xD, proshade_unsign yD,
                                                proshade_unsign zD );
    void findHighestValueInMap                ( fftw_complex* resIn, proshade_unsign xD, proshade_unsign yD, proshade_unsign zD, proshade_double* trsX,
                                                proshade_double* trsY, proshade_double* trsZ, proshade_double* mapPeak );
    void freeTranslationFunctionMemory        ( fftw_complex*& tmpIn1, fftw_complex*& tmpOut1, fftw_complex*& tmpIn2, fftw_complex*& tmpOut2, fftw_complex*& resOut,
                                                fftw_plan& forwardFourierObj1, fftw_plan& forwardFourierObj2, fftw_plan& inverseFourierCombo );
    void computeAngularThreshold              ( std::vector<proshade_double>* lonCO, std::vector<proshade_double>* latCO, proshade_unsign angRes );
    void initialiseInverseSHComputation       ( proshade_unsign shBand, double*& sigR, double*& sigI, double*& rcoeffs, double*& icoeffs, double*& weights, double*& workspace,
                                                fftw_plan& idctPlan, fftw_plan& ifftPlan );
}

#endif
