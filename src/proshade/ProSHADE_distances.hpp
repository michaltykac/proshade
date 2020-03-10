/*! \file ProSHADE_distances.hpp
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
#include "ProSHADE_wignerMatrices.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE_DISTANCES__
#define __PROSHADE_DISTANCES__

//============================================ ProSHADE_internal_distances Namespace

/*! \namespace ProSHADE_internal_distances
 \brief This namespace contains the functions used for computing distances between two structures as represented by the ProSHADE_data class.
 
 
 The ProSHADE_internal_distances namespace contains all the functionality required to compute distances between two or more ProSHADE_data
 objects and also all of the pre-computation functions.
 */
namespace ProSHADE_internal_distances
{
    proshade_double computeEnergyLevelsDescriptor                  ( ProSHADE_internal_data::ProSHADE_data* obj1,
                                                                     ProSHADE_internal_data::ProSHADE_data* obj2, ProSHADE_settings* settings );
    bool isBandWithinShell                                         ( proshade_unsign bandInQuestion, proshade_unsign shellInQuestion,
                                                                     ProSHADE_internal_spheres::ProSHADE_sphere** spheres );
    void computeRRPPearsonCoefficients                             ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                                     ProSHADE_settings* settings, proshade_unsign minCommonBands,
                                                                     proshade_unsign minCommonShells, std::vector<proshade_double>* bandDists );
    void allocateTrSigmaWorkspace                                  ( proshade_unsign minSpheres, proshade_unsign intOrder, proshade_double*& obj1Vals, proshade_double*& obj2Vals,
                                                                     proshade_double*& GLabscissas, proshade_double*& glWeights, proshade_complex*& radiiVals );
    void computeSphericalHarmonicsMagnitude                        ( ProSHADE_internal_data::ProSHADE_data* obj, proshade_unsign band, proshade_unsign order,
                                                                     proshade_unsign radius, proshade_double* result );
    void computeEMatricesForLM                                     ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                                     proshade_unsign bandIter, proshade_unsign orderIter, proshade_complex* radiiVals, proshade_unsign integOrder,
                                                                     proshade_double* abscissas, proshade_double* weights, proshade_double integRange, proshade_double sphereDist );
    proshade_double computeWeightsForEMatricesForLM                ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                                     proshade_unsign bandIter, proshade_unsign orderIter, proshade_double* obj1Vals, proshade_double* obj2Vals,
                                                                     proshade_unsign integOrder, proshade_double* abscissas, proshade_double* weights, proshade_double sphereDist );
    void releaseTrSigmaWorkspace                                   ( proshade_double*& obj1Vals, proshade_double*& obj2Vals, proshade_double*& GLabscissas,
                                                                     proshade_double*& glWeights, proshade_complex*& radiiVals );
    void computeEMatrices                                          ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                                     ProSHADE_settings* settings );
    void normaliseEMatrices                                        ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                                     ProSHADE_settings* settings );
    proshade_double computeTraceSigmaDescriptor                    ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                                     ProSHADE_settings* settings );
    void generateSO3CoeffsFromEMatrices                            ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                                     ProSHADE_settings* settings );
    void allocateInvSOFTWorkspaces                                 ( proshade_complex*& work1, proshade_complex*& work2, proshade_double*& work3, proshade_unsign band );
    void prepareInvSOFTPlan                                        ( fftw_plan* inverseSO3, proshade_unsign band, fftw_complex* work1, proshade_complex* invCoeffs );
    void releaseInvSOFTMemory                                      ( proshade_complex*& work1, proshade_complex*& work2, proshade_double*& work3 );
    void computeInverseSOFTTransform                               ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                                     ProSHADE_settings* settings );
    proshade_double computeRotationunctionDescriptor               ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                                     ProSHADE_settings* settings );
}

#endif
