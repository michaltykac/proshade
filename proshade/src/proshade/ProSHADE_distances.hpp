/*! \file ProSHADE_distances.hpp
    \brief This is the header file containing declarations of functions required for computation of shape distances.
 
    This header file contains the declarations of functions required to compute shape distances (specifically energy levels distances, trace sigma distances and the
    full rotation function distances) between two shapes. The user should not need to access these functions directly, as there are automated functions
    available in the higher levels of the ProSHADE organisation, which will call these in correct order and parse the results properly.
 
    Copyright by Michal Tykac and individual contributors. All rights reserved.
    
    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.
    
    This software is provided by the copyright holders and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are  disclaimed. In     no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or  services, loss of use, data     or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this  software, even if advised of the possibility     of such damage.
    
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.6.2
    \date      DEC 2021
*/

//============================================ ProSHADE
#include "ProSHADE_wignerMatrices.hpp"

//============================================ Overinclusion protection
#ifndef PROSHADE_DISTANCES
#define PROSHADE_DISTANCES

//============================================ ProSHADE_internal_distances Namespace

/*! \namespace ProSHADE_internal_distances
    \brief This namespace contains the functions used for computing distances between two structures as represented by the ProSHADE_data class.
 
 
    The ProSHADE_internal_distances namespace contains all the functionality required to compute distances between two or more ProSHADE_data
    objects and also all of the pre-computation functions.
 */
namespace ProSHADE_internal_distances
{
    proshade_double computeEnergyLevelsDescriptor     ( ProSHADE_internal_data::ProSHADE_data* obj1,
                                                        ProSHADE_internal_data::ProSHADE_data* obj2, ProSHADE_settings* settings );
    bool isBandWithinShell                            ( proshade_unsign bandInQuestion, proshade_unsign shellInQuestion,
                                                        ProSHADE_internal_spheres::ProSHADE_sphere** spheres );
    void computeRRPPearsonCoefficients                ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                        ProSHADE_settings* settings, proshade_unsign minCommonBands,
                                                        proshade_unsign minCommonShells, std::vector<proshade_double>* bandDists );
    void allocateTrSigmaWorkspace                     ( proshade_unsign minSpheres, proshade_unsign intOrder, proshade_double*& obj1Vals, proshade_double*& obj2Vals,
                                                        proshade_double*& GLabscissas, proshade_double*& glWeights, proshade_complex*& radiiVals );
    void computeSphericalHarmonicsMagnitude           ( ProSHADE_internal_data::ProSHADE_data* obj, proshade_unsign band, proshade_unsign order,
                                                        proshade_unsign radius, proshade_double* result );
    void computeEMatricesForLM                        ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                        proshade_unsign bandIter, proshade_unsign orderIter, proshade_complex* radiiVals, proshade_unsign integOrder,
                                                        proshade_double* abscissas, proshade_double* weights, proshade_double integRange, proshade_double sphereDist );
    proshade_double computeWeightsForEMatricesForLM   ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                        proshade_unsign bandIter, proshade_unsign orderIter, proshade_double* obj1Vals, proshade_double* obj2Vals,
                                                        proshade_unsign integOrder, proshade_double* abscissas, proshade_double* weights, proshade_single sphereDist );
    void releaseTrSigmaWorkspace                      ( proshade_double*& obj1Vals, proshade_double*& obj2Vals, proshade_double*& GLabscissas,
                                                        proshade_double*& glWeights, proshade_complex*& radiiVals );
    void computeEMatrices                             ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                        ProSHADE_settings* settings );
    void normaliseEMatrices                           ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                        ProSHADE_settings* settings );
    proshade_double computeTraceSigmaDescriptor       ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                        ProSHADE_settings* settings );
    void generateSO3CoeffsFromEMatrices               ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                        ProSHADE_settings* settings );
    void allocateInvSOFTWorkspaces                    ( proshade_complex*& work1, proshade_complex*& work2, proshade_double*& work3, proshade_unsign band );
    void prepareInvSOFTPlan                           ( fftw_plan* inverseSO3, int band, fftw_complex* work1, proshade_complex* invCoeffs );
    void releaseInvSOFTMemory                         ( proshade_complex*& work1, proshade_complex*& work2, proshade_double*& work3 );
    void computeInverseSOFTTransform                  ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                        ProSHADE_settings* settings );
    proshade_double computeRotationFunctionDescriptor ( ProSHADE_internal_data::ProSHADE_data* obj1, ProSHADE_internal_data::ProSHADE_data* obj2,
                                                        ProSHADE_settings* settings );
}

#endif
