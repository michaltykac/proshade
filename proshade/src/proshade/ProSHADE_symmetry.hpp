/*! \file ProSHADE_symmetry.hpp
    \brief This header file declares all the functions required for symmetry detection and construction.
 
    This header file declares the ProSHADE_internal_symmetry namespace, which groups all the function required to detect symmetry axes and subsequently to detect symmetry type from
    the combination of all detected symmetry axes. Currently, this code can detect C, D, T, O and I symmetries with the C and D symmetries having their fold detected as well.
 
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
#include "ProSHADE_distances.hpp"

//==================================================== Overinclusion protection
#ifndef PROSHADE_SYMMETRY
#define PROSHADE_SYMMETRY

//==================================================== ProSHADE_internal_symmetry Namespace
/*! \namespace ProSHADE_internal_symmetry
    \brief This namespace contains the symmetry detection related code.
 
    The ProSHADE_internal_symmetry namespace contains the functions related to the symmetry detection task.
 */
namespace ProSHADE_internal_symmetry
{
    bool isSymmetrySame                               ( std::vector< proshade_double* >* ret, proshade_double* sym, proshade_double simThres, proshade_signed* matchedPos );
    bool isSymmetrySame                               ( std::vector< proshade_double* >* ret, proshade_double* sym, proshade_double simThres, proshade_signed* matchedPos, proshade_double fscVal );
    void saveDSymmetry                                ( std::vector< proshade_double* >* ret, std::vector< proshade_double* >* CSymList,
                                                        proshade_unsign axisOne, proshade_unsign axisTwo );
    bool detectTetrahedralSymmetry                    ( std::vector< proshade_double* >* CSymList, proshade_double axErr, proshade_double minPeakHeight );
    void findTetra4C3s                                ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret,
                                                        proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj,
                                                        proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight );
    bool testGroupAgainstSymmetry                     ( std::vector< proshade_double* >* CSymList, std::vector< proshade_unsign >* grp,
                                                        proshade_double* sym, proshade_double axErr, proshade_double angle, bool improve,
                                                        proshade_unsign pos = 0 );
    bool findMissingAxes                              ( std::vector< std::vector< proshade_unsign > >* possibilities,
                                                        std::vector< proshade_double* >* CSymList, proshade_unsign requiredNoAxes,
                                                        proshade_double axErr, proshade_double angle, proshade_unsign fold,
                                                        ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_double minPeakHeight );
    proshade_double missingAxisHeight                 ( proshade_double xVal, proshade_double yVal, proshade_double zVal,
                                                        ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_unsign fold, proshade_double axErr );
    std::vector < proshade_double* > findMissingAxisPoints ( proshade_double xVal, proshade_double yVal, proshade_double zVal,
                                                             ProSHADE_internal_data::ProSHADE_data* dataObj, proshade_double axErr );
    bool sortArrVecHlp                                ( const proshade_double* a, const proshade_double* b );
    void saveMissingAxisNewOnly                       ( std::vector< proshade_double* >* axVec, proshade_double axX, proshade_double axY,
                                                        proshade_double axZ, proshade_double height, proshade_unsign fold, proshade_double axErr );
    void searchMissingSymmetrySpace                   ( ProSHADE_internal_data::ProSHADE_data* dataObj, std::vector< proshade_double* >* CSymList,
                                                        std::vector< proshade_unsign >* grp, std::vector< proshade_double* >* hlpVec,
                                                        proshade_double axErr, proshade_double angle, proshade_unsign fold,
                                                        proshade_double minPeakHeight );
    void findTetra3C2s                                ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret,
                                                        proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj,
                                                        proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight );
    bool testGroupAgainstGroup                        ( std::vector< proshade_double* >* CSymList, std::vector< proshade_unsign >* grp1,
                                                        std::vector< proshade_double* >* RetList, std::vector< proshade_unsign >* grp2,
                                                        proshade_double angle, proshade_double axErr );
    bool detectOctahedralSymmetry                     ( std::vector< proshade_double* >* CSymList, proshade_double axErr, proshade_double minPeakHeight );
    void findOcta3C4s                                 ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret,
                                                        proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj,
                                                        proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight );
    void findOcta4C3s                                 ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret,
                                                        proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj,
                                                        proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight );
    void findOcta6C2s                                 ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret,
                                                        proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj,
                                                        proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight );
    bool findMissingAxesDual                          ( std::vector< proshade_unsign >* possibilities,
                                                        std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, std::vector<
                                                        proshade_unsign >* retGroup, proshade_unsign requiredNoAxes, proshade_double axErr,
                                                        proshade_unsign noMatchesG1, proshade_double angle1, proshade_unsign noMatchesG2,
                                                        proshade_double angle2, proshade_unsign fold,  ProSHADE_internal_data::ProSHADE_data* dataObj );
    proshade_signed addAxisUnlessSame                 ( proshade_unsign fold, proshade_double axX, proshade_double axY, proshade_double axZ,
                                                        proshade_double axHeight, std::vector< proshade_double* >* prosp, proshade_double axErr );
    proshade_signed addAxisUnlessSame                 ( proshade_unsign fold, proshade_double axX, proshade_double axY, proshade_double axZ,
                                                       proshade_double axHeight, proshade_double averageFSC, std::vector< proshade_double* >* prosp, proshade_double axErr );
    bool checkFittingAxisDualAndSave                  ( std::vector< proshade_unsign >* retGroup, std::vector< proshade_double* >* ret,
                                                        proshade_unsign fold, proshade_double axX, proshade_double axY, proshade_double axZ,
                                                        std::vector< proshade_double* >* prosp, proshade_double axErr, proshade_unsign noMatchesG1,
                                                        proshade_double angle1, proshade_unsign noMatchesG2, proshade_double angle2,
                                                        ProSHADE_internal_data::ProSHADE_data* dataObj );
    bool detectIcosahedralSymmetry                    ( std::vector< proshade_double* >* CSymList, proshade_double axErr, proshade_double minPeakHeight );
    void findIcos6C5s                                 ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret,
                                                        proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj,
                                                        proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight );
    void findIcos10C3s                                ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret,
                                                        proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj,
                                                        proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight );
    void findIcos15C2s                                ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret,
                                                        proshade_double axErr, ProSHADE_internal_data::ProSHADE_data* dataObj,
                                                        proshade_signed verbose, proshade_signed messageShift, proshade_double minPeakHeight );
    bool findMissingAxesTriple                        ( std::vector< proshade_unsign >* possibilities,
                                                        std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, std::vector<
                                                        proshade_unsign >* retGroup, proshade_unsign requiredNoAxes, proshade_double axErr,
                                                        proshade_unsign noMatchesG1, proshade_double angle1, proshade_unsign noMatchesG2,
                                                        proshade_double angle2, proshade_unsign noMatchesG3, proshade_double angle3,
                                                        proshade_unsign fold, ProSHADE_internal_data::ProSHADE_data* dataObj );
    void checkFittingAxisTripleAndSave                ( std::vector< proshade_unsign >* retGroup, std::vector< proshade_double* >* ret,
                                                        proshade_unsign fold, proshade_double axX, proshade_double axY, proshade_double axZ,
                                                        std::vector< proshade_double* >* prosp, proshade_double axErr,
                                                        proshade_unsign noMatchesG1, proshade_double angle1, proshade_unsign noMatchesG2,
                                                        proshade_double angle2, proshade_unsign noMatchesG3, proshade_double angle3,
                                                        ProSHADE_internal_data::ProSHADE_data* dataObj );
    proshade_double findPredictedSingleAxisHeight     ( proshade_double* axis, proshade_double fold, ProSHADE_internal_data::ProSHADE_data* dataObj, ProSHADE_settings* settings );
    void findPredictedAxesHeights                     ( std::vector< proshade_double* >* ret, ProSHADE_internal_data::ProSHADE_data* dataObj, ProSHADE_settings* settings );
    void optimiseDGroupAngleFromAxesHeights           ( std::vector < std::vector< proshade_double > >* ret, ProSHADE_internal_data::ProSHADE_data* dataObj, ProSHADE_settings* settings );
    void optimiseDGroupAngleFromAxesHeights           ( std::vector < std::vector< proshade_double > >* allCs, std::vector< proshade_unsign > selection, ProSHADE_internal_data::ProSHADE_data* dataObj,
                                                        ProSHADE_settings* settings );
    void predictIcosAxes                              ( std::vector< proshade_double* >* CSymList, std::vector< std::vector< proshade_double* > >* ret, proshade_double axErr, proshade_double minPeakHeight );
    void predictOctaAxes                              ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, proshade_double minPeakHeight );
    void predictTetraAxes                             ( std::vector< proshade_double* >* CSymList, std::vector< proshade_double* >* ret, proshade_double axErr, proshade_double minPeakHeight );
    std::vector< proshade_unsign > findReliableUnphasedSymmetries ( std::vector < std::vector< proshade_double > >* allCs, proshade_signed verbose, proshade_signed messageShift, proshade_double tolerance );
    void allocateCentreOfMapFourierTransforms         ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, fftw_complex *&origMap, fftw_complex *&origCoeffs, fftw_complex *&rotMapComplex,
                                                        fftw_complex *&rotCoeffs, fftw_complex *&trFunc, fftw_complex *&trFuncCoeffs, fftw_plan *planForwardFourier, fftw_plan *planForwardFourierRot,
                                                        fftw_plan *planReverseFourierComb );
    void releaseCentreOfMapFourierTransforms          ( fftw_complex *origMap, fftw_complex *origCoeffs, fftw_complex *rotMapComplex, fftw_complex *rotCoeffs, fftw_complex *trFunc, fftw_complex *trFuncCoeffs,
                                                        fftw_plan planForwardFourier, fftw_plan planForwardFourierRot, fftw_plan planReverseFourierComb );
    std::vector< proshade_double > findTranslationBetweenRotatedAndOriginalMap ( ProSHADE_internal_data::ProSHADE_data* symStr, std::vector < proshade_double > symElem, fftw_complex *origCoeffs,
                                                                                 fftw_complex* rotMapComplex, fftw_complex* rotCoeffs, fftw_plan planForwardFourierRot, fftw_complex* trFuncCoeffs,
                                                                                 fftw_complex* trFunc, fftw_plan planReverseFourierComb );
    std::vector< proshade_double > findPointFromTranslations ( ProSHADE_internal_data::ProSHADE_data* symStr, std::vector < std::vector < proshade_double > > symElems, fftw_complex *origCoeffs,
                                                               fftw_complex* rotMapComplex, fftw_complex* rotCoeffs, fftw_plan planForwardFourierRot, fftw_complex* trFuncCoeffs,
                                                               fftw_complex* trFunc, fftw_plan planReverseFourierComb );
}

#endif
