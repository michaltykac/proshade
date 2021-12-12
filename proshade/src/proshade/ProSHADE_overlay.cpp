/*! \file ProSHADE_overlay.cpp
    \brief This source file contains the functions required for structure overlay computations.
 
    The function contained in this source file deal with computation of structure overlay, specifically finding the optimal rotation and
    translation including computing the inverse SOFT, structure padding for equalising the number of FT coefficients as similar.
 
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
#include "ProSHADE_overlay.hpp"

/*! \brief This function computes the overlay rotation function (i.e. the correlation function in SO(3) space).
 
    This function assumes it is called from the object to which the rotation function is to be assigned to (presumably the
    moving rather than static structure). It starts by computing the E matrices, normalising these using the Patterson-like
    normalisation, generating SO(3) coefficients from the E matrices and finally computing their inverse SOFT transform to
    get the rotation function.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] obj2 A pointer to the data class object of the other ( static ) structure.
 */
void ProSHADE_internal_data::ProSHADE_data::getOverlayRotationFunction ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* obj2 )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting rotation function computation.", settings->messageShift );
    
    //================================================ Compute un-weighted E matrices and their weights
    ProSHADE_internal_distances::computeEMatrices     ( obj2, this, settings );
    
    //================================================ Normalise E matrices by the magnitudes
    ProSHADE_internal_distances::normaliseEMatrices   ( obj2, this, settings );

    //================================================ Generate SO(3) coefficients
    ProSHADE_internal_distances::generateSO3CoeffsFromEMatrices ( obj2, this, settings );

    //================================================ Compute the inverse SO(3) Fourier Transform (SOFT) on the newly computed coefficients
    ProSHADE_internal_distances::computeInverseSOFTTransform ( obj2, this, settings );
    
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Rotation function obtained.", settings->messageShift );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds the optimal rotation between two structures as described by the settings object.
 
    This function takes the settings and two structure classes. It then reads in and processes both structures so that the
    globally optimal rotation overlay Euler angles are detected. However, this is only the case for rotation along the centre
    of the map of the second structure; therefore, either use Patterson data (usePhase = false), or be aware that better rotation
    may exist for different centre of rotation.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] obj1 A pointer to the data class object of the other ( static ) structure.
    \param[in] obj2 A pointer to the data class object of the first ( moving ) structure.
    \param[in] eulA The variable to which the best Euler alpha angle value will be saved to.
    \param[in] eulB The variable to which the best Euler beta angle value will be saved to.
    \param[in] eulG The variable to which the best Euler gamma angle value will be saved to.
 */
void ProSHADE_internal_overlay::getOptimalRotation ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* staticStructure, ProSHADE_internal_data::ProSHADE_data* movingStructure, proshade_double* eulA, proshade_double* eulB, proshade_double* eulG )
{
    //================================================ Read in the structures
    staticStructure->readInStructure                  ( settings->inputFiles.at(0), 0, settings );
    movingStructure->readInStructure                  ( settings->inputFiles.at(1), 1, settings );
    
    //================================================ Internal data processing  (COM, norm, mask, extra space)
    staticStructure->processInternalMap               ( settings );
    movingStructure->processInternalMap               ( settings );
    
    //================================================ Map to sphere
    staticStructure->mapToSpheres                     ( settings );
    movingStructure->mapToSpheres                     ( settings );
    
    //================================================ Get spherical harmonics
    staticStructure->computeSphericalHarmonics        ( settings );
    movingStructure->computeSphericalHarmonics        ( settings );
    
    //================================================ Get the rotation function of the pair
    movingStructure->getOverlayRotationFunction       ( settings, staticStructure );
    
    //================================================ Get inverse SO(3) map top peak Euler angle values
    ProSHADE_internal_peakSearch::getBestPeakEulerAngsNaive ( movingStructure->getInvSO3Coeffs (),
                                                              std::min ( staticStructure->getMaxBand(), movingStructure->getMaxBand() ) * 2,
                                                              eulA,  eulB,  eulG, settings );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds the optimal translation between two structures as described by the settings object given a rotation between the two objects.
 
    This function starts by loading and processing the structures according to the settings object (keeping phase is assumed, but callers
    responsibility). It then applies the required rotation to the second (moing) strucutre and then it follows with zero padding to make sure
    the structures have the same dimensions (again, it assumes map re-sampling was done,
    but setting to is callers responsibility). It then computes the translation function, finds the highest peak and returns the positions
    as well as height of this peak.
 
    \param[in] settings A pointer to settings class containing all the information required for map overlay computation.
    \param[in] staticStructure A pointer to the data class object of the other ( static ) structure.
    \param[in] movingStructure A pointer to the data class object of the first ( moving ) structure.
    \param[in] trsX The variable to which the best X-axis position value will be saved to.
    \param[in] trsY The variable to which the best Y-axis position value will be saved to.
    \param[in] trsZ The variable to which the best Z-axis position value will be saved to.
    \param[in] eulA The Euler alpha angle value, by which the moving structure is to be rotated by.
    \param[in] eulB The Euler beta angle value, by which the moving structure is to be rotated by.
    \param[in] eulG The Euler gamma angle value, by which the moving structure is to be rotated by.
 */
void ProSHADE_internal_overlay::getOptimalTranslation ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* staticStructure, ProSHADE_internal_data::ProSHADE_data* movingStructure, proshade_double* trsX, proshade_double* trsY, proshade_double* trsZ, proshade_double eulA, proshade_double eulB, proshade_double eulG )
{
    //================================================ Read in the structures
    staticStructure->readInStructure                  ( settings->inputFiles.at(0), 0, settings );
    movingStructure->readInStructure                  ( settings->inputFiles.at(1), 1, settings );
    
    //================================================ Determine spherical harmonics variables from the first structure (otherwise, they would be determined from the second structure)
    settings->determineAllSHValues                    ( staticStructure->xDimIndices, staticStructure->yDimIndices,
                                                       staticStructure->xDimSize,     staticStructure->yDimSize,    staticStructure->zDimSize );
    
    //================================================ Internal data processing  (COM, norm, mask, extra space)
    staticStructure->processInternalMap               ( settings );
    movingStructure->processInternalMap               ( settings );
    
    //================================================ Rotate map
    std::vector< proshade_double > rotCen             = movingStructure->rotateMapRealSpaceInPlace ( eulA, eulB, eulG );
    
    //================================================ Zero padding for smaller structure
    staticStructure->zeroPaddToDims                   ( std::max ( staticStructure->getXDim(), movingStructure->getXDim() ),
                                                        std::max ( staticStructure->getYDim(), movingStructure->getYDim() ),
                                                        std::max ( staticStructure->getZDim(), movingStructure->getZDim() ) );
    movingStructure->zeroPaddToDims                   ( std::max ( staticStructure->getXDim(), movingStructure->getXDim() ),
                                                        std::max ( staticStructure->getYDim(), movingStructure->getYDim() ),
                                                        std::max ( staticStructure->getZDim(), movingStructure->getZDim() ) );
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting translation function computation." );
    
    //================================================ Compute the translation function map
    movingStructure->computeTranslationMap            ( staticStructure );
    
    //================================================ Find highest peak in translation function
    proshade_double mapPeak                           = 0.0;
    proshade_unsign xDimS                             = staticStructure->getXDim();
    proshade_unsign yDimS                             = staticStructure->getYDim();
    proshade_unsign zDimS                             = staticStructure->getZDim();
    ProSHADE_internal_maths::findHighestValueInMap  ( movingStructure->getTranslationFnPointer(), xDimS, yDimS, zDimS, trsX, trsY, trsZ, &mapPeak );
    
    //================================================ Dont translate over half
    if ( *trsX > ( static_cast< proshade_double > ( xDimS ) / 2.0 ) ) { *trsX = *trsX - static_cast< proshade_double > ( xDimS ); }
    if ( *trsY > ( static_cast< proshade_double > ( yDimS ) / 2.0 ) ) { *trsY = *trsY - static_cast< proshade_double > ( yDimS ); }
    if ( *trsZ > ( static_cast< proshade_double > ( zDimS ) / 2.0 ) ) { *trsZ = *trsZ - static_cast< proshade_double > ( zDimS ); }
    
    //================================================ Convert the map positions onto translation in Angstroms and save the results
    computeTranslationsFromPeak                       ( staticStructure, movingStructure, trsX, trsY, trsZ );
    
    movingStructure->originalPdbTransX                = *trsX;
    movingStructure->originalPdbTransY                = *trsY;
    movingStructure->originalPdbTransZ                = *trsZ;
    
    proshade_single xMov                              = static_cast< proshade_single > ( -(*trsX) );
    proshade_single yMov                              = static_cast< proshade_single > ( -(*trsY) );
    proshade_single zMov                              = static_cast< proshade_single > ( -(*trsZ) );
    
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "Optimal map translation distances are " << *trsX << " ; " << *trsY << " ; " << *trsZ << " Angstroms with peak height " << mapPeak / ( static_cast< proshade_double > ( xDimS ) * static_cast< proshade_double > ( yDimS ) * static_cast< proshade_double > ( zDimS ) );

    //================================================ Save original from variables for PDB output
    movingStructure->mapMovFromsChangeX               = static_cast< proshade_double > ( movingStructure->xFrom );
    movingStructure->mapMovFromsChangeY               = static_cast< proshade_double > ( movingStructure->yFrom );
    movingStructure->mapMovFromsChangeZ               = static_cast< proshade_double > ( movingStructure->zFrom );
    
    //================================================ Move the map
    ProSHADE_internal_mapManip::moveMapByIndices      ( &xMov, &yMov, &zMov, movingStructure->getXDimSize(), movingStructure->getYDimSize(), movingStructure->getZDimSize(),
                                                        movingStructure->getXFromPtr(), movingStructure->getXToPtr(),
                                                        movingStructure->getYFromPtr(), movingStructure->getYToPtr(),
                                                        movingStructure->getZFromPtr(), movingStructure->getZToPtr(),
                                                        movingStructure->getXAxisOrigin(), movingStructure->getYAxisOrigin(), movingStructure->getZAxisOrigin() );

    ProSHADE_internal_mapManip::moveMapByFourier      ( movingStructure->getInternalMap(), xMov, yMov, zMov,
                                                        movingStructure->getXDimSize(), movingStructure->getYDimSize(), movingStructure->getZDimSize(),
                                                        static_cast< proshade_signed > ( movingStructure->getXDim() ), static_cast< proshade_signed > ( movingStructure->getYDim() ),
                                                        static_cast< proshade_signed > ( movingStructure->getZDim() ) );
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, hlpSS.str(), settings->messageShift );
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Translation function computation complete.", settings->messageShift );
    
    //================================================ Keep only the change in from and to variables
    movingStructure->mapMovFromsChangeX               = static_cast< proshade_double > ( movingStructure->xFrom ) - movingStructure->mapMovFromsChangeX;
    movingStructure->mapMovFromsChangeY               = static_cast< proshade_double > ( movingStructure->yFrom ) - movingStructure->mapMovFromsChangeY;
    movingStructure->mapMovFromsChangeZ               = static_cast< proshade_double > ( movingStructure->zFrom ) - movingStructure->mapMovFromsChangeZ;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the translation in Angstroms that corresponds to the translation function peak.

    \param[in] staticStructure A pointer to the data class object of the static structure.
    \param[in] movingStructure A pointer to the data class object of the moving structure.
    \param[in] trsX The variable to which the best X-axis position value will be saved to.
    \param[in] trsY The variable to which the best Y-axis position value will be saved to.
    \param[in] trsZ The variable to which the best Z-axis position value will be saved to.
*/
void ProSHADE_internal_overlay::computeTranslationsFromPeak ( ProSHADE_internal_data::ProSHADE_data* staticStructure, ProSHADE_internal_data::ProSHADE_data* movingStructure, proshade_double *trsX, proshade_double *trsY, proshade_double *trsZ )
{
    //================================================ Compute map movement
    proshade_single xSamplRate                        = ( static_cast< proshade_single > ( staticStructure->getXDimSize() ) / static_cast< proshade_single > ( staticStructure->getXDim() ) );
    proshade_single ySamplRate                        = ( static_cast< proshade_single > ( staticStructure->getYDimSize() ) / static_cast< proshade_single > ( staticStructure->getYDim() ) );
    proshade_single zSamplRate                        = ( static_cast< proshade_single > ( staticStructure->getZDimSize() ) / static_cast< proshade_single > ( staticStructure->getZDim() ) );
    
    proshade_single xCor                              = static_cast< proshade_single > ( static_cast< proshade_single > ( staticStructure->xFrom - movingStructure->xFrom ) * xSamplRate );
    proshade_single yCor                              = static_cast< proshade_single > ( static_cast< proshade_single > ( staticStructure->yFrom - movingStructure->yFrom ) * ySamplRate );
    proshade_single zCor                              = static_cast< proshade_single > ( static_cast< proshade_single > ( staticStructure->zFrom - movingStructure->zFrom ) * zSamplRate );
    
    proshade_single xMov                              = static_cast< proshade_single > ( 0.0 - static_cast< proshade_double > ( xCor ) - ( *trsX * static_cast< proshade_double > ( xSamplRate ) ) );
    proshade_single yMov                              = static_cast< proshade_single > ( 0.0 - static_cast< proshade_double > ( yCor ) - ( *trsY * static_cast< proshade_double > ( ySamplRate ) ) );
    proshade_single zMov                              = static_cast< proshade_single > ( 0.0 - static_cast< proshade_double > ( zCor ) - ( *trsZ * static_cast< proshade_double > ( zSamplRate ) ) );
    
    //================================================ Save translation vector back
   *trsX                                              = static_cast< proshade_double > ( -xMov );
   *trsY                                              = static_cast< proshade_double > ( -yMov );
   *trsZ                                              = static_cast< proshade_double > ( -zMov );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function gets the optimal translation vector and returns it as a standard library vector. It also applies the translation to the internal map.

    \param[in] staticStructure A pointer to the data class object of the other ( static ) structure.
    \param[out] X A vector of doubles with the optimal translation vector in Angstroms.
*/
std::vector< proshade_double > ProSHADE_internal_data::ProSHADE_data::getBestTranslationMapPeaksAngstrom ( ProSHADE_internal_data::ProSHADE_data* staticStructure )
{
    //================================================ Initialise local variables
    std::vector< proshade_double > ret;
    proshade_double mapPeak                           = 0.0;
    proshade_unsign xDimS                             = staticStructure->getXDim();
    proshade_unsign yDimS                             = staticStructure->getYDim();
    proshade_unsign zDimS                             = staticStructure->getZDim();
    proshade_double trsX = 0.0, trsY = 0.0, trsZ = 0.0;
    ProSHADE_internal_maths::findHighestValueInMap    ( this->getTranslationFnPointer(),
                                                        xDimS,
                                                        yDimS,
                                                        zDimS,
                                                       &trsX,
                                                       &trsY,
                                                       &trsZ,
                                                       &mapPeak );
    
    //================================================ Dont translate over half
    if ( trsX > ( static_cast< proshade_double > ( xDimS ) / 2.0 ) ) { trsX = trsX - static_cast< proshade_double > ( xDimS ); }
    if ( trsY > ( static_cast< proshade_double > ( yDimS ) / 2.0 ) ) { trsY = trsY - static_cast< proshade_double > ( yDimS ); }
    if ( trsZ > ( static_cast< proshade_double > ( zDimS ) / 2.0 ) ) { trsZ = trsZ - static_cast< proshade_double > ( zDimS ); }
    
    //================================================ Convert the map positions onto translation in Angstroms and save the results
    ProSHADE_internal_overlay::computeTranslationsFromPeak ( staticStructure, this, &trsX, &trsY, &trsZ );
    
    this->originalPdbTransX                           = trsX;
    this->originalPdbTransY                           = trsY;
    this->originalPdbTransZ                           = trsZ;
    
    proshade_single xMov                              = static_cast< proshade_single > ( -trsX );
    proshade_single yMov                              = static_cast< proshade_single > ( -trsY );
    proshade_single zMov                              = static_cast< proshade_single > ( -trsZ );
    
    //================================================ Save results as vector
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, static_cast<proshade_double> ( -xMov ) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, static_cast<proshade_double> ( -yMov ) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, static_cast<proshade_double> ( -zMov ) );

    //================================================ Save original from variables for PDB output
    this->mapMovFromsChangeX                          = static_cast<proshade_double> ( this->xFrom );
    this->mapMovFromsChangeY                          = static_cast<proshade_double> ( this->yFrom );
    this->mapMovFromsChangeZ                          = static_cast<proshade_double> ( this->zFrom );
    
    //================================================ Move the map
    ProSHADE_internal_mapManip::moveMapByIndices      ( &xMov, &yMov, &zMov, this->getXDimSize(), this->getYDimSize(), this->getZDimSize(),
                                                         this->getXFromPtr(), this->getXToPtr(),
                                                         this->getYFromPtr(), this->getYToPtr(),
                                                         this->getZFromPtr(), this->getZToPtr(),
                                                         this->getXAxisOrigin(), this->getYAxisOrigin(), this->getZAxisOrigin() );
    
    ProSHADE_internal_mapManip::moveMapByFourier      ( this->getInternalMap(), xMov, yMov, zMov,
                                                        this->getXDimSize(), this->getYDimSize(), this->getZDimSize(),
                                                        static_cast< proshade_signed > ( this->getXDim() ), static_cast< proshade_signed > ( this->getYDim() ), static_cast< proshade_signed > ( this->getZDim() ) );

    //================================================ Keep only the change in from and to variables
    this->mapMovFromsChangeX                          = static_cast<proshade_double> ( this->xFrom ) - this->mapMovFromsChangeX;
    this->mapMovFromsChangeY                          = static_cast<proshade_double> ( this->yFrom ) - this->mapMovFromsChangeY;
    this->mapMovFromsChangeZ                          = static_cast<proshade_double> ( this->zFrom ) - this->mapMovFromsChangeZ;
    
   //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function does the computation of the translation map and saves results internally.

    This function takes the static structure, the optimal translation to which should be found and then it
    proceeds to compute the Fourier transform of both this and the static structures. It then combines the
    coefficients for translation function and computes the inverse Fourier transform, thus obtaining the
    translation function. This function is then saved, while all other internal data are deleted.

    \param[in] staticStructure A pointer to the data class object of the other ( static ) structure.
*/
void ProSHADE_internal_data::ProSHADE_data::computeTranslationMap ( ProSHADE_internal_data::ProSHADE_data* staticStructure )
{
    //================================================ Do this using Fourier!
    fftw_complex *tmpIn1 = nullptr, *tmpOut1 = nullptr, *tmpIn2 = nullptr, *tmpOut2 = nullptr, *resOut = nullptr;
    fftw_plan forwardFourierObj1, forwardFourierObj2, inverseFourierCombo;
    proshade_unsign dimMult                           = staticStructure->getXDim() * staticStructure->getYDim() * staticStructure->getZDim();
    ProSHADE_internal_overlay::allocateTranslationFunctionMemory ( tmpIn1, tmpOut1, tmpIn2, tmpOut2, this->translationMap, resOut, forwardFourierObj1, forwardFourierObj2, inverseFourierCombo, staticStructure->getXDim(), staticStructure->getYDim(), staticStructure->getZDim() );
    
    //================================================ Fill in input data
    for ( proshade_unsign iter = 0; iter < dimMult; iter++ ) { tmpIn1[iter][0] = staticStructure->getMapValue ( iter ); tmpIn1[iter][1] = 0.0; }
    for ( proshade_unsign iter = 0; iter < dimMult; iter++ ) { tmpIn2[iter][0] = this->getMapValue            ( iter ); tmpIn2[iter][1] = 0.0; }
    
    //================================================ Calculate Fourier
    fftw_execute                                      ( forwardFourierObj1 );
    fftw_execute                                      ( forwardFourierObj2 );
    
    //================================================ Combine Fourier coeffs and invert
    ProSHADE_internal_maths::combineFourierForTranslation ( tmpOut1, tmpOut2, resOut, staticStructure->getXDim(), staticStructure->getYDim(), staticStructure->getZDim() );
    fftw_execute                                      ( inverseFourierCombo );
    
    //================================================ Free memory
    ProSHADE_internal_overlay::freeTranslationFunctionMemory ( tmpIn1, tmpOut1, tmpIn2, tmpOut2, resOut, forwardFourierObj1, forwardFourierObj2, inverseFourierCombo );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allocates the memory for the Fourier transforms required for translation function computation.
 
    \param[in] tmpIn1 Array to hold the static structure Fourier inputs.
    \param[in] tmpOut1 Array to hold the static structure Fourier outputs.
    \param[in] tmpIn2 Array to hold the moving structure Fourier inputs.
    \param[in] tmpOut2 Array to hold the moving structure Fourier outputs.
    \param[in] resIn Array to hold the final translation function values.
    \param[in] resOut Array to hold the combined Fourier coefficients of both structures.
    \param[in] forwardFourierObj1 FFTW plan for the forward Fourier of the static structure.
    \param[in] forwardFourierObj2 FFTW plan for the forward Fourier of the moving structure.
    \param[in] inverseFourierCombo FFTW plan for the backward Fourier of the combined Fourier factors of both structures.
    \param[in] xD The dimension of the X axis of the structures (assumes both structures have the same sizes and sampling).
    \param[in] yD The dimension of the Y axis of the structures (assumes both structures have the same sizes and sampling).
    \param[in] zD The dimension of the Z axis of the structures (assumes both structures have the same sizes and sampling).
 */
void ProSHADE_internal_overlay::allocateTranslationFunctionMemory ( fftw_complex*& tmpIn1, fftw_complex*& tmpOut1, fftw_complex*& tmpIn2, fftw_complex*& tmpOut2, fftw_complex*& resIn, fftw_complex*& resOut, fftw_plan& forwardFourierObj1, fftw_plan& forwardFourierObj2, fftw_plan& inverseFourierCombo, proshade_unsign xD, proshade_unsign yD, proshade_unsign zD )
{
    //================================================ Allocate memory
    tmpIn1                                            = new fftw_complex[xD * yD * zD];
    tmpOut1                                           = new fftw_complex[xD * yD * zD];
    tmpIn2                                            = new fftw_complex[xD * yD * zD];
    tmpOut2                                           = new fftw_complex[xD * yD * zD];
    resIn                                             = new fftw_complex[xD * yD * zD];
    resOut                                            = new fftw_complex[xD * yD * zD];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( tmpIn1 , __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( tmpOut1, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( tmpIn2,  __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( tmpOut2, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( resIn,   __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( resOut,  __FILE__, __LINE__, __func__ );
    
    //================================================ Get Fourier transforms of the maps
    forwardFourierObj1                                = fftw_plan_dft_3d ( static_cast< int > ( xD ), static_cast< int > ( yD ), static_cast< int > ( zD ), tmpIn1, tmpOut1, FFTW_FORWARD , FFTW_ESTIMATE );
    forwardFourierObj2                                = fftw_plan_dft_3d ( static_cast< int > ( xD ), static_cast< int > ( yD ), static_cast< int > ( zD ), tmpIn2, tmpOut2, FFTW_FORWARD , FFTW_ESTIMATE );
    inverseFourierCombo                               = fftw_plan_dft_3d ( static_cast< int > ( xD ), static_cast< int > ( yD ), static_cast< int > ( zD ), resOut, resIn  , FFTW_BACKWARD, FFTW_ESTIMATE );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function releases the memory for the Fourier transforms required for translation function computation.
 
    \param[in] tmpIn1 Array to hold the static structure Fourier inputs.
    \param[in] tmpOut1 Array to hold the static structure Fourier outputs.
    \param[in] tmpIn2 Array to hold the moving structure Fourier inputs.
    \param[in] tmpOut2 Array to hold the moving structure Fourier outputs.
    \param[in] resOut Array to hold the combined Fourier coefficients of both structures.
    \param[in] forwardFourierObj1 FFTW plan for the forward Fourier of the static structure.
    \param[in] forwardFourierObj2 FFTW plan for the forward Fourier of the moving structure.
    \param[in] inverseFourierCombo FFTW plan for the backward Fourier of the combined Fourier factors of both structures.
 */
void ProSHADE_internal_overlay::freeTranslationFunctionMemory ( fftw_complex*& tmpIn1, fftw_complex*& tmpOut1, fftw_complex*& tmpIn2, fftw_complex*& tmpOut2, fftw_complex*& resOut, fftw_plan& forwardFourierObj1, fftw_plan& forwardFourierObj2, fftw_plan& inverseFourierCombo )
{
    //================================================ Release memory
    fftw_destroy_plan                                 ( forwardFourierObj1 );
    fftw_destroy_plan                                 ( forwardFourierObj2 );
    fftw_destroy_plan                                 ( inverseFourierCombo );
    delete[] tmpIn1;
    delete[] tmpIn2;
    delete[] tmpOut1;
    delete[] tmpOut2;
    delete[] resOut;
    
    //======================================== Done
    return ;
    
}

/*! \brief This function changes the size of a structure to fit the supplied new limits.
 
    This function increases the map size by symetrically adding zeroes in each required dimension. The first zero is
    always added AFTER the structure, so for even size increases, there will be misplacement of centre of mass. The
    map position in the "real" world should not change.
 
    \param[in] xDim The X dimension size to which this structure should be padded into.
    \param[in] yDim The Y dimension size to which this structure should be padded into.
    \param[in] zDim The Z dimension size to which this structure should be padded into.
 */
void ProSHADE_internal_data::ProSHADE_data::zeroPaddToDims ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim )
{
    //================================================ Sanity check
    if ( ( this->xDimIndices > xDim  ) || ( this->yDimIndices > yDim  ) || ( this->zDimIndices > zDim  ) )
    {
        throw ProSHADE_exception ( "Cannot zero-pad in negative direction.", "EO00034", __FILE__, __LINE__, __func__, "The requested padded size of a structure is smaller than\n                    : the current size. If the user sees this error, there is\n                    : likely a considerable bug. Please report this error." );
    }
    
    //================================================ If done, do nothing
    if ( ( this->xDimIndices == xDim  ) && ( this->yDimIndices == yDim  ) && ( this->zDimIndices == zDim  ) ) { return ; }
    
    //================================================ Find out how many zeroes need to be added before and after
    proshade_unsign addXPre, addYPre, addZPre, addXPost, addYPost, addZPost;
    ProSHADE_internal_overlay::computeBeforeAfterZeroCounts ( &addXPre, &addYPre, &addZPre, &addXPost, &addYPost, &addZPost, xDim, yDim, zDim, this->xDimIndices, this->yDimIndices, this->zDimIndices );
    
    //================================================ Create a new map
    proshade_double* newMap                           = new proshade_double [xDim * yDim * zDim];
    
    //================================================ Do the hard work
    ProSHADE_internal_overlay::paddMapWithZeroes      ( this->internalMap, newMap, xDim, yDim, zDim, this->xDimIndices, this->yDimIndices, this->zDimIndices, addXPre, addYPre, addZPre );
    
    //================================================ Create a new internal map and copy
    delete[] this->internalMap;
    this->internalMap                                 = new proshade_double [xDim * yDim * zDim];
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( xDim * yDim * zDim ); iter++ ) { this->internalMap[iter] = newMap[iter]; }
    
    //================================================ Release memory
    delete[] newMap;
    
    //================================================ Update map related variables
    this->xDimSize                                    = static_cast< proshade_single > ( xDim ) * ( this->xDimSize / static_cast< proshade_single > ( this->xDimIndices ) );
    this->yDimSize                                    = static_cast< proshade_single > ( yDim ) * ( this->yDimSize / static_cast< proshade_single > ( this->yDimIndices ) );
    this->zDimSize                                    = static_cast< proshade_single > ( zDim ) * ( this->zDimSize / static_cast< proshade_single > ( this->zDimIndices ) );
    this->xDimIndices  = xDim    ; this->yDimIndices  = yDim    ; this->zDimIndices  = zDim;
    this->xGridIndices = xDim    ; this->yGridIndices = yDim    ; this->zGridIndices = zDim;
    this->xFrom       -= addXPre ; this->yFrom       -= addYPre ; this->zFrom       -= addZPre;
    this->xTo         += addXPost; this->yTo         += addYPost; this->zTo         += addZPost;
    this->xAxisOrigin -= addXPre ; this->yAxisOrigin -= addYPre ; this->zAxisOrigin -= addZPre ;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds the number of zeroes to be added after and before the structure along each dimension.
 
    \param[in] addXPre Variable pointer where the number of zeroes to be added before the map along X axis is saved.
    \param[in] addYPre Variable pointer where the number of zeroes to be added before the map along Y axis is saved.
    \param[in] addZPre Variable pointer where the number of zeroes to be added before the map along Z axis is saved.
    \param[in] addXPost Variable pointer where the number of zeroes to be added after the map along X axis is saved.
    \param[in] addYPost Variable pointer where the number of zeroes to be added after the map along Y axis is saved.
    \param[in] addZPost Variable pointer where the number of zeroes to be added after the map along Z axis is saved.
    \param[in] xDim The X dimension size in indices of the new map.
    \param[in] yDim The Y dimension size in indices of the new map.
    \param[in] zDim The Z dimension size in indices of the new map.
    \param[in] xDimIndices The X dimension size in indices of the old map.
    \param[in] yDimIndices The Y dimension size in indices of the old map.
    \param[in] zDimIndices The Z dimension size in indices of the old map.
 */
void ProSHADE_internal_overlay::computeBeforeAfterZeroCounts ( proshade_unsign* addXPre, proshade_unsign* addYPre, proshade_unsign* addZPre, proshade_unsign* addXPost, proshade_unsign* addYPost, proshade_unsign* addZPost, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_unsign xDimIndices, proshade_unsign yDimIndices, proshade_unsign zDimIndices )
{
    //================================================ Compute
   *addXPre                                           = ( xDim - xDimIndices ) / 2;
   *addYPre                                           = ( yDim - yDimIndices ) / 2;
   *addZPre                                           = ( zDim - zDimIndices ) / 2;
   *addXPost                                          = *addXPre; if ( ( xDim - xDimIndices ) % 2 == 1 ) { *addXPost += 1; }
   *addYPost                                          = *addYPre; if ( ( yDim - yDimIndices ) % 2 == 1 ) { *addYPost += 1; }
   *addZPost                                          = *addZPre; if ( ( zDim - zDimIndices ) % 2 == 1 ) { *addZPost += 1; }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function adds zeroes before and after the central map and copies the central map values into a new map.
 
    \param[in] oldMap The original, unpadded map.
    \param[in] newMap The array to which the new map will be saved into.
    \param[in] xDim The X dimension size of the new map.
    \param[in] yDim The Y dimension size of the new map.
    \param[in] zDim The Z dimension size of the new map.
    \param[in] xDimIndices The X dimension size in indices of the old map.
    \param[in] yDimIndices The Y dimension size in indices of the old map.
    \param[in] zDimIndices The Z dimension size in indices of the old map.
    \param[in] addXPre How many zeroes are to be added before the central map along the X axis?
    \param[in] addYPre How many zeroes are to be added before the central map along the Y axis?
    \param[in] addZPre How many zeroes are to be added before the central map along the Z axis?
 */
void ProSHADE_internal_overlay::paddMapWithZeroes ( proshade_double* oldMap, proshade_double*& newMap, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_unsign xDimIndices, proshade_unsign yDimIndices, proshade_unsign zDimIndices, proshade_unsign addXPre, proshade_unsign addYPre, proshade_unsign addZPre )
{
    //================================================ Initialise local variables
    proshade_unsign newMapIndex                       = 0;
    proshade_unsign oldMapIndex                       = 0;
    
    //================================================ Copy and padd as appropriate
    for ( proshade_unsign xIt = 0; xIt < xDim; xIt++ )
    {
        for ( proshade_unsign yIt = 0; yIt < yDim; yIt++ )
        {
            for ( proshade_unsign zIt = 0; zIt < zDim; zIt++ )
            {
                //==================================== Find map location
                newMapIndex                           = zIt + zDim * ( yIt + yDim * xIt );
                
                //==================================== If less than needed, add zero
                if ( xIt < addXPre ) { newMap[newMapIndex] = 0.0; continue; }
                if ( yIt < addYPre ) { newMap[newMapIndex] = 0.0; continue; }
                if ( zIt < addZPre ) { newMap[newMapIndex] = 0.0; continue; }
                
                //==================================== If more than needed, add zero
                if ( xIt >= ( addXPre + xDimIndices ) ) { newMap[newMapIndex] = 0.0; continue; }
                if ( yIt >= ( addYPre + yDimIndices ) ) { newMap[newMapIndex] = 0.0; continue; }
                if ( zIt >= ( addZPre + zDimIndices ) ) { newMap[newMapIndex] = 0.0; continue; }
                
                //==================================== If not padding, copy original values
                oldMapIndex                           = (zIt-addZPre) + zDimIndices * ( (yIt-addYPre) + yDimIndices * (xIt-addXPre) );
                newMap[newMapIndex]                   = oldMap[oldMapIndex];
            }
        }
    }
    
    //======================================== Done
    return ;
    
}

/*! \brief This function rotates a map based on the given Euler angles.
 
    This function starts by computing the Wigner D matrices for the given Euler angles and then it proceeds to multiply the
    spherical harmonics coefficients with these, thus producing spherical harmonics coefficients of a rotated structure. Then,
    it computes the inverse spherical harmonics decomposition, thus obtaining the sphere mapped values for the rotated structure.
    Finally, it interpolates these sphere mapped values back to Cartesian grid, thus obtaining a map rotated by the given Euler angles.
 
    \param[in] settings The settings object specifying how exactly the rotation is to be done.
    \param[in] eulerAlpha The rotation expressed as a pointer to Euler alpha angle.
    \param[in] eulerBeta The rotation expressed as a pointer to Euler beta angle.
    \param[in] eulerGamma The rotation expressed as a pointer to Euler gamma angle.
 */
void ProSHADE_internal_data::ProSHADE_data::rotateMapReciprocalSpace ( ProSHADE_settings* settings, proshade_double eulerAlpha, proshade_double eulerBeta, proshade_double eulerGamma )
{
    //================================================ Set maximum comparison bandwidth to maximum object bandwidth
    this->maxCompBand                                 = this->spheres[this->noSpheres-1]->getLocalBandwidth();
    
    //================================================ Save map COM after processing but before rotation
    this->findMapCOM                                  ( );
    this->mapCOMProcessChangeX                       += ( this->xCom - this->originalMapXCom );
    this->mapCOMProcessChangeY                       += ( this->yCom - this->originalMapYCom );
    this->mapCOMProcessChangeZ                       += ( this->zCom - this->originalMapZCom );
    
    //================================================ Compute the Wigner D matrices for the Euler angles
    ProSHADE_internal_wigner::computeWignerMatricesForRotation ( settings, this, -eulerAlpha, eulerBeta, -eulerGamma );
    
    //================================================ Initialise rotated Spherical Harmonics memory
    this->allocateRotatedSHMemory                     ( );
    
    //================================================ Multiply SH coeffs by Wigner
    this->computeRotatedSH                            ( );
    
    //================================================ Inverse the SH coeffs to shells
    this->invertSHCoefficients                        ( );
    
    //================================================ Find spherical cut-offs
    std::vector<proshade_double> lonCO, latCO;
    ProSHADE_internal_overlay::computeAngularThreshold ( &lonCO, &latCO, settings->maxBandwidth * 2 );
    
    //================================================ Allocate memory for the rotated map
    proshade_double *densityMapRotated                = new proshade_double [this->xDimIndices * this->yDimIndices * this->zDimIndices];
    ProSHADE_internal_misc::checkMemoryAllocation     ( densityMapRotated, __FILE__, __LINE__, __func__ );
    for ( unsigned int iter = 0; iter < static_cast<unsigned int> ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ ) { densityMapRotated[iter] = 0.0; }
    
    //================================================ Interpolate onto cartesian grid
    this->interpolateMapFromSpheres                   ( densityMapRotated );
    
    //================================================ Copy map
    for ( proshade_unsign iter = 0; iter < ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ )
    {
        this->internalMap[iter]                       = densityMapRotated[iter];
    }
    
    //================================================ Release rotated map (original is now rotated)
    delete[] densityMapRotated;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function rotates a map based on the given angle-axis rotation.
 
    This function takes the axis and angle of the required rotation as well as a pointer to which the rotated map should be
    saved into and proceeds to rotate the map in real space using tri-linear interpolation.
 
    \param[in] axX The x-axis element of the angle-axis rotation representation.
    \param[in] axY The y-axis element of the angle-axis rotation representation.
    \param[in] axZ The z-axis element of the angle-axis rotation representation.
    \param[in] axAng The angle about the axis by which the rotation is to be done.
    \param[in] map A pointer which will be set to point to the rotated map.
    \param[out] ret The rotation centre about which the rotation was done in Angstroms.
 */
std::vector< proshade_double > ProSHADE_internal_data::ProSHADE_data::rotateMapRealSpace ( proshade_double axX, proshade_double axY, proshade_double axZ, proshade_double axAng, proshade_double*& map )
{
    //================================================ Initialise local variables
    bool withinBounds                                 = true;
    proshade_double c000, c001, c010, c011, c100, c101, c110, c111, c00, c01, c10, c11, c0, c1;
    size_t arrPos                                     = 0;
    proshade_double xCOM, yCOM, zCOM;
    std::vector< proshade_double > ret;
    
    //================================================ Store sampling rates
    proshade_single xSampRate                         = this->xDimSize / static_cast< proshade_single > ( this->xTo - this->xFrom );
    proshade_single ySampRate                         = this->yDimSize / static_cast< proshade_single > ( this->yTo - this->yFrom );
    proshade_single zSampRate                         = this->zDimSize / static_cast< proshade_single > ( this->zTo - this->zFrom );
    
    //================================================ Compute map COM
    ProSHADE_internal_mapManip::findMAPCOMValues      ( this->internalMap, &xCOM, &yCOM, &zCOM, this->xDimSize, this->yDimSize, this->zDimSize, this->xFrom, this->xTo, this->yFrom, this->yTo, this->zFrom, this->zTo );
    
    //================================================ Allocate local variables
    proshade_single *mins                             = new proshade_single[3];
    proshade_single *maxs                             = new proshade_single[3];
    proshade_single *rotMat                           = new proshade_single[9];
    proshade_single *rotVec;
    proshade_single *interpMins                       = new proshade_single[3];
    proshade_single *interpMaxs                       = new proshade_single[3];
    proshade_single *interpDiff                       = new proshade_single[3];
    proshade_single *movs                             = new proshade_single[3];
    map                                               = new proshade_double[ this->xDimIndices * this->yDimIndices * this->zDimIndices ];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( mins,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( maxs,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat,     __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( interpMins, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( interpMaxs, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( interpDiff, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( movs,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( map,        __FILE__, __LINE__, __func__ );
    
    //================================================ Fill map with zeroes
    for ( size_t iter = 0; iter < static_cast< size_t > ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ ) { map[iter] = 0.0; }
    
    //================================================ Determine map max's and min's in terms of the hkl system
    mins[0]                                           = std::floor ( static_cast< proshade_single > ( this->xDimIndices ) / -2.0f );
    mins[1]                                           = std::floor ( static_cast< proshade_single > ( this->yDimIndices ) / -2.0f );
    mins[2]                                           = std::floor ( static_cast< proshade_single > ( this->zDimIndices ) / -2.0f );
        
    maxs[0]                                           = -mins[0];
    maxs[1]                                           = -mins[1];
    maxs[2]                                           = -mins[2];
    
    if ( this->xDimIndices % 2 == 0 ) { maxs[0] -= 1.0f; }
    if ( this->yDimIndices % 2 == 0 ) { maxs[1] -= 1.0f; }
    if ( this->zDimIndices % 2 == 0 ) { maxs[2] -= 1.0f; }
    
    //================================================ Rotate about COM instead of map midpoint
    movs[0]                                           = 0.0; //( mins[0] - static_cast< proshade_single > ( this->xFrom ) );
    movs[1]                                           = 0.0; //( mins[1] - static_cast< proshade_single > ( this->yFrom ) );
    movs[2]                                           = 0.0; //( mins[2] - static_cast< proshade_single > ( this->zFrom ) );
    
    //================================================ Save rotation centre
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, static_cast< proshade_double > ( ( -mins[0] * xSampRate ) + ( static_cast< proshade_single > ( this->xFrom ) * xSampRate ) ) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, static_cast< proshade_double > ( ( -mins[1] * ySampRate ) + ( static_cast< proshade_single > ( this->yFrom ) * ySampRate ) ) );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, static_cast< proshade_double > ( ( -mins[2] * zSampRate ) + ( static_cast< proshade_single > ( this->zFrom ) * zSampRate ) ) );
    
    //================================================ Get rotation matrix from Euler angles
    ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( rotMat, axX, axY, axZ, axAng );
    
    //================================================ For each point
    for ( proshade_single xIt = mins[0]; xIt <= maxs[0]; xIt += 1.0f )
    {
        for ( proshade_single yIt = mins[1]; yIt <= maxs[1]; yIt += 1.0f )
        {
            for ( proshade_single zIt = mins[2]; zIt <= maxs[2]; zIt += 1.0f )
            {
                //==================================== Compute new point position
                rotVec                                = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( rotMat, xIt - movs[0], yIt - movs[1], zIt - movs[2] );
                
                //==================================== Find surrounding grid points indices and check for boundaries
                withinBounds                          = true;
                for ( size_t posIt = 0; posIt < 3; posIt++ )
                {
                    //================================ Determine surrounding points indices in this dimension
                    interpMins[posIt]                 = std::floor ( rotVec[posIt] );
                    interpMaxs[posIt]                 = interpMins[posIt] + 1.0f;
                    
                    //================================ Check for boundaries
                    if ( ( maxs[posIt] < interpMins[posIt] ) || ( interpMins[posIt] < mins[posIt] ) || ( maxs[posIt] < interpMaxs[posIt] ) || ( interpMaxs[posIt] < mins[posIt] ) )
                    {
                        withinBounds                  = false;
                        break;
                    }
                    
                    //================================ Compute the difference from position to min index along this axis
                    interpDiff[posIt]                 = rotVec[posIt] - interpMins[posIt];
                }
                if ( !withinBounds ) { continue; }
                
                //==================================== Make sure interpolation max's are within bounds
                for ( size_t posIt = 0; posIt < 3; posIt++ )
                {
                    interpMaxs[posIt]                 = std::min ( maxs[posIt], std::max ( mins[posIt], interpMaxs[posIt] ) );
                }
                
                //==================================== Release memory
                delete[] rotVec;
                
                //==================================== Find the surrounding points values from their indices
                arrPos                                = static_cast< size_t > ( ( interpMins[2] - mins[2] ) + static_cast< proshade_single > ( this->zDimIndices ) * ( ( interpMins[1] - mins[1] ) + static_cast< proshade_single > ( this->yDimIndices ) * ( interpMins[0] - mins[0] ) ) );
                c000                                  = this->internalMap[arrPos];
                
                arrPos                                = static_cast< size_t > ( ( interpMaxs[2] - mins[2] ) + static_cast< proshade_single > ( this->zDimIndices ) * ( ( interpMins[1] - mins[1] ) + static_cast< proshade_single > ( this->yDimIndices ) * ( interpMins[0] - mins[0] ) ) );
                c001                                  = this->internalMap[arrPos];
                
                arrPos                                = static_cast< size_t > ( ( interpMins[2] - mins[2] ) + static_cast< proshade_single > ( this->zDimIndices ) * ( ( interpMaxs[1] - mins[1] ) + static_cast< proshade_single > ( this->yDimIndices ) * ( interpMins[0] - mins[0] ) ) );
                c010                                  = this->internalMap[arrPos];
                
                arrPos                                = static_cast< size_t > ( ( interpMaxs[2] - mins[2] ) + static_cast< proshade_single > ( this->zDimIndices ) * ( ( interpMaxs[1] - mins[1] ) + static_cast< proshade_single > ( this->yDimIndices ) * ( interpMins[0] - mins[0] ) ) );
                c011                                  = this->internalMap[arrPos];
                
                arrPos                                = static_cast< size_t > ( ( interpMins[2] - mins[2] ) + static_cast< proshade_single > ( this->zDimIndices ) * ( ( interpMins[1] - mins[1] ) + static_cast< proshade_single > ( this->yDimIndices ) * ( interpMaxs[0] - mins[0] ) ) );
                c100                                  = this->internalMap[arrPos];
                
                arrPos                                = static_cast< size_t > ( ( interpMaxs[2] - mins[2] ) + static_cast< proshade_single > ( this->zDimIndices ) * ( ( interpMins[1] - mins[1] ) + static_cast< proshade_single > ( this->yDimIndices ) * ( interpMaxs[0] - mins[0] ) ) );
                c101                                  = this->internalMap[arrPos];
                
                arrPos                                = static_cast< size_t > ( ( interpMins[2] - mins[2] ) + static_cast< proshade_single > ( this->zDimIndices ) * ( ( interpMaxs[1] - mins[1] ) + static_cast< proshade_single > ( this->yDimIndices ) * ( interpMaxs[0] - mins[0] ) ) );
                c110                                  = this->internalMap[arrPos];
                
                arrPos                                = static_cast< size_t > ( ( interpMaxs[2] - mins[2] ) + static_cast< proshade_single > ( this->zDimIndices ) * ( ( interpMaxs[1] - mins[1] ) + static_cast< proshade_single > ( this->yDimIndices ) * ( interpMaxs[0] - mins[0] ) ) );
                c111                                  = this->internalMap[arrPos];

                //==================================== Interpolate along x-axis
                c00                                   = ( c000 * ( 1.0 - static_cast< proshade_double > ( interpDiff[0] ) ) ) + ( c100 * static_cast< proshade_double > ( interpDiff[0] ) );
                c01                                   = ( c001 * ( 1.0 - static_cast< proshade_double > ( interpDiff[0] ) ) ) + ( c101 * static_cast< proshade_double > ( interpDiff[0] ) );
                c10                                   = ( c010 * ( 1.0 - static_cast< proshade_double > ( interpDiff[0] ) ) ) + ( c110 * static_cast< proshade_double > ( interpDiff[0] ) );
                c11                                   = ( c011 * ( 1.0 - static_cast< proshade_double > ( interpDiff[0] ) ) ) + ( c111 * static_cast< proshade_double > ( interpDiff[0] ) );
                
                //==================================== Interpolate along y-axis
                c0                                    = ( c00  * ( 1.0 - static_cast< proshade_double > ( interpDiff[1] ) ) ) + ( c10  * static_cast< proshade_double > ( interpDiff[1] ) );
                c1                                    = ( c01  * ( 1.0 - static_cast< proshade_double > ( interpDiff[1] ) ) ) + ( c11  * static_cast< proshade_double > ( interpDiff[1] ) );

                //==================================== Interpolate along z-axis
                arrPos                                = static_cast< size_t > ( ( zIt - mins[2] ) + static_cast< proshade_single > ( this->zDimIndices ) * ( ( yIt - mins[1] ) + static_cast< proshade_single > ( this->yDimIndices ) * ( xIt - mins[0] ) ) );
                map[arrPos]                           = ( c0   * ( 1.0 - static_cast< proshade_double > ( interpDiff[2] ) ) ) + ( c1   * static_cast< proshade_double > ( interpDiff[2] ) );
            }
        }
    }
    
    //================================================ Release memory
    delete[] mins;
    delete[] maxs;
    delete[] rotMat;
    delete[] interpMins;
    delete[] interpMaxs;
    delete[] interpDiff;
    delete[] movs;
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function rotates a map based on the given Euler angles in place.
 
    This function takes the Euler angles of the required rotation and proceeds to make use of the
    rotateMapRealSpace () function to rotate the map in real space using trilinear interpolation,
    replacing the original map with the rotated one.
 
    \param[in] eulerAlpha The rotation expressed as a pointer to Euler alpha angle.
    \param[in] eulerBeta The rotation expressed as a pointer to Euler beta angle.
    \param[in] eulerGamma The rotation expressed as a pointer to Euler gamma angle.
    \param[out] ret The rotation centre about which the rotation was done in Angstroms.
 */
std::vector< proshade_double > ProSHADE_internal_data::ProSHADE_data::rotateMapRealSpaceInPlace ( proshade_double eulA, proshade_double eulB, proshade_double eulG )
{
    //================================================ Initialise local variables
    proshade_double axX, axY, axZ, axAng, tmp, *rMat, *map;
    
    //================================================ Allocate local memory
    rMat =                                            new proshade_double[9];
    
    //================================================ Check local memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( rMat, __FILE__, __LINE__, __func__ );
    
    //================================================ Convert Euler angles to rotation matrix
    ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( eulG, eulB, eulA, rMat );
    
    //================================================ Transpose the rotation matrix
    tmp     = rMat[1];
    rMat[1] = rMat[3];
    rMat[3] = tmp;
    
    tmp     = rMat[2];
    rMat[2] = rMat[6];
    rMat[6] = tmp;
    
    tmp     = rMat[5];
    rMat[5] = rMat[7];
    rMat[7] = tmp;
 
    //================================================ Convert rotation matrix to angle-axis
    ProSHADE_internal_maths::getAxisAngleFromRotationMatrix ( rMat, &axX, &axY, &axZ, &axAng );
    
    //================================================ Rotate the internal map
    std::vector< proshade_double > ret                = this->rotateMapRealSpace ( axX, axY, axZ, axAng, map );
    
    //================================================ Copy the rotated map in place of the internal map
    for ( size_t iter = 0; iter < static_cast< size_t > ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ )
    {
        this->internalMap[iter]                       = map[iter];
    }
    
    //================================================ Release local memory
    delete[] rMat;
    delete[] map;
    
    //================================================ Save rotation centre for co-ordinates
    this->originalPdbRotCenX                          = ret.at(0);
    this->originalPdbRotCenY                          = ret.at(1);
    this->originalPdbRotCenZ                          = ret.at(2);
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function simply translates the map by a given number of Angstroms along the three axes.
 
    This function calls the internal functions to first provide the maximum possible movement by changing the frame
    of the map and secondly, it make the precise movement within this new frame using the Fourier translation approach.
 
    \param[in] settings The settings object specifying how exactly the rotation is to be done.
    \param[in] trsX The translation expressed as a number of angstroms to move by along the x-axis.
    \param[in] trsY The translation expressed as a number of angstroms to move by along the y-axis.
    \param[in] trsZ The translation expressed as a number of angstroms to move by along the z-axis.
 */
void ProSHADE_internal_data::ProSHADE_data::translateMap ( proshade_double trsX, proshade_double trsY, proshade_double trsZ )
{
    //================================================ Initialise local variables
    proshade_single xMov                              = static_cast< proshade_single > ( -trsX );
    proshade_single yMov                              = static_cast< proshade_single > ( -trsY );
    proshade_single zMov                              = static_cast< proshade_single > ( -trsZ );
    
    //================================================ Move the whole map frame to minimise the Fourier movement
    ProSHADE_internal_mapManip::moveMapByIndices      ( &xMov, &yMov, &zMov, this->getXDimSize(), this->getYDimSize(), this->getZDimSize(),
                                                        this->getXFromPtr(), this->getXToPtr(), this->getYFromPtr(), this->getYToPtr(),
                                                        this->getZFromPtr(), this->getZToPtr(), this->getXAxisOrigin(), this->getYAxisOrigin(), this->getZAxisOrigin() );
    
    //================================================ Finalise the movement by in-frame Fourier movement
    ProSHADE_internal_mapManip::moveMapByFourier      ( this->getInternalMap(), xMov, yMov, zMov, this->getXDimSize(), this->getYDimSize(), this->getZDimSize(),
                                                        static_cast< proshade_signed > ( this->getXDim() ), static_cast< proshade_signed > ( this->getYDim() ),
                                                        static_cast< proshade_signed > ( this->getZDim() ) );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allocates the memory required for storing the rotated Spherical Harmonics coefficients.
 */
void ProSHADE_internal_data::ProSHADE_data::allocateRotatedSHMemory ( )
{
    //================================================ Allocate the main pointer and check
    this->rotSphericalHarmonics                       = new proshade_complex* [this->noSpheres];
    ProSHADE_internal_misc::checkMemoryAllocation ( this->rotSphericalHarmonics, __FILE__, __LINE__, __func__ );
    
    //================================================ For each sphere
    for ( proshade_unsign iter = 0; iter < this->noSpheres; iter++ )
    {
        //============================================ Allocate the sphere storage space
        this->rotSphericalHarmonics[iter]             = new proshade_complex [static_cast<proshade_unsign> ( pow ( (this->spheres[iter]->getLocalBandwidth() * 2), 2) )];
        ProSHADE_internal_misc::checkMemoryAllocation ( this->rotSphericalHarmonics[iter],  __FILE__, __LINE__, __func__ );
        
        //============================================ Set values to zeroes
        for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( pow ( (this->spheres[iter]->getLocalBandwidth() * 2), 2) ); it++ )
        {
            this->rotSphericalHarmonics[iter][it][0]  = 0.0;
            this->rotSphericalHarmonics[iter][it][1]  = 0.0;
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function multiplies the objects spherical harmonics with the Wigner D matrices, obtaining rotated spherical harmonics coefficients.
 */
void ProSHADE_internal_data::ProSHADE_data::computeRotatedSH ( )
{
    //================================================ Initialise variables
    proshade_double WigDR, WigDI, *ShR, *ShI, retR, retI;
    proshade_unsign arrPos;
    
    //================================================ Compute
    for ( proshade_signed shell = 0; shell < static_cast<proshade_signed> ( this->noSpheres ); shell++ )
    {
        //============================================ For each band
        for ( proshade_signed bandIter = 0; bandIter < static_cast<proshade_signed> ( this->spheres[shell]->getLocalBandwidth() ); bandIter++ )
        {
            //======================================== For each order1
            for ( proshade_signed order1 = 0; order1 < ( ( bandIter * 2 ) + 1 ); order1++ )
            {
                //==================================== Get Spherical Harmonics value
                ShR                                   = getRealSphHarmValue  ( static_cast< proshade_unsign > ( bandIter ), static_cast< proshade_unsign > ( order1 ), static_cast< proshade_unsign > ( shell ) );
                ShI                                   = getImagSphHarmValue  ( static_cast< proshade_unsign > ( bandIter ), static_cast< proshade_unsign > ( order1 ), static_cast< proshade_unsign > ( shell ) );
                
                //==================================== For each order2
                for ( proshade_signed order2 = 0; order2 < ( ( bandIter * 2 ) + 1 ); order2++ )
                {
                    //================================ Get Wigner D values
                    this->getWignerMatrixValue        ( static_cast< proshade_unsign > ( bandIter ), static_cast< proshade_unsign > ( order1 ), static_cast< proshade_unsign > ( order2 ), &WigDR, &WigDI );
                    
                    //================================ Multiply SH and Wigner
                    ProSHADE_internal_maths::complexMultiplication ( ShR, ShI, &WigDR, &WigDI, &retR, &retI );

                    //================================ Save
                    arrPos                            = static_cast<proshade_unsign> ( seanindex ( static_cast< int > ( order2-bandIter ), static_cast< int > ( bandIter ),
                                                                                                   static_cast< int > ( this->spheres[shell]->getLocalBandwidth() ) ) );
                    this->rotSphericalHarmonics[shell][arrPos][0] += retR;
                    this->rotSphericalHarmonics[shell][arrPos][1] += retI;
                }
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function initialises internal variables for inverse Spherical Harmonics computation.
 
    \param[in] shBand The bandwidth for this particular shell.
    \param[in] sigR Pointer to be initialised for the real signal values.
    \param[in] sigI Pointer to be initialised for the imaginary signal values.
    \param[in] rcoeffs Pointer to be initialised for the real coefficient values.
    \param[in] icoeffs Pointer to be initialised for the imaginary coefficient values.
    \param[in] weights Pointer to be initialised for the transform weight values.
    \param[in] workspace Pointer to be initialised for the computation screatch space.
    \param[in] idctPlan Pointer reference to the cosine/sine transform plan to be created.
    \param[in] ifftPlan Pointer reference to the discrete 3D Fourier transform plan to be created.
 */
void ProSHADE_internal_overlay::initialiseInverseSHComputation ( proshade_unsign shBand, double*& sigR, double*& sigI, double*& rcoeffs, double*& icoeffs, double*& weights, double*& workspace, fftw_plan& idctPlan, fftw_plan& ifftPlan )
{
    //================================================ Initialise internal variables
    proshade_unsign oneDim                            = shBand * 2;
    
    //================================================ Allocate memory
    sigR                                              = new double [(oneDim * oneDim * 4)];
    sigI                                              = sigR + (oneDim * oneDim * 2);
    rcoeffs                                           = new double [(oneDim * oneDim * 2)];
    icoeffs                                           = rcoeffs + (oneDim * oneDim);
    weights                                           = new double [4 * shBand];
    workspace                                         = new double [( 20 * shBand * shBand ) + ( 24 * shBand )];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( sigR,      __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( rcoeffs,   __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( weights,   __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( workspace, __FILE__, __LINE__, __func__ );
    
    //================================================ Create the cosine/sine transform plan
    idctPlan                                          = fftw_plan_r2r_1d ( static_cast< int > ( oneDim ), weights, workspace, FFTW_REDFT01, FFTW_ESTIMATE );
    
    //================================================ Set up the discrete Fourier transform
    int rank, howmany_rank;
    fftw_iodim dims[1], howmany_dims[1];
    
    rank                                              = 1;
    dims[0].n                                         = 2 * static_cast< int > ( shBand );
    dims[0].is                                        = 2 * static_cast< int > ( shBand );
    dims[0].os                                        = 1;
    howmany_rank                                      = 1;
    howmany_dims[0].n                                 = 2 * static_cast< int > ( shBand );
    howmany_dims[0].is                                = 1;
    howmany_dims[0].os                                = 2 * static_cast< int > ( shBand );
    
    //================================================ Create the discrete Fourier transform
    ifftPlan                                          = fftw_plan_guru_split_dft ( rank, dims, howmany_rank, howmany_dims, sigR, sigI, rcoeffs, icoeffs, FFTW_ESTIMATE );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the shell mapped data from inverting the Spherical Harmonics coefficients.
 
 */
void ProSHADE_internal_data::ProSHADE_data::invertSHCoefficients ( )
{
    //================================================ Initialise local variables
    double *sigR =  nullptr, *sigI =  nullptr, *rcoeffs =  nullptr, *icoeffs =  nullptr, *weights =  nullptr, *workspace =  nullptr;
    fftw_plan idctPlan, ifftPlan;
    
    //================================================ For each shell
     for ( int shell = 0; shell < static_cast<int> ( this->noSpheres ); shell++ )
     {
         //=========================================== Initialise internal variables
         proshade_unsign oneDim                       = this->spheres[shell]->getLocalBandwidth() * 2;
         
         //=========================================== Allocate memory
         ProSHADE_internal_overlay::initialiseInverseSHComputation ( this->spheres[shell]->getLocalBandwidth(), sigR, sigI, rcoeffs, icoeffs, weights, workspace, idctPlan, ifftPlan );
         
         //=========================================== Compute weights for the transform using the appropriate shell related band
         makeweights                                  ( static_cast< int > ( this->spheres[shell]->getLocalBandwidth() ), weights );
         
         //============================================ Allocate rotated shell mapped data memory
         this->spheres[shell]->allocateRotatedMap     ( );
         
         //============================================ Load SH coeffs to arrays
         for ( unsigned int iter = 0; iter < static_cast<unsigned int> ( oneDim * oneDim ); iter++ )
         {
             rcoeffs[iter]                            = this->rotSphericalHarmonics[shell][iter][0];
             icoeffs[iter]                            = this->rotSphericalHarmonics[shell][iter][1];
             sigR[iter]                               = 0.0;
             sigI[iter]                               = 0.0;
         }
         
         //============================================ Get inverse spherical harmonics transform for the shell
         InvFST_semi_fly                              ( rcoeffs,
                                                        icoeffs,
                                                        sigR,
                                                        sigI,
                                                        static_cast< int > ( this->spheres[shell]->getLocalBandwidth() ),
                                                        workspace,
                                                        0,
                                                        static_cast< int > ( this->spheres[shell]->getLocalBandwidth() ),
                                                        &idctPlan,
                                                        &ifftPlan );
         
         //=========================================== Copy the results to the rotated shells array
         for ( unsigned int iter = 0; iter < static_cast<unsigned int> ( oneDim * oneDim ); iter++ )
         {
             this->spheres[shell]->setRotatedMappedData ( iter, sigR[iter] );
         }
         
         //=========================================== Release the plans
         fftw_destroy_plan ( idctPlan );
         fftw_destroy_plan ( ifftPlan );
         
         //=========================================== Release the memory
         delete[] sigR;
         delete[] rcoeffs;
         delete[] weights;
         delete[] workspace;
     }
    
     //================================================ Done
     return ;
     
}

/*! \brief This function computes the angular thresholds for longitude and lattitude angles.

    \param[in] lonCO Pointer to vector where longitude thresholds are to be stored at.
    \param[in] latCO Pointer to vector where lattitude thresholds are to be stored at.
    \param[in] angRes The angular resolution of the computation.
*/
void ProSHADE_internal_overlay::computeAngularThreshold ( std::vector<proshade_double>* lonCO, std::vector<proshade_double>* latCO, proshade_unsign angRes )
{
    //================================================ Longitude angular thresholds
    for ( proshade_unsign iter = 0; iter <= angRes; iter++ )
    {
        ProSHADE_internal_misc::addToDoubleVector ( lonCO, static_cast<proshade_double> ( iter ) * ( ( static_cast<proshade_double> ( M_PI ) * 2.0 ) / static_cast<proshade_double> ( angRes ) ) - ( static_cast<double> ( M_PI ) ) );
    }
    
    //================================================ Lattitude angular thresholds
    for ( proshade_unsign iter = 0; iter <= angRes; iter++ )
    {
        ProSHADE_internal_misc::addToDoubleVector ( latCO, static_cast<proshade_double> ( iter ) * ( static_cast<proshade_double> ( M_PI ) / static_cast<proshade_double> ( angRes ) ) - ( static_cast<proshade_double> ( M_PI ) / 2.0 ) );
    }
    
    //================================================ Done
    return ;
     
}

/*! \brief This function interpolates the density map from the sphere mapped data.
 
    \param[in] densityMapRotated The pointer to allocated memory where the new map values will be held.
 */
void ProSHADE_internal_data::ProSHADE_data::interpolateMapFromSpheres ( proshade_double*& densityMapRotated )
{
    //================================================ Initialise variables
    proshade_double rad = 0.0, lon = 0.0, lat = 0.0, newU = 0.0, newV = 0.0, newW = 0.0;
    proshade_unsign lowerLonL = 0, upperLonL = 0, lowerLonU = 0, upperLonU = 0, lowerLatL = 0, upperLatL = 0, lowerLatU = 0, upperLatU = 0, lowerShell = 0, upperShell = 0;
    proshade_double x00 = 0.0, x01 = 0.0, x10 = 0.0, x11 = 0.0, distLLon = 0.0, distLLat = 0.0, distLRad = 0.0, valLLon = 0.0, valULon = 0.0;
    proshade_double lowerShellValue = 0.0, upperShellValue = 0.0;
    proshade_double xSamplingRate                     = static_cast<proshade_double> ( this->xDimSize ) / static_cast<proshade_double> ( this->xDimIndices );
    proshade_double ySamplingRate                     = static_cast<proshade_double> ( this->yDimSize ) / static_cast<proshade_double> ( this->yDimIndices );
    proshade_double zSamplingRate                     = static_cast<proshade_double> ( this->zDimSize ) / static_cast<proshade_double> ( this->zDimIndices );
    proshade_signed arrPos;
    std::vector<proshade_double> lonCOU, latCOU, lonCOL, latCOL;
    
    for ( proshade_signed uIt = 0; uIt < static_cast<proshade_signed> (this->xDimIndices); uIt++ )
    {
        for ( proshade_signed vIt = 0; vIt < static_cast<proshade_signed> (this->yDimIndices); vIt++ )
        {
            for ( proshade_signed wIt = 0; wIt < static_cast<proshade_signed> (this->zDimIndices); wIt++ )
            {
                //==================================== Convert to centered coords
                newU                                  = static_cast<proshade_double> ( uIt - ( static_cast<proshade_signed> (this->xDimIndices) / 2 ) );
                newV                                  = static_cast<proshade_double> ( vIt - ( static_cast<proshade_signed> (this->yDimIndices) / 2 ) );
                newW                                  = static_cast<proshade_double> ( wIt - ( static_cast<proshade_signed> (this->zDimIndices) / 2 ) );
                
                //==================================== Deal with 0 ; 0 ; 0
                if ( ( newU == 0.0 ) && ( newV == 0.0 ) && ( newW == 0.0 ) )
                {
                    arrPos                            = wIt + static_cast< proshade_signed > ( this->zDimIndices ) * ( vIt + static_cast< proshade_signed > ( this->yDimIndices ) * uIt );
                    densityMapRotated[arrPos]         = this->internalMap[arrPos];
                    continue;
                }
                
                //==================================== Convert to spherical coords
                rad                                   = sqrt  ( pow( ( newU * xSamplingRate ), 2.0 ) +
                                                                pow( ( newV * ySamplingRate ), 2.0 ) +
                                                                pow( ( newW * zSamplingRate ), 2.0 ) );
                lon                                   = atan2 ( ( newV * ySamplingRate ), ( newU * xSamplingRate ) );
                lat                                   = asin  ( ( newW * zSamplingRate ) / rad );
                
                //==================================== Deal with nan's
                if ( rad   != rad ) { rad   = 0.0; }
                if ( lon   != lon ) { lon   = 0.0; }
                if ( lat   != lat ) { lat   = 0.0; }
                                
                //==================================== Find shells above and below
                lowerShell                            = 0;
                upperShell                            = 0;
                for ( proshade_unsign iter = 0; iter < (this->noSpheres-1); iter++ )
                {
                    if ( ( static_cast< proshade_double > ( this->spherePos.at(iter) ) <= rad ) && ( static_cast< proshade_double > ( this->spherePos.at(iter+1) ) > rad ) )
                    {
                        lowerShell                    = iter;
                        upperShell                    = iter+1;
                        break;
                    }
                }
                
                if ( upperShell == 0 )
                {
                    arrPos                            = wIt + static_cast< proshade_signed > ( this->zDimIndices ) * ( vIt + static_cast< proshade_signed > ( this->yDimIndices ) * uIt );
                    densityMapRotated[arrPos]         = 0.0;
                    continue;
                }
                
                //==================================== Get the longitude and lattitude cut-offs for this shell resolution
                lonCOL.clear(); latCOL.clear(); lonCOU.clear(); latCOU.clear();
                ProSHADE_internal_overlay::computeAngularThreshold ( &lonCOL, &latCOL, this->spheres[lowerShell]->getLocalAngRes() );
                ProSHADE_internal_overlay::computeAngularThreshold ( &lonCOU, &latCOU, this->spheres[upperShell]->getLocalAngRes() );
                
                //==================================== Find the angle cutoffs around the point
                for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( lonCOL.size() ); iter++ )
                {
                    if ( iter == ( static_cast<proshade_unsign> ( lonCOL.size() ) - 1 ) )
                    {
                        lowerLonL                     = 0;
                        upperLonL                     = 1;
                        break;
                    }
                    if ( ( std::floor(10000. * lonCOL.at(iter)) <= std::floor(10000. * lon) ) && ( std::floor(10000. * lonCOL.at(iter+1)) > std::floor(10000. * lon) ) )
                    {
                        lowerLonL                     = iter;
                        upperLonL                     = iter+1;
                        break;
                    }
                }
                if ( upperLonL == this->spheres[lowerShell]->getLocalAngRes() ) { upperLonL = 0; }
                
                for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( lonCOU.size() ); iter++ )
                {
                    if ( iter == ( static_cast<proshade_unsign> ( lonCOU.size() ) - 1 ) )
                    {
                        lowerLonU                     = 0;
                        upperLonU                     = 1;
                        break;
                    }
                    if ( ( std::floor(10000. * lonCOU.at(iter)) <= std::floor(10000. * lon) ) && ( std::floor(10000. * lonCOU.at(iter+1)) > std::floor(10000. * lon) ) )
                    {
                        lowerLonU                     = iter;
                        upperLonU                     = iter+1;
                        break;
                    }
                }
                if ( upperLonU == this->spheres[upperShell]->getLocalAngRes() ) { upperLonU = 0; }
                
                for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( latCOL.size() ); iter++ )
                {
                    if ( iter == ( static_cast<proshade_unsign> ( latCOL.size() ) - 1 ) )
                    {
                        lowerLatL                     = 0;
                        upperLatL                     = 1;
                        break;
                    }
                    if ( ( std::floor(10000. * latCOL.at(iter)) <=  std::floor(10000. * lat) ) && ( std::floor(10000. * latCOL.at(iter+1)) > std::floor(10000. * lat) ) )
                    {
                        lowerLatL                     = iter;
                        upperLatL                     = iter+1;
                        break;
                    }
                }
                if ( upperLatL == this->spheres[lowerShell]->getLocalAngRes() ) { upperLatL = 0; }
                
                for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( latCOU.size() ); iter++ )
                {
                    if ( iter == ( static_cast<proshade_unsign> ( latCOU.size() ) - 1 ) )
                    {
                        lowerLatU                     = 0;
                        upperLatU                     = 1;
                        break;
                    }
                    if ( ( std::floor(10000. * latCOU.at(iter)) <=  std::floor(10000. * lat) ) && ( std::floor(10000. * latCOU.at(iter+1)) > std::floor(10000. * lat) ) )
                    {
                        lowerLatU                     = iter;
                        upperLatU                     = iter+1;
                        break;
                    }
                }
                if ( upperLatU == this->spheres[upperShell]->getLocalAngRes() ) { upperLatU = 0; }
                
                //==================================== Interpolate lower shell
                x00                                   = this->spheres[lowerShell]->getRotatedMappedData ( lowerLatL * this->spheres[lowerShell]->getLocalAngRes() + lowerLonL );
                x01                                   = this->spheres[lowerShell]->getRotatedMappedData ( lowerLatL * this->spheres[lowerShell]->getLocalAngRes() + upperLonL );
                x10                                   = this->spheres[lowerShell]->getRotatedMappedData ( upperLatL * this->spheres[lowerShell]->getLocalAngRes() + lowerLonL );
                x11                                   = this->spheres[lowerShell]->getRotatedMappedData ( upperLatL * this->spheres[lowerShell]->getLocalAngRes() + upperLonL );
                        
                distLLon                              = std::abs ( lon - lonCOL.at(lowerLonL) ) / ( std::abs( lon - lonCOL.at(lowerLonL) ) + std::abs( lon - lonCOL.at(upperLonL) ) );
                valLLon                               = ( ( 1.0 - distLLon ) * x00 ) + ( distLLon * x01 );
                valULon                               = ( ( 1.0 - distLLon ) * x10 ) + ( distLLon * x11 );
                        
                distLLat                              = std::abs ( lat - latCOL.at(lowerLatL) ) / ( std::abs( lat - latCOL.at(lowerLatL) ) + std::abs( lat - latCOL.at(upperLatL) ) );
                lowerShellValue                       = ( ( 1.0 - distLLat ) * valLLon ) + ( distLLat * valULon );
                
                //==================================== Interpolate upper shell
                x00                                   = this->spheres[upperShell]->getRotatedMappedData ( lowerLatU * this->spheres[upperShell]->getLocalAngRes() + lowerLonU );
                x01                                   = this->spheres[upperShell]->getRotatedMappedData ( lowerLatU * this->spheres[upperShell]->getLocalAngRes() + upperLonU );
                x10                                   = this->spheres[upperShell]->getRotatedMappedData ( upperLatU * this->spheres[upperShell]->getLocalAngRes() + lowerLonU );
                x11                                   = this->spheres[upperShell]->getRotatedMappedData ( upperLatU * this->spheres[upperShell]->getLocalAngRes() + upperLonU );
                        
                distLLon                              = std::abs ( lon - lonCOU.at(lowerLonU) ) / ( std::abs( lon - lonCOU.at(lowerLonU) ) + std::abs( lon - lonCOU.at(upperLonU) ) );
                valLLon                               = ( ( 1.0 - distLLon ) * x00 ) + ( distLLon * x01 );
                valULon                               = ( ( 1.0 - distLLon ) * x10 ) + ( distLLon * x11 );
                        
                distLLat                              = std::abs ( lat - latCOU.at(lowerLatU) ) / ( std::abs( lat - latCOU.at(lowerLatU) ) + std::abs( lat - latCOU.at(upperLatU) ) );
                upperShellValue                       = ( ( 1.0 - distLLat ) * valLLon ) + ( distLLat * valULon );
                
                //==================================== Interpolate between shells
                distLRad                              = std::abs ( rad - static_cast< proshade_double > ( this->spherePos.at(lowerShell) ) ) / ( std::abs( rad - static_cast< proshade_double > ( this->spherePos.at(lowerShell) ) ) +
                                                                                                              std::abs( rad - static_cast< proshade_double > ( this->spherePos.at(upperShell) ) ) );
                
                arrPos                                = wIt + static_cast< proshade_signed > ( this->zDimIndices ) * ( vIt + static_cast< proshade_signed > ( this->yDimIndices ) * uIt );
                densityMapRotated[arrPos]             = ( ( 1.0 - distLRad ) * lowerShellValue ) + ( distLRad * upperShellValue );
            }
     
        }

    }
    
    //================================================ Done
    return ;
     
}

/*! \brief This function returns a vector of three floats, the three Euler angles of the best peak in the rotation map.

    \param[in] settings A pointer to settings class containing all the information required for map overlay.
    \param[out] val A vector of the Euler angles of the best peak in the rotation function map.
*/
std::vector< proshade_double > ProSHADE_internal_data::ProSHADE_data::getBestRotationMapPeaksEulerAngles ( ProSHADE_settings* settings )
{
    //================================================ Initialise local variables
    std::vector < proshade_double > ret;
    proshade_double eulA, eulB, eulG;
    
    //================================================ Get inverse SO(3) map top peak Euler angle values
    ProSHADE_internal_peakSearch::getBestPeakEulerAngsNaive ( this->getInvSO3Coeffs (),
                                                              this->getMaxBand() * 2,
                                                             &eulA, &eulB, &eulG, settings );
    
    //================================================ Re-format to vector
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, eulA );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, eulB );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, eulG );
     
    //================================================ Done
    return                                            ( ret );
    
}
