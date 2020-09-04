/*! \file ProSHADE.hpp
    \brief This is the main header file providing the main access class and its functions.

    This file contains the header declarations for the main access class (ProSHADE_run), which the user can use to
    start a run of the ProSHADE tool. This is generally done by firstly creating and setting a ProSHADE_settings class
    instance, which is then supplied to the constructor of the ProSHADE_run class defined here. Once the class constructor
    is run, ProSHADE run with the supplied settings will also be complete and the resulting instance of ProSHADE_run class
    will contain the results. To access these results, the ProSHADE_run class provides accessor functions.

    Copyright by Michal Tykac and individual contributors. All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
 
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.4
    \date      SEP 2020
 */

//==================================================== ProSHADE library code
#include "ProSHADE_tasks.hpp"

//==================================================== Overinclusion protection
#ifndef __PROSHADE__
#define __PROSHADE__
 
/*! \class ProSHADE_run
    \brief This class provides the access point to the library.
 
    This class codes the object that the user of the library needs to create (and presumably delete) in order to get access to the ProSHADE library.   
 */

class ProSHADE_run
{
private:
    //================================================ General variables
    proshade_unsign noStructures;                     //!< This variable hold the total number of structures used. This is needed for the python numpy calls.
    proshade_signed verbose;                          //!< This variable holds the verbose value for the object.
    
    //================================================ Variables regarding distances computation
    std::vector < proshade_double > enLevs;           //!< Vector holding energy levels distances from the first to all other supplied structures.
    std::vector < proshade_double > trSigm;           //!< Vector holding trace sigma distances from the first to all other supplied structures.
    std::vector < proshade_double > rotFun;           //!< Vector holding full rotation function distances from the first to all other supplied structures.
    
    //================================================ Variables regarding symmetry detection
    std::vector< proshade_double* > RecomSymAxes;     //!< Vector holding the recommended symmetry axes information.
    
    //================================================ Variables regarding re-boxing task
    std::vector < proshade_signed* > originalBounds;  //!< Original boundaries of the map.
    std::vector < proshade_signed* > reboxedBounds;   //!< Re-boxed boundaries of the map.
    std::vector < proshade_double* > manipulatedMaps; //!< The map (in XYZ format) after all manipulations are done. It will have the dimensions of reboxedBounds, but the rest of map information will not be available in the simpleAccess ProSHADE run.
    
    //================================================ Variables regarding overlay optimisation
    std::vector < proshade_double > eulerAngles;      //!< Vector of three Euler angles (ZXZ convention) specifying the rotation required to best overlay two structures.
    std::vector < proshade_double > translation;      //!< Vector of three translation vectors specifying the translation required to best overlat two structures.
    
    //================================================ Variables regarding symmetry detection
    std::string symRecommType;                        //!< The resulting recommended symmetry type for the symmetry detection task.
    proshade_unsign symRecommFold;                    //!< The resulting recommended symmetry fold foe the symmetry detection task.
    
private:
    //================================================ Mutator functions
    void setRecommendedSymmetry                       ( std::string val );
    void setRecommendedFold                           ( proshade_unsign val );
    
    //================================================ Task completion functions
    void setSymmetryResults                           ( ProSHADE_settings* settings );
    
public:
    //================================================ Constructors / Destructors
    ProSHADE_run                                      ( ProSHADE_settings* settings );
   ~ProSHADE_run                                      ( void );
    
public:
    //================================================ General accessor functions
    proshade_unsign getNoStructures                   ( void );
    proshade_signed getVerbose                        ( void );
    
    //================================================ Distances accessor functions
    proshade_double getEnergyLevelsVectorValue        ( proshade_unsign pos = 0 );
    proshade_unsign getEnergyLevelsLength             ( void );
    proshade_double getTraceSigmaVectorValue          ( proshade_unsign pos = 0 );
    proshade_unsign getTraceSigmaLength               ( void );
    proshade_double getRotationFunctionVectorValue    ( proshade_unsign pos = 0 );
    proshade_unsign getRotationFunctionLength         ( void );
    
    //================================================ Symmetry accessor functions
    proshade_unsign getNoSymmetryAxes                 ( void );
    
public:
    //================================================ Distances results accessor functions
    std::vector< proshade_double > getEnergyLevelsVector ( void );
    std::vector< proshade_double > getTraceSigmaVector ( void );
    std::vector< proshade_double > getRotationFunctionVector  ( void );

    //================================================ Symmetry results accessor functions
    std::string getSymmetryType                       ( void );
    proshade_unsign getSymmetryFold                   ( void );
    std::vector< std::string > getSymmetryAxis        ( proshade_unsign axisNo );
    
    //================================================ Re-boxing results accessor functions
    std::vector< proshade_signed > getOriginalBounds  ( proshade_unsign strNo );
    std::vector< proshade_signed > getReBoxedBounds   ( proshade_unsign strNo );
    proshade_double getMapValue                       ( proshade_unsign strNo, proshade_unsign mapIndex );
    
    //================================================ Overlay results accessor functions
    std::vector< proshade_double > getEulerAngles     ( void );
    std::vector< proshade_double > getOptimalRotMat   ( void );
    std::vector< proshade_double > getTranslation     ( void );
};

//==================================================== These functions should be in ProSHADE_run class, but I cannot make them work with Numpy from there, so they are here.
void getEnergyLevelsVectorNumpy                       ( ProSHADE_run* run, int verbose, double *enLevVec, int len );
void getTraceSigmaVectorNumpy                         ( ProSHADE_run* run, int verbose, double *trSigVec, int len );
void getRotationFunctionVectorNumpy                   ( ProSHADE_run* run, int verbose, double *rotFnVec, int len );
                  
void getOriginalBoundsVectorNumpy                     ( ProSHADE_run* run, proshade_unsign strNo, int *boundsVec,   int len );
void getReBoxedBoundsVectorNumpy                      ( ProSHADE_run* run, proshade_unsign strNo, int *reboxVec,    int len );
void getReBoxedMap                                    ( ProSHADE_run* run, proshade_unsign strNo, double *reboxMap, int len );
                  
void getOptimalEulerAngles                            ( ProSHADE_run* run, double *eulerAngs, int len );
void getOptimalTranslation                            ( ProSHADE_run* run, double *translate, int len );

#endif
