/*! \file ProSHADE.hpp
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

//============================================ ProSHADE library code
#include "ProSHADE_tasks.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE__
#define __PROSHADE__
 
/*! \class ProSHADE_run
    \brief This class provides the access point to the library.
 
 This class codes the object that the user of the library needs to create (and presumably delete) in order to get access to the ProSHADE library.
 */

class ProSHADE_run
{
private:
    //======================================== General variables
    proshade_unsign noStructures;             //!< This variable hold the total number of structures used. This is needed for the python numpy calls.
    proshade_signed verbose;                  //!< This variable holds the verbose value for the object.
    
    //======================================== Variables regarding distances computation
    std::vector < proshade_double > enLevs;   //!< Vector holding energy levels distances from the first to all other supplied structures.
    std::vector < proshade_double > trSigm;   //!< Vector holding trace sigma distances from the first to all other supplied structures.
    std::vector < proshade_double > rotFun;   //!< Vector holding full rotation function distances from the first to all other supplied structures.
    
    //======================================== Variables regarding symmetry detection
    std::vector< proshade_double* > RecomSymAxes; //!< Vector holding the recommended symmetry axes information.
    
    //======================================== Variables regarding re-boxing task
    std::vector < proshade_signed* > originalBounds; //!< Original boundaries of the map.
    std::vector < proshade_signed* > reboxedBounds;  //!< Re-boxed boundaries of the map.
    std::vector < proshade_double* > manipulatedMaps; //!< The map (in XYZ format) after all manipulations are done. It will have the dimensions of reboxedBounds, but the rest of map information will not be available in the simpleAccess ProSHADE run.
    
    //======================================== Variables regarding overlay optimisation
    std::vector < proshade_double > eulerAngles; //!< Vector of three Euler angles (ZXZ convention) specifying the rotation required to best overlay two structures.
    std::vector < proshade_double > translation; //!< Vector of three translation vectors specifying the translation required to best overlat two structures.
    
    //======================================== Variables regarding symmetry detection
    std::string symRecommType;                //!< The resulting recommended symmetry type for the symmetry detection task.
    proshade_unsign symRecommFold;            //!< The resulting recommended symmetry fold foe the symmetry detection task.
    
private:
    //======================================== Mutator functions
    void                                      setRecommendedSymmetry                ( std::string val );
    void                                      setRecommendedFold                    ( proshade_unsign val );
    
    //======================================== Task completion functions
    void                                      setSymmetryResults                    ( ProSHADE_settings* settings );
    
public:
    //======================================== Constructors / Destructors
    ProSHADE_run                                                                    ( ProSHADE_settings* settings );
   ~ProSHADE_run                                                                    ( void );
    
public:
    //======================================== General accessor functions
    proshade_unsign                           getNoStructures                       ( void );
    proshade_signed                           getVerbose                            ( void );
    
    //======================================== Distances accessor functions
    proshade_double                           getEnergyLevelsVectorValue            ( proshade_unsign pos = 0 );
    proshade_unsign                           getEnergyLevelsLength                 ( void );
    proshade_double                           getTraceSigmaVectorValue              ( proshade_unsign pos = 0 );
    proshade_unsign                           getTraceSigmaLength                   ( void );
    proshade_double                           getRotationFunctionVectorValue        ( proshade_unsign pos = 0 );
    proshade_unsign                           getRotationFunctionLength             ( void );
    
    //======================================== Symmetry accessor functions
    proshade_unsign                           getNoSymmetryAxes                     ( void );
    
public:
    //======================================== Distances results accessor functions
    std::vector< proshade_double >            getEnergyLevelsVector                 ( void );
    std::vector< proshade_double >            getTraceSigmaVector                   ( void );
    std::vector< proshade_double >            getRotationFunctionVector             ( void );

    //======================================== Symmetry results accessor functions
    std::string                               getSymmetryType                       ( void );
    proshade_unsign                           getSymmetryFold                       ( void );
    std::vector< std::string >                getSymmetryAxis                       ( proshade_unsign axisNo );
    
    //======================================== Re-boxing results accessor functions
    std::vector< proshade_signed >            getOriginalBounds                     ( proshade_unsign strNo );
    std::vector< proshade_signed >            getReBoxedBounds                      ( proshade_unsign strNo );
    proshade_double                           getMapValue                           ( proshade_unsign strNo, proshade_unsign mapIndex );
    
    //======================================== Overlay results accessor functions
    std::vector< proshade_double >            getEulerAngles                        ( void );
    std::vector< proshade_double >            getTranslation                        ( void );
};

//============================================ These functions should be in ProSHADE_run class, but I cannot make them work with Numpy from there, so they are here.
void getEnergyLevelsVectorNumpy               ( ProSHADE_run* run, int verbose, double *enLevVec, int len );
void getTraceSigmaVectorNumpy                 ( ProSHADE_run* run, int verbose, double *trSigVec, int len );
void getRotationFunctionVectorNumpy           ( ProSHADE_run* run, int verbose, double *rotFnVec, int len );
          
void getOriginalBoundsVectorNumpy             ( ProSHADE_run* run, proshade_unsign strNo, int *boundsVec,   int len );
void getReBoxedBoundsVectorNumpy              ( ProSHADE_run* run, proshade_unsign strNo, int *reboxVec,    int len );
void getReBoxedMap                            ( ProSHADE_run* run, proshade_unsign strNo, double *reboxMap, int len );
          
void getOptimalEulerAngles                    ( ProSHADE_run* run, double *eulerAngs, int len );
void getOptimalTranslation                    ( ProSHADE_run* run, double *translate, int len );

#endif
