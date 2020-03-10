/*! \file ProSHADE_tasks.hpp
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
#include "ProSHADE_overlay.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE_TASKS__
#define __PROSHADE_TASKS__

//============================================ ProSHADE_internal_tasks Namespace
/*! \namespace ProSHADE_internal_tasks
 \brief This namespace contains the main driving functions for each task.
 
 The ProSHADE_internal_tasks namespace contains the driving functions for all the tasks that can be
 accomplished by the ProSHADE tool. The user should not need to access
 this namespace when using the library.
 */
namespace ProSHADE_internal_tasks
{
    void MapManipulationTask                  ( ProSHADE_settings* settings, std::vector < proshade_signed* >* originalBounds,
                                                std::vector < proshade_signed* >* reboxedBounds, std::vector < proshade_double* >* manipulatedMaps );
    void DistancesComputationTask             ( ProSHADE_settings* settings, std::vector< proshade_double >* enLevs, std::vector< proshade_double >* trSigm,
                                                std::vector< proshade_double >* rotFun );
    void SymmetryDetectionTask                ( ProSHADE_settings* settings, std::vector< proshade_double* >* axes );
    void MapOverlayTask                       ( ProSHADE_settings* settings, std::vector < proshade_double >* eulerAngles,
                                                std::vector < proshade_double >* translation );
    
    void ReportDistancesResults               ( ProSHADE_settings* settings, std::string str1, std::string str2, proshade_double enLevDist,
                                                proshade_double trSigmDist, proshade_double rotFunDist );
    
    void checkMapManipulationSettings         ( ProSHADE_settings* settings );
    void checkDistancesSettings               ( ProSHADE_settings* settings );
    void checkSymmetrySettings                ( ProSHADE_settings* settings );
    void checkOverlaySettings                 ( ProSHADE_settings* settings );
}



#endif
