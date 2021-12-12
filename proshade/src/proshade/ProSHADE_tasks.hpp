/*! \file ProSHADE_tasks.hpp
    \brief This header declares all the taks functions.
 
    The ProSHADE_internal_tasks namespace declared in this header file contains all the task functions, which drive the progression towards computation of a particular task results. It also contains the functions
    required to check that all the input information required was supplied.
 
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

//==================================================== Overinclusion protection
#ifndef PROSHADE_TASKS
#define PROSHADE_TASKS

//==================================================== ProSHADE_internal_tasks Namespace
/*! \namespace ProSHADE_internal_tasks
    \brief This namespace contains the main driving functions for each task.
 
    The ProSHADE_internal_tasks namespace contains the driving functions for all the tasks that can be
    accomplished by the ProSHADE tool. The user should not need to access
    this namespace when using the library.
 */
namespace ProSHADE_internal_tasks
{
    void MapManipulationTask                          ( ProSHADE_settings* settings, std::vector < proshade_signed* >* originalBounds,
                                                        std::vector < proshade_signed* >* reboxedBounds, std::vector < proshade_double* >* manipulatedMaps );
    void DistancesComputationTask                     ( ProSHADE_settings* settings, std::vector< proshade_double >* enLevs, std::vector< proshade_double >* trSigm,
                                                        std::vector< proshade_double >* rotFun );
    void SymmetryDetectionTask                        ( ProSHADE_settings* settings, std::vector< proshade_double* >* axes, std::vector < std::vector< proshade_double > >* allCs,
                                                        std::vector< proshade_double >* mapCOMShift );
    void MapOverlayTask                               ( ProSHADE_settings* settings, std::vector < proshade_double >* rotationCentre, std::vector < proshade_double >* eulerAngles,
                                                        std::vector < proshade_double >* finalTranslation );
    void SymmetryCentreDetectionTask                  ( ProSHADE_settings* settings, std::vector < std::vector< proshade_double > >* allCs, std::vector< proshade_double* >* axes,
                                                        proshade_unsign strIndex = 0 );

    void ReportDistancesResults                       ( ProSHADE_settings* settings, std::string str1, std::string str2, proshade_double enLevDist,
                                                        proshade_double trSigmDist, proshade_double rotFunDist );
            
    void checkMapManipulationSettings                 ( ProSHADE_settings* settings );
    void checkDistancesSettings                       ( ProSHADE_settings* settings );
    void checkSymmetrySettings                        ( ProSHADE_settings* settings );
    void checkOverlaySettings                         ( ProSHADE_settings* settings );
}



#endif
