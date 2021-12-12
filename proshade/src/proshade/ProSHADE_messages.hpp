/*! \file ProSHADE_messages.hpp
    \brief This header file contains all message function declarations.
 
    The functions declared in this headerr file (using the ProSHADE_internal_messages namespace) are used by ProSHADE to report various things to the user either using stdout
    (for progress messages) or the stderr (for error messages).
 
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

//============================================ ProSHADE
#include "ProSHADE_exceptions.hpp"

//============================================ Overinclusion protection
#ifndef PROSHADE_MESSAGES
#define PROSHADE_MESSAGES

//============================================ ProSHADE_internal_messages Namespace
/*! \namespace ProSHADE_internal_messages
    \brief This namespace contains all the functions used for reporting back to the user.
 
    The ProSHADE_internal_messages namespace wraps all functions for reporting back to the user. The user should not need to access this namespace when using
    the  library.
 */
namespace ProSHADE_internal_messages
{
    void printWellcomeMessage                         ( proshade_signed verbose );
    void printTerminateMessage                        ( proshade_signed verbose );
    void printProgressMessage                         ( proshade_signed verbose, proshade_signed messageLevel, std::string message, proshade_signed messageShift = 0 );
    void printWarningMessage                          ( proshade_signed verbose, std::string message, std::string warnCode );
    void printHelp  [[noreturn]]                      ( void );
}

#endif
