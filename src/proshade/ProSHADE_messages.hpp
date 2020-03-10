/*! \file ProSHADE_messages.hpp
 \brief This header contains all message functions.
 
 The functions declared in here are used by ProSHADE to report various things to the user either using stdout (for progress messages) or the stderr
 (for error messages).
 
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
#include "ProSHADE_exceptions.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE_MESSAGES__
#define __PROSHADE_MESSAGES__

//============================================ ProSHADE_internal_messages Namespace
/*! \namespace ProSHADE_internal_messages
 \brief This namespace contains all the functions used for reporting back to the user.
 
 The ProSHADE_internal_messages namespace wraps all functions for reporting back to the user. The user should not need to access this namespace when using
 the  library.
 */
namespace ProSHADE_internal_messages
{
    void printWellcomeMessage                 ( proshade_signed verbose );
    void printTerminateMessage                ( proshade_signed verbose );
    void printProgressMessage                 ( proshade_signed verbose, proshade_signed messageLevel, std::string message );
    void printWarningMessage                  ( proshade_signed verbose, std::string message, std::string warnCode );
    void printHelp                            ( void );
}

#endif
