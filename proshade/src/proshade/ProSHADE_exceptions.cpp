/*! \file ProSHADE_exceptions.cpp
    \brief This source  file contains the ProSHADE_exception class functions..
 
    This file simply declares the ProSHADE_exceptions class functions in order for them to be declared in a translation unit and
    therefore to avoid the weak vtable warning.
 
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
#include "ProSHADE_exceptions.hpp"

/*! \brief This function returns the exception error code.
 
    \param[out] errc Error code of the exception.
 */
std::string ProSHADE_exception::get_errc ( ) { return ( this->errc ); }

/*! \brief This function returns the exception location file name.
 
    \param[out] file The file from where the exception originated.
 */
std::string ProSHADE_exception::get_file ( ) { return ( this->file ); }

/*! \brief This function returns the exception  location line.
 
    \param[out] line The line of code that caused the exception.
 */
int long ProSHADE_exception::get_line ( ) { return ( this->line ); }

/*! \brief This function returns the exception causing function name.
 
    \param[out] func The name of the function from which the exception has originated.
 */
std::string ProSHADE_exception::get_func ( ) { return ( this->func ); }

/*! \brief This function returns the exception description.
 
    \param[out] info The information provided with the exception.
 */
std::string ProSHADE_exception::get_info ( ) { return ( this->info ); }
