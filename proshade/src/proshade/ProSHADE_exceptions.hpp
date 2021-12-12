/*! \file ProSHADE_exceptions.hpp
    \brief This header  file contains the ProSHADE_exception class, which is used to report errors occuring during run-time.
 
    This file declares the ProSHADE_exception class. This class is derived from runtime_error and is therefore throwable; it extends the standard
    exceptions class by containing a bit more information in a specific format, so that the error messages resulting from throwing this class will be
    a bit more informative.
 
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
#include "ProSHADE_precomputedValues.hpp"

//==================================================== Overinclusion protection
#ifndef PROSHADE_EXCEPTIONS
#define PROSHADE_EXCEPTIONS

/*! \class ProSHADE_exception
    \brief This class is the representation of ProSHADE exception.
 
    An object of this class is thrown whenever ProSHADE encounters exceptional case and needs handling it in unusual
    manner. It is a slight expansion on the usual C++ extension class.
 */
class ProSHADE_exception : public std::runtime_error
{
    std::string errc;
    std::string file;
    std::string func;
    std::string info;
    int long line;
    
public:
    ProSHADE_exception ( const char* msg, std::string errc_, std::string file_, unsigned int line_, std::string func_, std::string info_): std::runtime_error ( msg )
    {
        this->errc                                    = errc_;
        this->file                                    = file_;
        this->line                                    = line_;
        this->func                                    = func_;
        this->info                                    = info_;
    }
    
    virtual std::string get_errc                      ( void );
    virtual std::string get_file                      ( void );
    virtual int long get_line                         ( void );
    virtual std::string get_func                      ( void );
    virtual std::string get_info                      ( void );
};

#endif
