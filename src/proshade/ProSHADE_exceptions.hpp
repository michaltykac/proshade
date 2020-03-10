/*! \file ProSHADE_exceptions.hpp
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
#include "ProSHADE_precomputedValues.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE_EXCEPTIONS__
#define __PROSHADE_EXCEPTIONS__

/*! \class ProSHADE_exception
 \brief This class is the representation of ProSHADE exception.
 
 An object of this class is thrown whenever ProSHADE encounters exceptional case and needs handling it in unusual
 manner. It is a slight expansion on the usual C++ extension class.
 */
class ProSHADE_exception : public std::runtime_error
{
    std::string errc;
    std::string file;
    unsigned int line;
    std::string func;
    std::string info;
    
public:
    ProSHADE_exception ( const char* msg, std::string errc_, std::string file_, unsigned int line_, std::string func_, std::string info_): std::runtime_error ( msg )
    {
        this->errc                            = errc_;
        this->file                            = file_;
        this->line                            = line_;
        this->func                            = func_;
        this->info                            = info_;
    }
   ~ProSHADE_exception ( ) throw ( ) { }
    
    std::string get_errc () const { return errc; }
    std::string get_file () const { return file; }
    unsigned int get_line() const { return line; }
    std::string get_func () const { return func; }
    std::string get_info () const { return info; }
};

#endif
