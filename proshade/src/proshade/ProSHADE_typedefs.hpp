/*! \file ProSHADE_typedefs.hpp
    \brief This header contains the type definitions for ProSHADE.
 
    The types deined in this header file are used throughout ProSHADE instead of C++ native typees. The use of this is that should
    different type prove to be more appropriate or have different size on different platforms, all types can be changed in a single file.
 
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
#include "ProSHADE_version.hpp"

//==================================================== Standard library
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>
#include <exception>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <utility>

//==================================================== Do not use the following flags for the included files - this causes a lot of warnings that have nothing to do with ProSHADE
#if defined ( __GNUC__ )
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wpedantic"
    #pragma GCC diagnostic ignored "-Wshadow"
    #pragma GCC diagnostic ignored "-Wall"
    #pragma GCC diagnostic ignored "-Wextra"
    #pragma GCC diagnostic ignored "-Wdouble-promotion"
    #pragma GCC diagnostic ignored "-Wconversion"
#endif

//==================================================== Do not use the following flags for the included files - this causes a lot of warnings that have nothing to do with ProSHADE
#if defined ( __clang__ )
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wpedantic"
    #pragma clang diagnostic ignored "-Wshadow"
    #pragma clang diagnostic ignored "-Wall"
    #pragma clang diagnostic ignored "-Wextra"
    #pragma clang diagnostic ignored "-Wdouble-promotion"
    #pragma clang diagnostic ignored "-Weverything"
#endif

//==================================================== Remove MSVC C4996 Warnings caused by Gemmi code
#if defined ( _MSC_VER )
    #pragma warning ( disable:4996 )
#endif

//==================================================== Gemmi
#ifndef __PROSHADE_GEMMI_INCLUDE__
    #define __PROSHADE_GEMMI_INCLUDE__
    #include <gemmi/mmread.hpp>
    #include <gemmi/ccp4.hpp>
    #include <gemmi/it92.hpp>
    #include <gemmi/dencalc.hpp>
    #include <gemmi/fprime.hpp>
    #include <gemmi/gz.hpp>
#endif

//==================================================== Enable MSVC C4996 Warnings for the rest of the code
#if defined ( _MSC_VER )
    #pragma warning ( default:4996 )
#endif

//==================================================== FFTW3
#ifdef __cplusplus
extern "C" {
#endif
    
#include <fftw3.h>
    
#ifdef __cplusplus
}
#endif

//==================================================== SOFT
#ifdef __cplusplus
extern "C" {
#endif
    
#include <wrap_fftw.h>
#include <makeweights.h>
#include <s2_primitive.h>
#include <s2_cospmls.h>
#include <s2_legendreTransforms.h>
#include <s2_semi_fly.h>
#include <rotate_so3_utils.h>
#include <utils_so3.h>
#include <soft_fftw.h>
#include <rotate_so3_fftw.h>
    
#ifdef __cplusplus
}
#endif

//==================================================== Now the flags can be restored and used as per the CMakeLists.txt file.
#if defined ( __GNUC__ )
    #pragma GCC diagnostic pop
#endif

//==================================================== Now the flags can be restored and used as per the CMakeLists.txt file.
#if defined ( __clang__ )
    #pragma clang diagnostic pop
#endif

//==================================================== GetOpt port (BSD License, works on Windows as well as linux)
#include <getopt_port/getopt_port.h>


//==================================================== Define maths constants in case cmath fails to have them
#ifndef M_PI
    #define M_PI           3.14159265358979323846
#endif

#ifndef M_SQRT1_2
    #define M_SQRT1_2      0.707106781186547524401
#endif

//==================================================== Overinclusion protection
#ifndef PROSHADE_TYPEDEFS
#define PROSHADE_TYPEDEFS

//==================================================== The Task data type
enum ProSHADE_Task { NA, Distances, Symmetry, OverlayMap, MapManip };

//==================================================== ProSHADE Typedefs
typedef float                                         proshade_single;
typedef double                                        proshade_double;
typedef signed long int                               proshade_signed;
typedef unsigned long int                             proshade_unsign;
typedef double                                        proshade_complex[2];
typedef double                                        proshade_triplet[3];

#endif

