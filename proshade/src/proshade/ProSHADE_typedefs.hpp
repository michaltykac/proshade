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
    \version   0.7.5.4
    \date      MAR 2021
 */

//==================================================== ProSHADE
#include "ProSHADE_version.hpp"

//==================================================== MSVC Specific definition to allow M_PI 
#define _USE_MATH_DEFINES

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

//==================================================== Overinclusion protection
#ifndef __PROSHADE_TYPEDEFS__
#define __PROSHADE_TYPEDEFS__

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

