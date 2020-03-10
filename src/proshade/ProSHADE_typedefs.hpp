/*! \file ProSHADE_typedefs.hpp
  \brief This header contains the type definitions for ProSHADE.
 
 The types deined here are used throughout ProSHADE. The use of this is that should
 different type prove to be more appropriate or have different size on different
 platforms, all types can be changed in a single file.
 
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
#include "ProSHADE_version.hpp"

//============================================ Standard library
#include <iostream>
#include <ctime>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>
#include <exception>
#include <complex>
#include <algorithm>
#include <tgmath.h>
#include <getopt.h>

//============================================ Overinclusion protection
#ifndef __PROSHADE_TYPEDEFS__
#define __PROSHADE_TYPEDEFS__

//======================================== The Task data type
enum ProSHADE_Task { NA, Distances, Symmetry, OverlayMap, MapManip };

//============================================ ProSHADE Typedefs
typedef float                                 proshade_single;
typedef double                                proshade_double;
typedef signed long int                       proshade_signed;
typedef unsigned long int                     proshade_unsign;
typedef double                                proshade_complex[2];
typedef double                                proshade_triplet[3];

#endif

