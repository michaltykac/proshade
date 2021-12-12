/*! \file pyProSHADE_bindings.cpp
    \brief This file contains and combines the PyBind11 bindings for the ProSHADE python module.
    
    This file defines the python module for PyBind11, sets the main values and then includes all of the bindings from the other files.
    
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

//==================================================== Include getopt_port for python
#include <getopt_port/getopt_port.h>
#include <getopt_port/getopt_port.c>

//==================================================== Include ProSHADE
#include "ProSHADE.hpp"

//==================================================== Include full ProSHADE - including cpp files looks horrible, but it is the only way I can find to stop PyBind11 from complaining on Windows10. I believe this has to do with windows not exporting symbols unless the __cdecl_dllexport is used, which for most functions is not yet there...
#include "ProSHADE_precomputedValues.cpp"
#include "ProSHADE_exceptions.cpp"
#include "ProSHADE_misc.cpp"
#include "ProSHADE_maths.cpp"
#include "ProSHADE_tasks.cpp"
#include "ProSHADE_io.cpp"
#include "ProSHADE_data.cpp"
#include "ProSHADE_symmetry.cpp"
#include "ProSHADE_overlay.cpp"
#include "ProSHADE_wignerMatrices.cpp"
#include "ProSHADE_spheres.cpp"
#include "ProSHADE_mapManip.cpp"
#include "ProSHADE_messages.cpp"
#include "ProSHADE_distances.cpp"
#include "ProSHADE_peakSearch.cpp"
#include "ProSHADE_sphericalHarmonics.cpp"
#include "ProSHADE.cpp"

//==================================================== Include PyBind11 header
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

//==================================================== Forward declarations of pyProSHADE functions
void    add_settingsClass                             ( pybind11::module& pyProSHADE );
void    add_dataClass                                 ( pybind11::module& pyProSHADE );
void    add_distancesClass                            ( pybind11::module& pyProSHADE );

//==================================================== Remove the bindings that are not modifyable in python
PYBIND11_MAKE_OPAQUE                                  ( std::vector < std::string > )

//==================================================== Enable MSVC C4996 Warnings for the rest of the code
#if defined ( _MSC_VER )
    #pragma warning ( default:4996 )
#endif

//==================================================== Now the flags can be restored and used as per the CMakeLists.txt file.
#if defined ( __GNUC__ )
    #pragma GCC diagnostic pop
#endif

//==================================================== Now the flags can be restored and used as per the CMakeLists.txt file.
#if defined ( __clang__ )
    #pragma clang diagnostic pop
#endif

//==================================================== Include the other codes
#include "pyProSHADE.cpp"
#include "pyProSHADE_data.cpp"
#include "pyProSHADE_distances.cpp"
#include "pyProSHADE_maths.cpp"

//==================================================== Declare the exported functions
PYBIND11_MODULE ( proshade, pyProSHADE )
{
    //================================================ Create new, modifyable bindings
    pybind11::bind_vector < std::vector < std::string > > ( pyProSHADE, "<VectorOfStrings class> (Use append to add entries and [] to access them)", pybind11::module_local ( true ) );
    
    //================================================ Set the module description
    pyProSHADE.doc ( )                                = "Protein Shape Description and Symmetry Detection (ProSHADE) python module"; // Module docstring
    
    //================================================ Set the module version
    pyProSHADE.attr ( "__version__" )                 = PROSHADE_VERSION;

    //================================================ Export the ProSHADE_Task enum
    pybind11::enum_ < ProSHADE_Task > ( pyProSHADE, "ProSHADE_Task" )
        .value                                        ( "NA",          NA          )
        .value                                        ( "Distances",   Distances   )
        .value                                        ( "Symmetry",    Symmetry    )
        .value                                        ( "OverlayMap",  OverlayMap  )
        .value                                        ( "MapManip",    MapManip    )
        .export_values                                ( );
    
    //================================================ Export the ProSHADE_Settings class
    add_settingsClass                                 ( pyProSHADE );
    add_dataClass                                     ( pyProSHADE );
    add_distancesClass                                ( pyProSHADE );
    add_mathsNamespace                                ( pyProSHADE );
}
