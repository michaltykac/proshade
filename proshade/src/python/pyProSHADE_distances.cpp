/*! \file pyProSHADE_distances.cpp
    \brief This file contains the PyBind11 bindings for the ProSHADE_internal_distances namespace.
    
    This file provides the bindings for the ProSHADE_internal_distances namespace functions.
    
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

//==================================================== Include PyBind11 header
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

//==================================================== Add the ProSHADE_settings and ProSHADE_run classes to the PyBind11 module
void add_distancesClass ( pybind11::module& pyProSHADE )
{
    pyProSHADE.def                                    ( "computeEnergyLevelsDescriptor",     &ProSHADE_internal_distances::computeEnergyLevelsDescriptor,     "This function computes the energy levels descriptor value between two objects.",         pybind11::arg ( "obj1" ), pybind11::arg ( "obj2" ), pybind11::arg ( "settings" ) );
    pyProSHADE.def                                    ( "computeTraceSigmaDescriptor",       &ProSHADE_internal_distances::computeTraceSigmaDescriptor,       "This function computes the trace sigma descriptor value between two objects.",         pybind11::arg ( "obj1" ), pybind11::arg ( "obj2" ), pybind11::arg ( "settings" ) );
    pyProSHADE.def                                    ( "computeRotationFunctionDescriptor", &ProSHADE_internal_distances::computeRotationFunctionDescriptor, "This function computes the rotation function descriptor value between two objects.",   pybind11::arg ( "obj1" ), pybind11::arg ( "obj2" ), pybind11::arg ( "settings" ) );
}
