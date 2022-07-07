/*! \file pyProSHADE_mapManip.cpp
    \brief This file contains the PyBind11 bindings for selected ProSHADE_internal_mapManip namespace functions.
    
    This file provides the bindings for some selected ProSHADE_internal_mapManip namespace members and functions. These are typically functions
    that process the internal map and/or manipulate it in some way.
    
    Copyright by Michal Tykac and individual contributors. All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
 
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.6.6
    \date      JUL 2022
 */

//==================================================== Include PyBind11 header
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

//==================================================== Add the ProSHADE_settings and ProSHADE_run classes to the PyBind11 module
void add_mapManipNamespace ( pybind11::module& pyProSHADE )
{
    pyProSHADE.def                                    ( "findMAPCOMValues",
    [] ( ProSHADE_internal_data::ProSHADE_data* dataObj, ProSHADE_settings* settings ) -> pybind11::array_t < proshade_double >
    {
        //== Run the detection
        proshade_double xMapCOM = 0.0, yMapCOM = 0.0, zMapCOM = 0.0;
        ProSHADE_internal_mapManip::findMAPCOMValues  ( dataObj->internalMap,
                                                       &xMapCOM, &yMapCOM, &zMapCOM,
                                                        dataObj->xDimSize, dataObj->yDimSize, dataObj->zDimSize,
                                                        dataObj->xFrom, dataObj->xTo,
                                                        dataObj->yFrom, dataObj->yTo,
                                                        dataObj->zFrom, dataObj->zTo, settings->removeNegativeDensity );
        
        //== Allocate memory for the numpy values
        proshade_double* npVals                       = new proshade_double[3];
        ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );

        //== Copy values
        npVals[0]                                     = xMapCOM;
        npVals[1]                                     = yMapCOM;
        npVals[2]                                     = zMapCOM;

        //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
        pybind11::capsule pyCapsuleMapCOM ( npVals, []( void *f ) { proshade_double* foo = reinterpret_cast< proshade_double* > ( f ); delete foo; } );

        //== Copy the value
        pybind11::array_t < proshade_double > retArr = pybind11::array_t< proshade_double > ( { static_cast<int> ( 3 ) },                                    // Shape
                                                                                              { sizeof(proshade_double) },                                   // C-stype strides
                                                                                              npVals,                                                        // Data
                                                                                              pyCapsuleMapCOM );                                             // Capsule

        //== Done
        return                                        ( retArr );
    }, "This function takes the proshade map object and procceds to compute the Centre of Mass of this object in Angstroms in real space.", pybind11::arg ( "dataObj" ), pybind11::arg ( "settings" ) );
    
}
