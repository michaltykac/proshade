/*! \file pyProSHADE_maths.cpp
    \brief This file contains the PyBind11 bindings for selected ProSHADE_internal_maths namespace functions.
    
    This file provides the bindings for some selected ProSHADE_internal_maths namespace members and functions. These are useful for de-bugging,
    but they may not be directly useful to the user.
    
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
void add_mathsNamespace ( pybind11::module& pyProSHADE )
{
    pyProSHADE.def                                    ( "getAngleAxisFromEulerAngles",
    [] ( proshade_double eulerAlpha, proshade_double eulerBeta, proshade_double eulerGamma ) -> pybind11::array_t < proshade_double >
    {
        //== Allocate output memory
        proshade_double* retVals = new proshade_double[4];
        ProSHADE_internal_misc::checkMemoryAllocation ( retVals, __FILE__, __LINE__, __func__ );
        
        //== Allocate workspace memory
        proshade_double* rMat = new proshade_double[9];
        ProSHADE_internal_misc::checkMemoryAllocation ( rMat, __FILE__, __LINE__, __func__ );
        proshade_double aX, aY, aZ, aA;

        //== Compute the conversion
        ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( eulerAlpha, eulerBeta, eulerGamma, rMat );
        ProSHADE_internal_maths::getAxisAngleFromRotationMatrix ( rMat, &aX, &aY, &aZ, &aA );
        
        //== Release workspace memory
        delete[] rMat;
        
        //== Save output
        retVals[0] = aX;
        retVals[1] = aY;
        retVals[2] = aZ;
        retVals[3] = aA;

        //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
        pybind11::capsule pyCapsuleAAFromEA ( retVals, []( void *f ) { proshade_double* foo = reinterpret_cast< proshade_double* > ( f ); delete foo; } );

        //== Copy the value
        pybind11::array_t < proshade_double > retArr = pybind11::array_t < proshade_double > ( { 4 },                          // Shape
                                                                                               { sizeof(proshade_double) },    // C-stype strides
                                                                                               retVals,                        // Data
                                                                                               pyCapsuleAAFromEA );            // Capsule

        //== Done
        return ( retArr );
    }, "This function converts the ZXZ Euler angles to angle-axis representation, returning numpy.ndarray with the following four numbers: [0] = x-axis; [1] = y-axis; [2] = z-axis; [3] = angle." );
}
