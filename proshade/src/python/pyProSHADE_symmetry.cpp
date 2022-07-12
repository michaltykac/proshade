/*! \file pyProSHADE_symmetry.cpp
    \brief This file contains the PyBind11 bindings for selected ProSHADE_internal_symmetry namespace functions.
    
    This file provides the bindings for some selected ProSHADE_internal_symmetry namespace members and functions. These are typically functions
    that process one ProSHADE_data object and find/process symmetry detection.
    
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
void add_symmetryNamespace ( pybind11::module& pyProSHADE )
{
//    pyProSHADE.def                                    ( "findReliableUnphasedSymmetries",
//    [] ( ProSHADE_settings* settings, proshade_signed verbose, proshade_signed messageShift, proshade_double tolerance ) -> pybind11::array_t < proshade_unsign >
//    {
//        //== Run the detection
//        std::vector< proshade_unsign > rels           = ProSHADE_internal_symmetry::findReliableUnphasedSymmetries ( &settings->allDetectedCAxes, verbose, messageShift, tolerance );
//        
//        //== Allocate memory for the numpy values
//        proshade_unsign* npVals                       = new proshade_unsign[static_cast<unsigned int> ( rels.size() )];
//        ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );
//
//        //== Copy values
//        for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( rels.size() ); iter++ ) { npVals[iter] = static_cast< proshade_unsign > ( rels.at(iter) ); }
//
//        //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
//        pybind11::capsule pyCapsuleRelSyms ( npVals, []( void *f ) { proshade_unsign* foo = reinterpret_cast< proshade_unsign* > ( f ); delete foo; } );
//
//        //== Copy the value
//        pybind11::array_t < proshade_unsign > retArr = pybind11::array_t< proshade_unsign > ( { static_cast<int> ( rels.size() ) },                          // Shape
//                                                                                              { sizeof(proshade_unsign) },                                   // C-stype strides
//                                                                                              npVals,                                                        // Data
//                                                                                              pyCapsuleRelSyms );                                            // Capsule
//
//        //== Done
//        return ( retArr );
//    }, "This function checks the list of detected axes (presumably from phaseless symmetry detection) and returns the best dihedral (or cyclic, if no dihedral is found) point group, or empty vector if nothing is found.", pybind11::arg ( "settings" ), pybind11::arg ( "verbose" ), pybind11::arg ( "messageShift" ) = 1, pybind11::arg ( "tolerance" ) = 0.1 );
    
    pyProSHADE.def                                    ( "optimiseDGroupAngleFromAxesHeights",
    [] ( pybind11::array_t < proshade_unsign > selection, ProSHADE_internal_data::ProSHADE_data* dataObj, ProSHADE_settings* settings )
    {
        //== Sanity check
        pybind11::buffer_info selection_buf = selection.request();
        if ( selection_buf.ndim != 1 ) { std::cerr << "!!! ProSHADE PYTHON MODULE ERROR !!! The first argument to optimiseDGroupAngleFromAxesHeights() must be a 1D numpy array!" << std::endl; exit ( EXIT_FAILURE ); }
        if ( selection_buf.shape.at(0) != 2 ) { std::cerr << "!!! ProSHADE PYTHON MODULE ERROR !!! The first argument to optimiseDGroupAngleFromAxesHeights() must be array of length 2!" << std::endl; exit ( EXIT_FAILURE ); }
        
        //== Copy data to C++ format
        proshade_unsign* arrStart = static_cast< proshade_unsign* > ( selection_buf.ptr );
        std::vector< proshade_unsign > sel;
        ProSHADE_internal_misc::addToUnsignVector     ( &sel, arrStart[0] );
        ProSHADE_internal_misc::addToUnsignVector     ( &sel, arrStart[1] );
        
        //== Convert data type
        std::vector< std::vector< proshade_double > > allCs;
        std::vector< proshade_double > hlpVec;
        for ( size_t it1 = 0; it1 < dataObj->getCyclicAxes()->size(); it1++ )
        {
            hlpVec.clear                                  ( );
            for ( size_t it2 = 0; it2 < 7; it2++ ) { ProSHADE_internal_misc::addToDoubleVector ( &hlpVec, dataObj->getCyclicAxes()->at(it1)[it2] ); }
            ProSHADE_internal_misc::addToDoubleVectorVector ( &allCs, hlpVec );
        }
        
        //== Call the C++ function
        ProSHADE_internal_symmetry::optimiseDGroupAngleFromAxesHeights ( &allCs, sel, dataObj, settings );
        
        //== Done (all changes are in the data object, so nothing needs to be passed to Python)
        return ;
    }, "This function takes two axes with almost dihedral angle and optimises their relative positions as well as orientation with respect to the optimal angle and the rotation function.", pybind11::arg ( "selection" ), pybind11::arg ( "dataObj" ), pybind11::arg ( "settings" ) );
    
    pyProSHADE.def                                    ( "findPointFromTranslations",
    [] ( ProSHADE_settings* settings, ProSHADE_internal_data::ProSHADE_data* dataObj, pybind11::array_t < float > allCs, proshade_unsign whichAxis ) -> pybind11::array_t < proshade_double >
    {
        //== Sanity check
        pybind11::buffer_info allCs_buf = allCs.request();
        if ( allCs_buf.ndim != 2 ) { std::cerr << "!!! ProSHADE PYTHON MODULE ERROR !!! The third argument to findPointFromTranslations() must be a 2D numpy array!" << std::endl; exit ( EXIT_FAILURE ); }
        
        //== Copy to C++ data format
        std::vector< std::vector < proshade_double > > allCAxes;
        std::vector < proshade_double > hlpVec;
        float* arrStart = static_cast< float* > ( allCs_buf.ptr );
        for ( proshade_unsign axIt = 0; axIt < static_cast< proshade_unsign > ( allCs_buf.shape.at(0) ); axIt++ )
        {
            //== Clean up
            hlpVec.clear                              ( );
            
            //== Copy vals
            for ( proshade_unsign vIt = 0; vIt < 7; vIt++ )
            {
                ProSHADE_internal_misc::addToDoubleVector ( &hlpVec, static_cast< proshade_double > ( arrStart[(axIt*7)+vIt] ) );
            }
            
            ProSHADE_internal_misc::addToDoubleVectorVector ( &allCAxes, hlpVec );
        }
        
        //== Prepare all required memory
        fftw_complex *origMap = nullptr, *origCoeffs = nullptr, *rotMapComplex = nullptr, *rotCoeffs = nullptr, *trFunc = nullptr, *trFuncCoeffs = nullptr;
        fftw_plan planForwardFourier, planForwardFourierRot, planReverseFourierComb;
        ProSHADE_internal_symmetry::allocateCentreOfMapFourierTransforms ( dataObj->getXDim(), dataObj->getYDim(), dataObj->getZDim(), origMap, origCoeffs, rotMapComplex, rotCoeffs, trFunc, trFuncCoeffs, &planForwardFourier, &planForwardFourierRot, &planReverseFourierComb );

        //== Compute Fourier for the original map
        for ( size_t it = 0; it < static_cast< size_t > ( dataObj->getXDim() * dataObj->getYDim() * dataObj->getZDim() ); it++ ) { origMap[it][0] = dataObj->getMapValue( it ); origMap[it][1] = 0.0; }
        fftw_execute                                  ( planForwardFourier );
        
        //== Run C++ code
        std::vector< proshade_unsign > axLst;
        std::vector< std::vector< proshade_double > > symElems;
        ProSHADE_internal_misc::addToUnsignVector     ( &axLst, whichAxis );
        symElems                                      = dataObj->getAllGroupElements ( &allCAxes, axLst, "C", settings->axisErrTolerance );
        std::vector< proshade_double > pointPos       = ProSHADE_internal_symmetry::findPointFromTranslations ( dataObj,
                                                                                                                symElems,
                                                                                                                origCoeffs, rotMapComplex,
                                                                                                                rotCoeffs, planForwardFourierRot,
                                                                                                                trFuncCoeffs, trFunc,
                                                                                                                planReverseFourierComb );
        
        //== Release required memory
        ProSHADE_internal_symmetry::releaseCentreOfMapFourierTransforms ( origMap, origCoeffs, rotMapComplex, rotCoeffs, trFunc, trFuncCoeffs, planForwardFourier, planForwardFourierRot, planReverseFourierComb );
        
        //== Allocate memory for the numpy values
        proshade_double* npVals                       = new proshade_double[3];
        ProSHADE_internal_misc::checkMemoryAllocation ( npVals, __FILE__, __LINE__, __func__ );

        //== Copy values
        for ( proshade_unsign iter = 0; iter < 3; iter++ ) { npVals[iter] = static_cast< proshade_double > ( pointPos.at(iter) ); }

        //== Create capsules to make sure memory is released properly from the allocating language (C++ in this case)
        pybind11::capsule pyCapsuleSymPoint ( npVals, []( void *f ) { proshade_double* foo = reinterpret_cast< proshade_double* > ( f ); delete foo; } );

        //== Copy the value
        pybind11::array_t < proshade_double > retArr = pybind11::array_t< proshade_double > ( { static_cast<int> ( 3 ) },                                    // Shape
                                                                                              { sizeof(proshade_double) },                                   // C-stype strides
                                                                                              npVals,                                                        // Data
                                                                                              pyCapsuleSymPoint );                                           // Capsule

        //== Done
        return ( retArr );
        
        
    }, "This function computes the average of optimal translations for a cyclic point group.", pybind11::arg ( "settings" ), pybind11::arg ( "dataObj" ), pybind11::arg ( "allCs" ), pybind11::arg ( "whichAxis" ) );
}
