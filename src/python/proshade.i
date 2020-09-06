/*! \file proshade.i
    \brief This is the SWIG interface file for creating ProSHADE python bindings.

    This file contains the SWIG C++ to python interface, which allows SWIG to take the C++ code of
    ProSHADE and produce a working Python module. This file should work automatically with the CMake
    ProSHADE installation and the user should not need to do anything with it.

    Copyright by Michal Tykac and individual contributors. All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no
    event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or
    profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility 
    of such damage.
 
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.4
    \date      SEP 2020
 */

//============================================ Module name
%module proshade

//============================================ C++ includes and definitions
%{
//============================================ This definition is required! Do not forget this
#define SWIG_FILE_WITH_INIT
    
//============================================ Include the standard modules
#include <vector>
#include <string>
#include <memory>

//============================================ Include Python and ProSHADE files which are to be pythonised
#include "Python.h"
#include "ProSHADE_typedefs.hpp"
#include "ProSHADE_settings.hpp"
#include "ProSHADE_data.hpp"
#include "ProSHADE_distances.hpp"
#include "ProSHADE.hpp"
%}

//============================================ Include into SWIG the typemap files
%include "std_string.i"
%include "std_vector.i"
%include "std_shared_ptr.i"
%include "stl.i"
%include "typemaps.i"
%include "numpy.i"

//============================================ Initialise array work for numpy type conversions
%init %{
import_array();
%}

//============================================ Apply the numpy typemaps for distances
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *enLevVec, int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *trSigVec, int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *rotFnVec, int len ) }

//============================================ Apply the numpy typemaps for reboxing
%apply ( int*    ARGOUT_ARRAY1, int DIM1 ) { ( int *boundsVec,   int len ) }
%apply ( int*    ARGOUT_ARRAY1, int DIM1 ) { ( int *reboxVec,    int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *reboxMap, int len ) }

//============================================ Apply the numpy typemaps for overlay
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *eulerAngs, int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *translate, int len ) }

//============================================ Apply the numpy typemaps for ProSHADE_data functions
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *mapArrayPython,     int len ) }
%apply ( double* IN_ARRAY1,     int DIM1 ) { ( double *mapChangedInPython, int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *sphericalHarmsReal, int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *sphericalHarmsImag, int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *eMatsLMReal,        int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *eMatsLMImag,        int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *so3CoefsReal,       int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *so3CoefsImag,       int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *rotFunReal,         int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *rotFunImag,         int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *rotMat,             int len ) }
%apply ( double* IN_ARRAY1,     int DIM1 ) { ( double *mapVals,            int len ) }
%apply ( int*    ARGOUT_ARRAY1, int DIM1 ) { ( int* reBoxBounds,           int len ) }
%apply ( int*    IN_ARRAY1,     int DIM1 ) { ( int* newBounds,             int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *trsFunReal,         int len ) }
%apply ( double* ARGOUT_ARRAY1, int DIM1 ) { ( double *trsFunImag,         int len ) }

//============================================ Include the pythonised ProSHADE code to SWIG
%include "ProSHADE_typedefs.hpp"
%include "ProSHADE_settings.hpp"
%include "ProSHADE_data.hpp"
%include "ProSHADE_distances.hpp"
%include "ProSHADE.hpp"

//============================================ Create templates for non-numpy related type conversions
%template ( _string_list ) std::vector< std::string > ;
%template ( _float_list )  std::vector< proshade_single > ;
%template ( _double_list ) std::vector< proshade_double > ;

//============================================ Define some python functions for easier usage
//============================================ Distances accessors
%pythoncode %{

def getEnergyLevelsDescrNumpy ( pRun ):
    return ( getEnergyLevelsVectorNumpy       ( pRun, pRun.getVerbose(), pRun.getNoStructures( ) - 1 ) )
    
def getTraceSigmaDescrNumpy ( pRun ):
    return ( getTraceSigmaVectorNumpy         ( pRun, pRun.getVerbose(), pRun.getNoStructures( ) - 1 ) )
    
def getRotationFunctionDescrNumpy ( pRun ):
    return ( getRotationFunctionVectorNumpy   ( pRun, pRun.getVerbose(), pRun.getNoStructures( ) - 1 ) )
%}

//============================================ Symmetry
%pythoncode %{
def getDetectedSymmetryType ( pRun ):
    return                                    ( pRun.getSymmetryType ( ) )
    
def getDetectedSymmetryFold ( pRun ):
    return                                    ( pRun.getSymmetryFold ( ) )
    
def getDetectedSymmetryAxes ( pRun ):
    retArr                                    = []
    for iter in range( 0, pRun.getNoSymmetryAxes ( ) ):
        hlpArr                                = pRun.getSymmetryAxis ( iter )
        hlpTlp                                = ( hlpArr[0], float ( hlpArr[1] ), float ( hlpArr[2] ), float ( hlpArr[3] ), float ( hlpArr[4] ), float ( hlpArr[5] ) )
        retArr.append                         ( hlpTlp )
    return                                    ( retArr )
%}
    
//============================================ Reboxing
%pythoncode %{
def getOrigBounds ( pRun ):
    import numpy
    ret                                       = numpy.empty ( ( 0, 6 ) )
    for iter in range ( 0, pRun.getNoStructures ( ) ):
        ret                                   = numpy.append ( ret, [getOriginalBoundsVectorNumpy ( pRun, iter, 6 )], axis = 0 )
    return                                    ( ret )
        
def getReboxBounds ( pRun ):
    import numpy
    ret                                       = numpy.empty ( ( 0, 6 ) )
    for iter in range ( 0, pRun.getNoStructures ( ) ):
        ret                                   = numpy.append ( ret, [getReBoxedBoundsVectorNumpy ( pRun, iter, 6 )], axis = 0 )
    return                                    ( ret )
    
def getReboxMap ( pRun ):
    ret                                       = []
    reboxedBounds                             = getReboxBounds ( pRun )
    for iter in range ( 0, pRun.getNoStructures ( ) ):
        ret.append                            ( getReBoxedMap ( pRun, iter, int ( ( reboxedBounds[iter][1] - reboxedBounds[iter][0] + 1 ) *
                                                                                  ( reboxedBounds[iter][3] - reboxedBounds[iter][2] + 1 ) *
                                                                                  ( reboxedBounds[iter][5] - reboxedBounds[iter][4] + 1 ) ) ) )
    return                                    ( ret )
%}

//============================================ Overlay
%pythoncode %{
def getEulerAngles ( pRun ):
    return                                    ( getOptimalEulerAngles ( pRun, 3 ) )
    
def getTranslation ( pRun ):
    return                                    ( getOptimalTranslation ( pRun, 3 ) )
    
def getRotationMat ( pRun ):
    import numpy
    eAngs                                     = getEulerAngles ( pRun )
    
    ret                                       = numpy.empty ( ( 3, 3 ) )
    
    ret[0][0]                                 =  numpy.cos ( eAngs[0] ) * numpy.cos ( eAngs[1]  ) * numpy.cos ( eAngs[2] ) - numpy.sin ( eAngs[0] ) * numpy.sin ( eAngs[2] );
    ret[0][1]                                 =  numpy.sin ( eAngs[0] ) * numpy.cos ( eAngs[1]  ) * numpy.cos ( eAngs[2] ) + numpy.cos ( eAngs[0] ) * numpy.sin ( eAngs[2] );
    ret[0][2]                                 = -numpy.sin ( eAngs[1] ) * numpy.cos ( eAngs[2] );
    
    ret[1][0]                                 = -numpy.cos ( eAngs[0] ) * numpy.cos ( eAngs[1]  ) * numpy.sin ( eAngs[2] ) - numpy.sin ( eAngs[0] ) * numpy.cos ( eAngs[2] );
    ret[1][1]                                 = -numpy.sin ( eAngs[0] ) * numpy.cos ( eAngs[1]  ) * numpy.sin ( eAngs[2] ) + numpy.cos ( eAngs[0] ) * numpy.cos ( eAngs[2] );
    ret[1][2]                                 =  numpy.sin ( eAngs[1] ) * numpy.sin ( eAngs[2] );
    
    ret[2][0]                                 =  numpy.cos ( eAngs[0] ) * numpy.sin ( eAngs[1]  );
    ret[2][1]                                 =  numpy.sin ( eAngs[0] ) * numpy.sin ( eAngs[1]  );
    ret[2][2]                                 =  numpy.cos ( eAngs[1] );
    
    return                                    ( ret )
%}

//============================================ Spherical harmonics access
%pythoncode %{
def getSphericalHarmonics ( pStruct, verbose = -1 ):
    import numpy
    ret                                       = numpy.empty ( [ 0, pStruct.getSphericalHarmonicsLenForShell ( pStruct.noSpheres - 1, verbose ) ], dtype = complex )
    for sph in range ( 0, pStruct.noSpheres ):
        realSH                                = pStruct.getRealSphericalHarmonicsForShell ( sph, verbose, pStruct.getSphericalHarmonicsLenForShell ( sph, verbose ) )
        imagSH                                = pStruct.getImagSphericalHarmonicsForShell ( sph, verbose, pStruct.getSphericalHarmonicsLenForShell ( sph, verbose ) )

        if len( realSH ) < ret.shape[1]:
            for iter in range( len( realSH ), ret.shape[1] ):
                realSH                        = numpy.append ( realSH, [ 0.0 ] )
            
        if len( imagSH ) < ret.shape[1]:
            for iter in range( len( imagSH ), ret.shape[1] ):
                imagSH                        = numpy.append ( imagSH, [ 0.0 ] )
                
        ret                                   = numpy.append ( ret, [ realSH + 1j * imagSH ], axis = 0 )
                
    return                                    ( ret )
%}

//============================================ E Matrices access
%pythoncode %{
def getEMatrix ( pStruct ):
    import numpy
    ret                                       = numpy.empty ( [ pStruct.getMaxBand(),
                                                                len ( range( -pStruct.getMaxBand(), ( pStruct.getMaxBand() + 1 ) ) ),
                                                                len ( range( -pStruct.getMaxBand(), ( pStruct.getMaxBand() + 1 ) ) ) ], dtype = complex )
    for bnd in range ( 0, pStruct.getMaxBand() ):
        for ord1 in range ( -bnd, ( bnd + 1 ) ):
            realEVals                         = pStruct.getRealEMatrixValuesForLM ( bnd, ord1 + bnd, len( range( -bnd, ( bnd + 1 ) ) ) )
            imagEVals                         = pStruct.getImagEMatrixValuesForLM ( bnd, ord1 + bnd, len( range( -bnd, ( bnd + 1 ) ) ) )

            if len( realEVals ) < ret.shape[2]:
                for iter in range( len( realEVals ), ret.shape[2] ):
                    realEVals                 = numpy.append ( realEVals, [ 0.0 ] )
                
            if len( imagEVals ) < ret.shape[2]:
                for iter in range( len( imagEVals ), ret.shape[2] ):
                    imagEVals                 = numpy.append ( imagEVals, [ 0.0 ] )
                    
            ret[bnd][ord1]                    = realEVals + 1j * imagEVals
                
    return                                    ( ret )
%}
    
//============================================ SO(3) Coefficients access
%pythoncode %{
def getSO3Coeffs ( pStruct ):
    import numpy

    realSO3Coeffs                             = pStruct.getRealSO3Coeffs ( int ( ( 4 * numpy.power ( pStruct.getMaxBand(), 3 )  - pStruct.getMaxBand() ) / 3.0 ) )
    imagSO3Coeffs                             = pStruct.getImagSO3Coeffs ( int ( ( 4 * numpy.power ( pStruct.getMaxBand(), 3 )  - pStruct.getMaxBand() ) / 3.0 ) )
    ret                                       = realSO3Coeffs + 1j * imagSO3Coeffs
                
    return                                    ( ret )
%}

//============================================ Rotation function access
%pythoncode %{
def getRotationFunction1D ( pStruct ):
    import numpy

    realRotFun                                = pStruct.getRealRotFunction ( int ( numpy.power ( pStruct.getMaxBand() * 2.0, 3.0 ) ) )
    imagRotFun                                = pStruct.getImagRotFunction ( int ( numpy.power ( pStruct.getMaxBand() * 2.0, 3.0 ) ) )
    
    return                                    (  realRotFun + 1j * imagRotFun )
    
def getRotationFunction3D ( pStruct ):
    import numpy

    realRotFun                                = pStruct.getRealRotFunction ( int ( numpy.power ( pStruct.getMaxBand() * 2.0, 3.0 ) ) )
    imagRotFun                                = pStruct.getImagRotFunction ( int ( numpy.power ( pStruct.getMaxBand() * 2.0, 3.0 ) ) )
    ret                                       = numpy.empty ( [ int ( pStruct.getMaxBand() * 2.0 ), int ( pStruct.getMaxBand() * 2.0 ), int ( pStruct.getMaxBand() * 2.0 ) ], dtype = complex )
                
    for eA in range ( 0, int ( pStruct.getMaxBand() * 2.0 ) ):
       for eB in range ( 0, int ( pStruct.getMaxBand() * 2.0 ) ):
           for eG in range ( 0, int ( pStruct.getMaxBand() * 2.0 ) ):
               index                          = int ( eG + ( pStruct.getMaxBand() * 2.0 ) * ( eB + ( pStruct.getMaxBand() * 2.0 ) * eA ) )
               ret[eA][eB][eG]                = realRotFun[index] + 1j * imagRotFun[index]
    
    return                                    ( ret )
               
def getRotationMatrixFromRotFunIndices ( pStruct, first, second, third ):
    import numpy
               
    oneDMat                                   = pStruct.getRotMatrixFromRotFunInds ( first, second, third, 9  )
    ret                                       = numpy.empty ( [ 3, 3 ] )
    
    ret[0][0] = oneDMat[0]; ret[0][1] = oneDMat[1]; ret[0][2] = oneDMat[2];
    ret[1][0] = oneDMat[3]; ret[1][1] = oneDMat[4]; ret[1][2] = oneDMat[5];
    ret[2][0] = oneDMat[6]; ret[2][1] = oneDMat[7]; ret[2][2] = oneDMat[8];
    
    return                                    ( ret )
%}

//============================================ Symmetry axes access
%pythoncode %{
def getSymmetryAxesPython ( pStruct, pSet ):
    retArr                                    = []
    for iter in range( 0, pStruct.getNoSymmetryAxes ( pSet ) ):
        hlpArr                                = pStruct.getSymmetryAxis ( pSet, iter )
        hlpTlp                                = ( hlpArr[0], float ( hlpArr[1] ), float ( hlpArr[2] ), float ( hlpArr[3] ), float ( hlpArr[4] ), float ( hlpArr[5] ) )
        retArr.append                         ( hlpTlp )
    return                                    ( retArr )
%}

//============================================ Rotation matrix from Euler angles
%pythoncode %{
def getRotationMatrixFromEulerZXZ ( eAngs ):
    import numpy
    
    ret                                       = numpy.empty ( ( 3, 3 ) )
    
    ret[0][0]                                 =  numpy.cos ( eAngs[0] ) * numpy.cos ( eAngs[1]  ) * numpy.cos ( eAngs[2] ) - numpy.sin ( eAngs[0] ) * numpy.sin ( eAngs[2] );
    ret[0][1]                                 =  numpy.sin ( eAngs[0] ) * numpy.cos ( eAngs[1]  ) * numpy.cos ( eAngs[2] ) + numpy.cos ( eAngs[0] ) * numpy.sin ( eAngs[2] );
    ret[0][2]                                 = -numpy.sin ( eAngs[1] ) * numpy.cos ( eAngs[2] );
    
    ret[1][0]                                 = -numpy.cos ( eAngs[0] ) * numpy.cos ( eAngs[1]  ) * numpy.sin ( eAngs[2] ) - numpy.sin ( eAngs[0] ) * numpy.cos ( eAngs[2] );
    ret[1][1]                                 = -numpy.sin ( eAngs[0] ) * numpy.cos ( eAngs[1]  ) * numpy.sin ( eAngs[2] ) + numpy.cos ( eAngs[0] ) * numpy.cos ( eAngs[2] );
    ret[1][2]                                 =  numpy.sin ( eAngs[1] ) * numpy.sin ( eAngs[2] );
    
    ret[2][0]                                 =  numpy.cos ( eAngs[0] ) * numpy.sin ( eAngs[1]  );
    ret[2][1]                                 =  numpy.sin ( eAngs[0] ) * numpy.sin ( eAngs[1]  );
    ret[2][2]                                 =  numpy.cos ( eAngs[1] );
    
    return                                    ( ret )
%}

//============================================ Access internal map 3D
%pythoncode %{
def getMapPython3D ( pStruct ):
    import numpy
    
    oneDMap                                   = pStruct.getMapPython ( pStruct.getMapArraySizePython() )
    ret                                       = numpy.empty ( [ int ( pStruct.getXDim() ), int ( pStruct.getYDim() ), int ( pStruct.getZDim() ) ] )
    arrPos                                    = 0
    
    for xIt in range ( 0, int ( pStruct.getXDim() ) ):
       for yIt in range ( 0, int ( pStruct.getYDim() ) ):
           for zIt in range ( 0, int ( pStruct.getZDim() ) ):
               arrPos                         = zIt + pStruct.getZDim() * ( yIt + pStruct.getYDim() * xIt )
               ret[xIt][yIt][zIt]             = oneDMap[arrPos]
               
    return                                    ( ret )
               
def getMapPython1D ( pStruct ):
    import numpy
               
    return                                    ( pStruct.getMapPython ( pStruct.getMapArraySizePython() ) )
%}
    
//============================================ Set internal map 1D and 3D
%pythoncode %{
def setMapPython1D ( pStruct, map ):
    pStruct.setMapPython                      ( map )
    
def setMapPython3D ( pStruct, map ):
    import numpy
    
    map1D                                     = numpy.empty ( int ( pStruct.getXDim() ) * int ( pStruct.getYDim() ) * int ( pStruct.getZDim() ) )
    arrPos                                    = 0
    
    for xIt in range ( 0, int ( pStruct.getXDim() ) ):
       for yIt in range ( 0, int ( pStruct.getYDim() ) ):
           for zIt in range ( 0, int ( pStruct.getZDim() ) ):
               arrPos                         = zIt + pStruct.getZDim() * ( yIt + pStruct.getYDim() * xIt )
               map1D[arrPos]                  = map[xIt][yIt][zIt]
    
    pStruct.setMapPython                      ( map1D )
               
def setNewMapPython1D ( pStruct, map ):
    pStruct.setNewMapPython                   ( map )

def setNewMapPython3D ( pStruct, map ):
    import numpy
    
    map1D                                     = numpy.empty ( int ( pStruct.getXDim() ) * int ( pStruct.getYDim() ) * int ( pStruct.getZDim() ) )
    arrPos                                    = 0
    
    for xIt in range ( 0, int ( pStruct.getXDim() ) ):
        for yIt in range ( 0, int ( pStruct.getYDim() ) ):
            for zIt in range ( 0, int ( pStruct.getZDim() ) ):
                arrPos                        = zIt + pStruct.getZDim() * ( yIt + pStruct.getYDim() * xIt )
                map1D[arrPos]                 = map[xIt][yIt][zIt]
                
    pStruct.setNewMapPython                   ( map1D )
%}

//============================================ Set internal map 1D and 3D
%pythoncode %{
def convert3Dto1DArray ( array3D ):
    import numpy
    assert                                    ( array3D.ndim == 3)

    array1D                                   = numpy.empty ( array3D.shape[0] * array3D.shape[1] * array3D.shape[2] )
    arrPos                                    = 0

    for xIt in range( 0, array3D.shape[0] ):
        for yIt in range( 0, array3D.shape[1] ):
            for zIt in range( 0, array3D.shape[2] ):
                arrPos                        = zIt + array3D.shape[2] * ( yIt + array3D.shape[1] * xIt )
                array1D[arrPos]               = array3D[xIt][yIt][zIt]  

    return                                    ( array1D )              
%}

//============================================ Translation function access
%pythoncode %{
def getTranslationFunction1D ( pStruct ):
    import numpy

    realTrsFun                                = pStruct.getRealTranslationFunction ( pStruct.getXDim() * pStruct.getYDim() * pStruct.getZDim() )
    imagTrsFun                                = pStruct.getImagTranslationFunction ( pStruct.getXDim() * pStruct.getYDim() * pStruct.getZDim() )
    
    return                                    (  realTrsFun + 1j * imagTrsFun )
    
def getTranslationFunction3D ( pStruct ):
    import numpy

    realTrsFun                                = pStruct.getRealTranslationFunction ( pStruct.getXDim() * pStruct.getYDim() * pStruct.getZDim() )
    imagTrsFun                                = pStruct.getImagTranslationFunction ( pStruct.getXDim() * pStruct.getYDim() * pStruct.getZDim() )
    ret                                       = numpy.empty ( [ int ( pStruct.getXDim()), int ( pStruct.getYDim() ), int ( pStruct.getZDim() ) ], dtype = complex )
                
    for xIt in range ( 0, int ( pStruct.getXDim() ) ):
       for yIt in range ( 0, int ( pStruct.getYDim() ) ):
           for zIt in range ( 0, int ( pStruct.getZDim() ) ):
               index                          = int ( zIt + ( pStruct.getZDim() ) * ( yIt + ( pStruct.getYDim() ) * xIt ) )
               ret[xIt][yIt][zIt]             = realTrsFun[index] + 1j * imagTrsFun[index]
    
    return                                    ( ret )
 %}
