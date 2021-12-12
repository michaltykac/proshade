######################################################
######################################################
#   \file howto_accessInterimResults.py
#   \brief This code demonstrates how interim results of internal ProSHADE computation can be accessed from python.
#
#   This file starts by reading in a mask from a file and creating a Fourier weights array as an array of ones. These
#   arrays are for demonstration purposes, but te user could obtain different arrays or obtain different arrays and
#   use such arrays instead.
#
#   Next, the settings object is created and set followed by the structure being read from file. The commented lines
#   show how the structure could be read if no mask, no weights or no arrays were required by the user. Note that the
#   arrays will be applied only if the input file format is density map (and will be ignored if the input file format
#   is co-ordinate). Also, the arrays are not required to have the same dimensions as the read in map, but it is strongly
#   recommended as the Fourier space re-sampling of input arrays may introduce artifacts that could cause extra issues
#   with the results.
#
#   Finally, the code shows how the spherical harmonics, normalised E-matrices, SO(3) coefficients and the (self-)rotation
#   function can be accessed.
#
#   Copyright by Michal Tykac and individual contributors. All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#   1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#   2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#   3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
#   This software is provided by the copyright holders and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In     no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data     or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility     of such damage.
#
#   \author    Michal Tykac
#   \author    Garib N. Murshudov
#   \version   0.7.6.2
#   \date      DEC 2021
######################################################
######################################################

######################################################
### Import system modules
import sys
import numpy

######################################################
### Import ProSHADE from system folder
import proshade

######################################################
### Create mask and weights array. Do not need to have
### the same shape of map, but preferred.
import mrcfile
with mrcfile.open('/Users/mysak/BioCEV/proshade/playground/emd_0116_mask.map') as mrc:
    maskArr = mrc.data
    mrc.close()

weightsArr = numpy.ones( ( maskArr.shape[0], maskArr.shape[1], int(maskArr.shape[2]/2) ), dtype='double' )

######################################################
### Create the settings object
pSet                                                  = proshade.ProSHADE_settings ( proshade.Symmetry )

######################################################
### Set basic settings values
pSet.requestedResolution                              = 20.0  ## This is to make the example run faster as here we are not really interested in accuracy of results
pSet.changeMapResolution                              = True  ## but rather in showing how to use ProSHADE.
pSet.verbose                                          = -1

######################################################
### Read in the map and process it for RF calculation
pStruct                                               = proshade.ProSHADE_data ( )
pStruct.readInStructure                               ( "/Users/mysak/BioCEV/proshade/playground/emd_0116.map.gz", 0, pSet, maskArr, weightsArr )

### No mask - the empty array is required in this case
#pStruct.readInStructure                               ( "/Users/mysak/BioCEV/proshade/playground/emd_0116.map.gz", 0, pSet, numpy.ndarray ( 0 ), weightsArr )

### No weights - the empty array is optional in this case
#pStruct.readInStructure                               ( "/Users/mysak/BioCEV/proshade/playground/emd_0116.map.gz", 0, pSet, maskArr )
#   OR
#pStruct.readInStructure                               ( "/Users/mysak/BioCEV/proshade/playground/emd_0116.map.gz", 0, pSet, maskArr, numpy.ndarray ( 0 ) )

### No arrays - no arrays are required, but if weights empty array is given, the mask empty array must also be given
#pStruct.readInStructure                               ( "/Users/mysak/BioCEV/proshade/playground/emd_0116.map.gz", 0, pSet )
#pStruct.readInStructure                               ( "/Users/mysak/BioCEV/proshade/playground/emd_0116.map.gz", 0, pSet, numpy.ndarray ( 0 ) )
#pStruct.readInStructure                               ( "/Users/mysak/BioCEV/proshade/playground/emd_0116.map.gz", 0, pSet, numpy.ndarray ( 0 ), numpy.ndarray ( 0 ) )

pStruct.processInternalMap                            ( pSet )
pStruct.mapToSpheres                                  ( pSet )
pStruct.computeSphericalHarmonics                     ( pSet )

######################################################
### Accessing spherical harmonics
### =============================
###
### In order to access the spherical harmonics
### values for each sphere, please use the
### following functions:
###
### 1) getSphericalHarmonics()
###     This function returns a 2D numpy.ndarray of
###     complex numbers. This 2D array is organised so
###     that the first dimension is the sphere number
###     and the second dimension is a 1D array of all
###     the spherical harmonics values. This 1D array
###     is complicated by the fact that each sphere
###     can have different number of bands (depending
###     on the settings.progressiveSphereMapping value).
###     This means that some values for some spheres
###     are empty. Best way of accessing this array is
###     to use the second function.
###
### 2) findSHIndex()
###     This function takes the shell, band and order
###     values and returns the correct index for these
###     in the 1D array. Please note that the order
###     value needs to be in range -band <= order <= band
###     and NOT in the 0 <= order <= (2 * band) + 1.
###
sphericalHarmonics                                    = pStruct.getSphericalHarmonics ( )

shell =  3
band  =  4
order = -2
Shell3Band4OrderMin2Value                             = sphericalHarmonics[shell][pStruct.findSHIndex(shell, band, order)]
print                                                 ( Shell3Band4OrderMin2Value )
# Expected output: (-347.47243493408905+99.22184548560784j)

######################################################
### Computing self-rotation function
### ================================
###
### This function computes the self-rotation
### function by firstly computing and normalising
### the E matrices, then combining these into
### the SO(3) coefficients and finally calcula-
### ting the inverse SO(3) Fourier Transform
### (SOFT) from them.
###
pStruct.computeRotationFunction                       ( pSet )

######################################################
### Accessing E Matrices
### ====================
###
### ProSHADE allows access to the weighted E matrices
### ( Integral _0 ^rMAX ( c1^lm * c2*^lm' ) of
### the structure combination. These are
### returned as a 3D Numpy array with indices being
### band of the E matrix, order1 of the E matrix
### and order2 of the E matrix. Note that because
### indices need to go from zero, the order of indexing
### goes like, for example, this:
###
### BAND = 2 || ORDER = -2  || INDEX [2][0]
### BAND = 2 || ORDER = -1  || INDEX [2][1]
### BAND = 2 || ORDER =  0  || INDEX [2][2]
### BAND = 2 || ORDER =  1  || INDEX [2][3]
### BAND = 2 || ORDER =  2  || INDEX [2][4]
### i.e. order index = ORDER + BAND
###
### NOTE: As Numpy arrays have single shape, the
### lower bands (which will have less orders)
### are padded with zeroes to have the same length
### as the largest band.
###
eMat                                                  = pStruct.getEMatrix ( )

Band4OrderOneMin2OrderTwo3EMatrixValue                = eMat[4][2][7] # Band = 4, Order1 = -2 and Order2 = 3

print ( Band4OrderOneMin2OrderTwo3EMatrixValue )
# Expected output: (6.990090606803223e-06-3.1783258844641986e-06j)

######################################################
### Accessing SO(3) coefficients
### ============================
###
### ProSHADE also allows access to the SO(3)
### coefficients computed by normalising the E
### matrix values and dealing with the signs.
### The inverse SO(3) Fourier Transform (SOFT)
### of these values then results in the rotation
### function. The complete complex array of
### these values can be accessed as shown.
###
### NOTE: Given that the SO(3) coefficients have the same
### indexing as the E matrices, ProSHADE provides them
### to the user in the same format as the E matrices.
### This does, however, mean that the same caveats do
### apply to the SO(3) coefficients array as did to the
### E matrix array.
###
so3Coeffs                                             = pStruct.getSO3Coefficients ( )


Band4OrderOneMin2OrderTwo3SO3CoeffsValue              = so3Coeffs[4][2][7] # Band = 4, Order1 = -2 and Order2 = 3

print ( Band4OrderOneMin2OrderTwo3SO3CoeffsValue )
# Expected output: (-2.0704102862098065e-05+9.413953229328724e-06j)

######################################################
### Accessing self-rotation function
### ================================
###
### ProSHADE also gives access to the self -
### rotation function as shown next. The returned 3D complex
### array has dimensions 2 * bandwidth. Please note that
### the indices have nothing to do with the angle values,
### if you want to know the Euler angle values for a
### particular index, you need to convert it
### yourself. Alternatively (and a recommended
### approach is), you can use the
### getRotationMatrixFromEulerIndices() function
### also demonstrated here.
###
selfRotationFunction                                  = pStruct.getRotationFunctionMap ( )

### Find the map maximum
rotMapMax                                             = numpy.where ( selfRotationFunction == numpy.amax ( selfRotationFunction ) )

### Find maximum value
print                                                 ( "Rotation map maximum is: " + str( selfRotationFunction[rotMapMax[0][0]][rotMapMax[1][0]][rotMapMax[2][0]] ) )

### Expected output: Rotation map maximum is: (0.9922541180641358+6.5382598050851006e-21j)

### Find rotation matrix for the maximum (uses only the real parts, using magnitudes would be better)
rotMatMaxVal                                          = pStruct.getRotationMatrixFromSOFTCoordinates ( rotMapMax[0][0], rotMapMax[1][0], rotMapMax[2][0] )
print                                                 ( rotMatMaxVal )

### Expected output: [[ 1.00000000e+00 -4.66973965e-16  0.00000000e+00]
### Expected output:  [ 4.66973965e-16  1.00000000e+00 -0.00000000e+00]
### Expected output:  [ 0.00000000e+00  0.00000000e+00  1.00000000e+00]]

######################################################
### Release ProSHADE memory
del pStruct
del pSet

######################################################
### Done
