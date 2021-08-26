######################################################
######################################################
#   \file plotRotationFunctionSpheres.py
#   \brief This file shows how ProSHADE can be used in conjunction with EMDA to plot the rotation function spheres
#
#   This code demonstrates how the self-rotation function computed for a particular structure can be mapped onto a set of concentric spheres with
#   the radius of a sphere is the angle of the axis-angle rotation representation, while the longitude and latitude represent the axis of the axis-
#   angle rotation representation. Furthermore, it shows how the
#
#   Copyright by Rangana Sanjeewa Warshamanage, Michal Tykac and individual contributors. All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#   1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#   2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#   3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
#   This software is provided by the copyright holders and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In     no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data     or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility     of such damage.
#
#   \author    Rangana Sanjeewa Warshamanage
#   \author    Michal Tykac
#   \author    Garib N. Murshudov
#   \version   0.1.1
#   \date      AUG 2021
######################################################
######################################################

######################################################
### Import modules
### ==============
###
### This is where Python modules are loaded.
###

### System modules
import sys

### Import numpy
import numpy

### Import ProSHADE
import proshade

######################################################
### Global settings
### ===============
###
### This is where all the settings are given.
###
structureFilename                                     = "/Users/mysak/BioCEV/proshade/playground/emd_6324.map"
structureFold                                         = 12
structureReSampleMap                                  = True
computationResolution                                 = 8

######################################################
### Local functions
### ===============
###
### This is where local functions are defined.
###
def plot3d ( data, shNo ):
    """
    This function plots a single sphere from the proshade supplied sphere mapped rotation
    function data.

    Parameters
    ----------
    matrix : data
        A square matrix indexed by longitude x latitude positions and containing the mapped
        rotation function values.
        
    int : shNo
        The number (index) of the shell being plotted.

    Returns
    -------
    none

    """
    ### Import matplotlib
    from matplotlib import cm
    import matplotlib.pyplot as plt

    ### Initialise local variables
    n, m                                              = data.shape

    ### Meshing a unit sphere according to n, m
    theta                                             = numpy.linspace ( 0,           2 * numpy.pi, num = n )
    phi                                               = numpy.linspace ( -numpy.pi/2, numpy.pi/2,   num = m )
    theta, phi                                        = numpy.meshgrid ( theta, phi )
    
    ### NOTE: phi angle follows MATLAB convention
    x, y, z                                           = numpy.cos ( phi ) * numpy.cos ( theta ), numpy.cos ( phi ) * numpy.sin ( theta ), numpy.sin ( phi )

    # Data normalization - to be done
    # normalized to min-max: does not give enough contrast in the figure
    #fmax, fmin                                        = data.max(), data.min()
    #data                                              = (data - fmin)/(fmax - fmin)
    #data                                              = data/numpy.mean(data) # normalize to mean for good contrast

    ### Prepare the projection
    fig                                               = plt.figure ( figsize = plt.figaspect( 1.0 ) )
    ax                                                = fig.gca ( projection = '3d' )
    ax.set_proj_type                                  ( 'ortho' )
    ax.plot_surface                                   ( x, y, z, rstride = 3, cstride = 3, facecolors=cm.Blues ( data ) )
    
    ### Colobar for facecolors
    mappable                                          = plt.cm.ScalarMappable ( cmap = cm.Blues )
    mappable.set_array                                ( data )
    plt.colorbar                                      ( mappable, fraction = 0.025, pad = 0.04 )
    
    ### Prepare the plot/box
    ax.view_init                                      ( azim = 0.0, elev = 90.0 )
    ax.set_box_aspect                                 ( ( 1, 1, 1 ) )
    ax.axis                                           ( 'off' )
    title                                             = 'ProShade Spherical harmonics: shell-%i' %shNo
    plt.title                                         ( title )
    
    ### Plot
    plt.show                                          ( )
    plt.close                                         ( 'all' )

######################################################
### Compute rotation function
### =========================
###
### This is the standard proshade code to compute the
### self-rotation function for a given structure.
###

### Create the settings object
pSet                                                  = proshade.ProSHADE_settings ( )

### Apply the user supplied settings
pSet.task                                             = proshade.Symmetry
pSet.setResolution                                    ( computationResolution )
pSet.setMapResolutionChange                           ( structureReSampleMap )
pSet.verbose                                          = 2

### Read in the structure
pStruct                                               = proshade.ProSHADE_data ( )
pStruct.readInStructure                               ( structureFilename, 0, pSet )
    
### Do all the required computations in correct order
pStruct.processInternalMap                            ( pSet )
pStruct.mapToSpheres                                  ( pSet )
pStruct.computeSphericalHarmonics                     ( pSet )
pStruct.computeRotationFunction                       ( pSet )

######################################################
### Obtain the sphere-mapped rotation function
### ==========================================
###
### This is how sphere mapped rotation function can
### be retrieved from proshade
###
sphereMappedData                                      = proshade.getRotationFunctionSpheres ( pStruct, structureFold )

######################################################
### Plot the spheres
### ==========================================
###
### This is where each sphere is individually plotted.
###
for sphereNo in range ( sphereMappedData.shape[0] ):
    plot3d                                            ( sphereMappedData[sphereNo,:,:], sphereNo )

######################################################
### Wait for matplotlib window closure
### ==================================
###
### This is to avoid error message when python tries
### to terminate while matplotlib window is open.
###


######################################################
### Done
### ====
###
### The code exits here
###
sys.exit                                              ( )
