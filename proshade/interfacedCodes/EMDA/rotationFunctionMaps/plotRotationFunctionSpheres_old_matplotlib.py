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
#   \version   0.1.2
#   \date      FEB 2022
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
#structureFilename                                     = "/Users/mysak/BioCEV/proshade/playground/emd_0020.map"
structureFilename                                     = "/Users/mysak/BioCEV/proshade/xx_EMDBSymmetry/EMD-0020/emd_0020.map.gz"
structureFold                                         = 2
structureReSampleMap                                  = True
computationResolution                                 = 20

######################################################
### Local functions
### ===============
###
### This is where local functions are defined.
###
def plot3dIndividually ( data, ncols = 6, stride = 5 ):
    """
    This function plots the rotation function mapped into axis-angle space data as returned by ProSHADE. This
    particular version plots all spheres separately as individual sub-plots of one large figure.

    Parameters
    ----------
    matrix : data
        A square matrix indexed by longitude x latitude positions and containing the mapped
        rotation function values.
        
    int : ncols
        The maximum allowed number of sub-plots in a row
        
    int : stride
        How fine sampled should the spheres surface be?

    Returns
    -------
    none

    """
    ### Import matplotlib
    from matplotlib import cm
    import matplotlib.pyplot as plt
    
    ### Figure how many rows and columns we need given that we can have up to 6 columns
    nrows                                             = divmod( data.shape[0], ncols )[0] + 1
    if data.shape[0] < ncols:
        ncols                                         = data.shape[0]

    ### Create figure with right aspect ratios
    fig                                               = plt.figure ( figsize = plt.figaspect ( nrows / ncols ) )

    ### For each sphere, create a sub-plot
    plotIndex                                         = 0
    for sphNo in range ( data.shape[0] ):

        ### Initialise local variables
        n, m                                          = data[plotIndex].shape
    
        ### Meshing a unit sphere according to n, m
        theta                                         = numpy.linspace ( 0,           2 * numpy.pi, num = n )
        phi                                           = numpy.linspace ( -numpy.pi/2, numpy.pi/2,   num = m )
        theta, phi                                    = numpy.meshgrid ( theta, phi )
        
        ### NOTE: phi angle follows MATLAB convention
        x, y, z                                       = numpy.cos ( phi ) * numpy.cos ( theta ), numpy.cos ( phi ) * numpy.sin ( theta ), numpy.sin ( phi )
    
        # Data normalization - to be done
        # normalized to min-max: does not give enough contrast in the figure
        #fmax, fmin                                    = data.max(), data.min()
        #data                                          = (data - fmin)/(fmax - fmin)
        #data                                          = data/numpy.mean(data) # normalize to mean for good contrast
    
        ### Prepare the projection
        ax                                            = fig.add_subplot ( nrows, ncols, plotIndex + 1, projection='3d' )
        ax.set_proj_type                              ( 'ortho' )
        ax.plot_surface                               ( x, y, z, rstride = stride, cstride = stride, facecolors=cm.Blues ( data[plotIndex] ) )
        
        ### Colobar for facecolors
        mappable                                      = plt.cm.ScalarMappable ( cmap = cm.Blues )
        mappable.set_array                            ( data[plotIndex] )
        plt.colorbar                                  ( mappable, fraction = 0.025, pad = 0.04 )
        
        ### Prepare the plot/box
        ax.view_init                                  ( azim = 0.0, elev = 90.0 )
        ax.set_box_aspect                             ( ( 1, 1, 1 ) )
        ax.axis                                       ( 'off' )
        
        ### Set subtitle
        ax.set_title                                  ( 'Angle ' + "{:.2f}".format( (plotIndex+1) * ((numpy.pi*2.0)/(data.shape[0]+1)) * ( 180.0 / numpy.pi ) ) + ' deg' )
        
        ### Shift iterator
        plotIndex                                     = plotIndex + 1
        
    ### Main title
    fig.suptitle                                      ( 'Self-rotation function mapped into axis-angle space for all axes with particular angles MASTER', fontsize=16 )
    
    ### Plot
    plt.tight_layout                                  ( )
    plt.show                                          ( )
    plt.close                                         ( 'all' )
    
def plot3dAverage ( data, ncols = 6, stride = 5 ):
    """
    This function plots the rotation function mapped into axis-angle space data as returned by ProSHADE. This
    particular version plots the average for each axis over all the angles appropriate for the required fold.

    Parameters
    ----------
    matrix : data
        A square matrix indexed by longitude x latitude positions and containing the mapped
        rotation function values.
        
    int : ncols
        The maximum allowed number of sub-plots in a row
        
    int : stride
        How fine sampled should the spheres surface be?

    Returns
    -------
    none

    """
    ### Import matplotlib
    from matplotlib import cm
    import matplotlib.pyplot as plt
    
    ### Sum the data
    summedData                                        = numpy.zeros ( [ data.shape[1], data.shape[2] ] )
    for sphNo in range ( data.shape[0] ):
        summedData                                   += data[sphNo]

    summedData                                       /= data.shape[0]

    ### Create figure with right aspect ratios
    fig                                               = plt.figure ( figsize = plt.figaspect ( 1.0 ) )

    ### Initialise local variables
    n, m                                              = summedData.shape
    
    ### Meshing a unit sphere according to n, m
    theta                                             = numpy.linspace ( 0,           2 * numpy.pi, num = n )
    phi                                               = numpy.linspace ( -numpy.pi/2, numpy.pi/2,   num = m )
    theta, phi                                        = numpy.meshgrid ( theta, phi )
        
    ### NOTE: phi angle follows MATLAB convention
    x, y, z                                           = numpy.cos ( phi ) * numpy.cos ( theta ), numpy.cos ( phi ) * numpy.sin ( theta ), numpy.sin ( phi )
    
    # Data normalization - to be done
    # normalized to min-max: does not give enough contrast in the figure
    #fmax, fmin                                    = data.max(), data.min()
    #data                                          = (data - fmin)/(fmax - fmin)
    #data                                          = data/numpy.mean(data) # normalize to mean for good contrast
    
    ### Prepare the projection
    ax                                                = fig.add_subplot ( 1, 1, 1, projection='3d' )
    ax.set_proj_type                                  ( 'ortho' )
    ax.plot_surface                                   ( x, y, z, rstride = stride, cstride = stride, facecolors=cm.Blues ( summedData ) )
        
    ### Colobar for facecolors
    mappable                                          = plt.cm.ScalarMappable ( cmap = cm.Blues )
    mappable.set_array                                ( summedData )
    plt.colorbar                                      ( mappable, fraction = 0.025, pad = 0.04 )
        
    ### Prepare the plot/box
    ax.view_init                                      ( azim = 0.0, elev = 90.0 )
    ax.set_box_aspect                                 ( ( 1, 1, 1 ) )
    ax.axis                                           ( 'off' )
        
    ### Set subtitle
    ax.set_title                                      ( 'Sum over all angles required for this fold' )
        
    ### Main title
    fig.suptitle                                      ( 'Sum of self-rotation function mapped into axis-angle space for all axes', fontsize=16 )
    
    ### Plot
    plt.show                                          ( )
    plt.close                                         ( 'all' )
    
def plot3dOverlay ( data, ncols = 6, stride = 5, alpha = 0.7 ):
    """
    This function plots the rotation function mapped into axis-angle space data as returned by ProSHADE. This
    particular version plots ???

    Parameters
    ----------
    matrix : data
        A square matrix indexed by longitude x latitude positions and containing the mapped
        rotation function values.
        
    int : ncols
        The maximum allowed number of sub-plots in a row
        
    int : stride
        How fine sampled should the spheres surface be?

    Returns
    -------
    none

    """
    ### Import matplotlib
    from matplotlib import cm
    import matplotlib.pyplot as plt

    ### Create figure with right aspect ratios
    fig                                               = plt.figure ( figsize = plt.figaspect ( 1.0 ) )

    ### For each sphere, create a sub-plot
    plotIndex                                         = 0
    for sphNo in range ( data.shape[0] ):

        ### Initialise local variables
        n, m                                          = data[plotIndex].shape
    
        ### Meshing a unit sphere according to n, m
        theta                                         = numpy.linspace ( 0,           2 * numpy.pi, num = n )
        phi                                           = numpy.linspace ( -numpy.pi/2, numpy.pi/2,   num = m )
        theta, phi                                    = numpy.meshgrid ( theta, phi )
        
        ### NOTE: phi angle follows MATLAB convention
        x, y, z                                       = numpy.cos ( phi ) * numpy.cos ( theta ), numpy.cos ( phi ) * numpy.sin ( theta ), numpy.sin ( phi )
    
        # Data normalization - to be done
        # normalized to min-max: does not give enough contrast in the figure
        #fmax, fmin                                    = data.max(), data.min()
        #data                                          = (data - fmin)/(fmax - fmin)
        #data                                          = data/numpy.mean(data) # normalize to mean for good contrast
    
        ### Prepare the projection
        ax                                            = fig.add_subplot ( projection='3d' )
        ax.set_proj_type                              ( 'ortho' )
        ax.plot_surface                               ( x, y, z, rstride = stride, cstride = stride, facecolors=cm.Blues ( data[plotIndex] ), alpha = alpha )
        
        ### Prepare the plot/box
        ax.view_init                                  ( azim = 0.0, elev = 90.0 )
        ax.set_box_aspect                             ( ( 1, 1, 1 ) )
        ax.axis                                       ( 'off' )
        
        ### Shift iterator
        plotIndex                                     = plotIndex + 1
        
    ### Main title
    fig.suptitle                                      ( 'Self-rotation function mapped into axis-angle space for all axes with particular angles overlay', fontsize=16 )
    
    ### Plot
    plt.tight_layout                                  ( )
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
pSet.verbose                                          = 3

### Read in the structure
pStruct                                               = proshade.ProSHADE_data ( )
pStruct.readInStructure                               ( structureFilename, 0, pSet )

pSet.addExtraSpace                                    = 10.0
    
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
### This is where each spheres are plotted.
###

plot3dIndividually                                    ( sphereMappedData, stride = 1 )
#plot3dAverage                                         ( sphereMappedData )
#plot3dOverlay                                         ( sphereMappedData )

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
