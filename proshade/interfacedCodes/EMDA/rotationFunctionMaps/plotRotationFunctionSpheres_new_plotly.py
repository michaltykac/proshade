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
#   \version   0.3.0
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

### Import plotting libraries
from turtle import color
import math
import plotly.graph_objects as go
import plotly.offline as offline

######################################################
### Global settings
### ===============
###
### This is where all the settings are given.
###
structureFilename                                     = "/Users/mysak/BioCEV/proshade/playground/emd_10575.map"
structureFold                                         = 2
structureReSampleMap                                  = True
computationResolution                                 = 4
verbosity                                             = 3

######################################################
### Compute rotation function
### =========================
###
### This is the standard proshade code to compute the
### self-rotation function for a given structure and
### to return the sphere mapped version of it.
###
def getSphereMappedRF ( strName, strFold, res = 8.0, reSample = True, verbosity = 3 ):
    """
    This function contains all the ProSHADE calls and works required to compute the sphere
    mapped rotation function values for a particular structure (first argument), for a particular
    symmetry fold (second argument) and given some other optional arguments.

    Parameters
    ----------
    strName : string
        This is the path to the structure file for the structure to have the sphere mapped RF computed.
        
    strFold : integer
        The fold (or rotation angles) for which the spheres should be computed.
       
    res : float
        The resolution to which the computation is to be done up to.
       
    reSample : boolean
        Should the density map obtained from the structure file be re-sampled to match the requested
        resolution (True) or should the original data be used in their full size (False)?
        
    verbosity : integer
        How loud should the ProSHADE run be?

    Returns
    -------
    sphereMappedData : numpy.ndarray
        A 3D array of real numbers, where the first index is the sphere ID (there will always be strFold-1
        spheres), while the second and third indices are the lattitude and longitude indices on the surface
        of a sphere. These three indices will then refer to a value of the rotation function when mapped onto
        a sphere with radius equal to the angle of rotation and the latitude and longitude representing the
        rotation axis.

    """
    ### Create the settings object
    pSet                                              = proshade.ProSHADE_settings ( )
    
    ### Apply the user supplied settings
    pSet.task                                         = proshade.Symmetry
    pSet.setResolution                                ( res )
    pSet.setMapResolutionChange                       ( reSample )
    pSet.verbose                                      = verbosity
    pSet.addExtraSpace                                = 10.0
    
    ### Read in the structure
    pStruct                                           = proshade.ProSHADE_data ( )
    pStruct.readInStructure                           ( strName, 0, pSet )
        
    ### Do all the required computations in correct order
    pStruct.processInternalMap                        ( pSet )
    pStruct.mapToSpheres                              ( pSet )
    pStruct.computeSphericalHarmonics                 ( pSet )
    pStruct.computeRotationFunction                   ( pSet )
    
    ##################################################
    ### Obtain the sphere-mapped rotation function
    ### ==========================================
    ###
    ### This is how sphere mapped rotation function can
    ### be retrieved from proshade
    ###
    sphereMappedData                                  = proshade.getRotationFunctionSpheres ( pStruct, strFold )
    return                                            ( sphereMappedData )

######################################################
### 3D Plotting code
### ================
###
### This is the code required to plot the interactive
### 3D plots of the sphere mapped rotation function.
###
def axes(axis, vmax, clr='black', wdth=3):
    vmin = -vmax
    if axis.lower() == 'x':
        xcrd=[i*0.1 for i in range(vmin, vmax, 1)]
        ycrd=[i*0. for i in range(vmin, vmax, 1)]
        zcrd=[i*0. for i in range(vmin, vmax, 1)]
    elif axis.lower() == 'y':
        xcrd=[i*0. for i in range(vmin, vmax, 1)]
        ycrd=[i*0.1 for i in range(vmin, vmax, 1)]
        zcrd=[i*0. for i in range(vmin, vmax, 1)]
    elif axis.lower() == 'z':
        xcrd=[i*0. for i in range(vmin, vmax, 1)]
        ycrd=[i*0. for i in range(vmin, vmax, 1)]
        zcrd=[i*0.1 for i in range(vmin, vmax, 1)]
    else:
        raise SystemExit("Invalid axis label")

    trace = go.Scatter3d(x=xcrd, y=ycrd, z=zcrd, marker=dict(size=0.1), line=dict(color=clr,width=wdth))
    return trace

def annot(xcrd, ycrd, zcrd, txt, xancr='left'):
    strng=dict(showarrow=False, x=xcrd, y=ycrd, z=zcrd, text=txt, xanchor=xancr, font=dict(color='black',size=16))
    return strng

def plot3d_plotly(data):
    '''
    plot data in spherical coordinates
    Inspiration:
    https://python.plainenglish.io/how-to-create-a-3d-model-of-the-solar-system-with-plotly-in-python-2a74e672b771
    '''
    from plotly.offline import iplot
    

    n, m = data.shape

    # Meshing a unit sphere according to n, m
    theta = numpy.linspace(0, 2 * numpy.pi, num=n)
    phi = numpy.linspace(-numpy.pi/2, numpy.pi/2, num=m)
    theta, phi = numpy.meshgrid(theta, phi)
    # NOTE: phi angle follows MATLAB convention
    x, y, z = numpy.cos(phi)*numpy.cos(theta), numpy.cos(phi)*numpy.sin(theta), numpy.sin(phi)

    trace1 = go.Surface(z=z, x=x, y=y, surfacecolor=data)
    vmax = 20
    trace_x = axes("X", vmax, clr="red")
    trace_y = axes("Y", vmax, clr="green")
    trace_z = axes("Z", vmax, clr="blue")
    layout = go.Layout(
        showlegend=False,
        scene = dict(xaxis=dict(#title='X',
                                #titlefont_color='black',
                                #range=[-1.5,1.5],
                                backgroundcolor='white',
                                color='white',
                                gridcolor='white',
                                showgrid=False
                                ),
                    yaxis=dict(#title='Y',
                               #titlefont_color='black',
                               #range=[-1.5,1.5],
                               backgroundcolor='white',
                               color='white',
                               gridcolor='white',
                               showgrid=False
                               ),
                    zaxis=dict(#title='Z',
                               #titlefont_color='black',
                               #range=[-1.5,1.5],
                               backgroundcolor='white',
                               color='white',
                               gridcolor='white',
                               showgrid=False
                               ),
                    annotations=[
                                annot(vmax/10, 0, 0, 'X'),
                                annot(0, vmax/10, 0, 'Y'),
                                annot(0, 0, vmax/10, 'Z'),
                                ]
                               ))

    fig = go.Figure(data=[trace1, trace_x, trace_y, trace_z], layout=layout)
    plot = offline.plot(fig, filename="plot.html")



######################################################
### Get data and plot them
### ======================
###
### This is where it all comes together.
###
if __name__ == "__main__":

    sphereMappedData                                  = getSphereMappedRF ( structureFilename, structureFold, computationResolution, structureReSampleMap, verbosity )
    for i in range ( sphereMappedData.shape[0] ):
        plot3d_plotly                                 ( sphereMappedData[i,:,:] )
    

######################################################
### Done
### ====
###
### The code exits here
###
sys.exit                                              ( )
