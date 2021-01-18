##############################################
##############################################
#   \file setup.py
#   \brief This file drives the python module installation using pip. It has been modified from https://www.benjack.io/2017/06/12/python-cpp-tests.html .
#
#   This file has been created by modification of the codes published in https://www.benjack.io/2017/06/12/python-cpp-tests.html . This file
#   is the driving code for pip installation of the ProSHADE modules. It defines the pip install Extension classes, which are given to setuptools
#   in order to allow building the module using CMake from pip. It also contains all the information supplied to the python module.
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
#   \version   0.7.5.1
#   \date      JAN 2021
##############################################
##############################################

##########################################################################################
##########################################################################################
##### Global settings
##########################################################################################
##########################################################################################
gl_version                                            = '0.7.5.1'
gl_download                                           = 'https://github.com/michaltykac/proshade/archive/v{0}.tar.gz'.format(gl_version)


##########################################################################################
##########################################################################################
##### Load required modules
##########################################################################################
##########################################################################################
import os
import re
import sys
import sysconfig
import platform
import subprocess

from distutils.version import LooseVersion
from distutils.dir_util import copy_tree
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

##########################################################################################
##########################################################################################
##### Create simple Extension module class
##########################################################################################
##########################################################################################
class CMakeExtension ( Extension ):
    def __init__ ( self, name, sourcedir='' ):
        Extension.__init__                            ( self, name, sources=[] )
        self.sourcedir                                = os.path.abspath ( sourcedir )

##########################################################################################
##########################################################################################
##### Create commands class to drive the CMake calls
##########################################################################################
##########################################################################################
class CMakeBuild ( build_ext ):
    """
    This is the function that will be called during the pip installation.
    
    Firstly, this function checks if cmake command can be found and then it checks if this cmake installation
    has at least the required version (3.4). If both conditions are met, then the function will call the
    building function to build using CMake.
    """
    def run ( self ):
        ### Check if CMake can be found
        try:
            out                                       = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError                        ( "Cannot find CMake on your system. Please install CMake version 3.4 or higher to allow for installation of ProSHADE." )

        ### Check for CMake version
        cmake_version                                 = LooseVersion ( re.search ( r'version\s*([\d.]+)', out.decode()).group(1) )
        if cmake_version < '3.4.0':
            raise RuntimeError                        ( "CMake version " + str( cmake_version ) + " found on your system, however, ProSHADE requires at least CMake version 3.4. Please consider updating your CMake." )
    
        ### Build the extensions
        for ext in self.extensions:
            self.build_extension                      ( ext )

    """
    """
    def build_extension ( self, ext ):
        ### Path to where the extension will be build
        extdir                                        = os.path.abspath( os.path.dirname ( self.get_ext_fullpath ( ext.name ) ) )
        
        ### CMake arguments to make sure the CMake output goes to the correct folder and that the correct python is used
        cmake_args                                    = ['-DBUILD_PYTHON=TRUE',
                                                         '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                                                         '-DPYTHON_EXECUTABLE=' + sys.executable ]

        ### Make arguments
        cfg                                           = 'Debug' if self.debug else 'Release'
        build_args                                    = ['--config', cfg]

### Windows specific arguments will be added when Windows is supported ...
#        if platform.system() == "Windows":
#            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
#                cfg.upper(),
#                extdir)]
#            if sys.maxsize > 2**32:
#                cmake_args += ['-A', 'x64']
#            build_args += ['--', '/m']
#        else:

        ### Platform spacifc CMake and Make commands
        if platform.system() != "Windows":
            cmake_args                               += [ '-DCMAKE_BUILD_TYPE=' + cfg ]
            build_args                               += [ '--', '-j2' ]

        ### Set the environment flags
        env                                           = os.environ.copy()
        env['CXXFLAGS']                               = '{} -DVERSION_INFO=\\"{}\\"'.format( env.get( 'CXXFLAGS', '' ), self.distribution.get_version() )
        
        ### If temp folder does not exist, create it
        if not os.path.exists ( self.build_temp ):
            os.makedirs                               ( self.build_temp )
    
        ### Run CMake
        subprocess.check_call                         ( ['cmake', os.path.join ( ext.sourcedir, 'proshade' )]  + cmake_args, cwd = self.build_temp, env = env )
        subprocess.check_call                         ( ['cmake', '--build', '.'] + build_args, cwd = self.build_temp )

##########################################################################################
##########################################################################################
##### Set README.md as long description
##########################################################################################
##########################################################################################
this_directory                                        = os.path.abspath ( os.path.dirname ( __file__ ) )
with open(os.path.join(this_directory, 'README.md')) as f:
    long_description                                  = f.read ( )

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

extra_files = package_files('proshade')

##########################################################################################
##########################################################################################
##### Module info
##########################################################################################
##########################################################################################
setup (
    name                                              = 'proshade',
    version                                           =  gl_version,
    author                                            = 'Michal Tykac, Garib N. Murshudov',
    author_email                                      = 'Michal.Tykac@gmail.com',
    url                                               = 'https://github.com/michaltykac/proshade',
    download_url                                      =  gl_download,
    description                                       = 'Protein Shape Description and Symmetry Detection (ProSHADE) python module',
    long_description                                  = long_description,
    long_description_content_type                     = 'text/markdown',
    ext_modules                                       = [ CMakeExtension ( 'proshade' ) ],
    cmdclass                                          = dict ( build_ext = CMakeBuild ),
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: C',
        'Programming Language :: C++',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords                                          = 'bioinformatics protein-shapes symmetry-detection computational-biology structural-biology',
    project_urls                                      = { 'Github': 'https://github.com/michaltykac/proshade',
                                                          'Bug Reports': 'https://github.com/michaltykac/proshade/issues' },
    setup_requires                                    = ['numpy','setuptools'],
    install_requires                                  = ['numpy','setuptools'],
    packages                                          = ['proshade'],
    package_data                                      = {'' : extra_files },
    zip_safe                                          = False,
)

##########################################################################################
##########################################################################################
##### Done
##########################################################################################
##########################################################################################
