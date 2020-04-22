# - Find the FFTW3 library
#
# This cmake module attempts to find the FFTW3 library and header files. 
#
# It sets the following variables:
#   FFTW_FOUND                  ===> true if FFTW3 is found on the system
#   FFTW_LIBRARIES              ===> full path to the folder containing libfftw3.a/so/dylib/...
#   FFTW_INCLUDES               ===> full path to the folder containing fftw3.h header file
#   FFTW_DEFINITIONS            ===> these should be compiler flags required to link against this 
#                                    library, but it will be empty as I am not sure how to make this work...
#
#  This file is part of the ProSHADE library for calculating
# shape descriptors and symmetry operators of protein structures.
# This is a prototype code, which is by no means complete or fully
# tested. Its use is at your own risk only. There is no quarantee
# that the results are correct.
# 
# Author    Michal Tykac
# Author    Garib N. Murshudov
#
##########################################################################################

##########################################################################################
### Set the FFTW_ROOT if not already defined and if available
if( NOT FFTW_ROOT AND DEFINED ENV{FFTWDIR} )
  set   ( FFTW_ROOT $ENV{FFTWDIR}                                                         )
endif()

##########################################################################################
### Save the list of library suffixes for this system
set     ( SEARCH_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES}                          )

##########################################################################################
### Set the hints and paths to the standard values
set     ( FFTW_HINTS ${FFTW_ROOT} $ENV{FFTW_ROOT} ${CUSTOM_FFTW3_LIB_PATH} ${CUSTOM_FFTW3_INC_PATH} )
set     ( FFTW_PATHS /usr /usr/local /opt /opt/local ${FFTW_ROOT}                         )


##########################################################################################
### Find the library and header
foreach    ( SUFFIX ${SEARCH_LIBRARY_SUFFIXES} )
  find_path   ( FFTW_LIBRARIES
               NAMES "libfftw3${SUFFIX}"
               HINTS ${FFTW_HINTS}
               PATH_SUFFIXES "lib" "lib64" "lib/x86_64" "lib/x64" "lib/x86" "lib/x86_64-linux-gnu"
               PATHS ${FFTW_PATHS}
               DOC "FFTW3 library full path (libfftw3.a/so/dyblib)"
               NO_DEFAULT_PATH
             )
  if    ( FFTW_LIBRARIES )
    break (                                                                               )
  endif ( FFTW_LIBRARIES )
endforeach ( SUFFIX ${SEARCH_LIBRARY_SUFFIXES} )

find_path    ( FFTW_INCLUDES
               NAMES "fftw3.h"
               HINTS ${FFTW_HINTS}
               PATH_SUFFIXES "include" "api" "inc" "include/x86_64" "include/x64"
               PATHS ${FFTW_PATHS}
               DOC "FFTW3 include full path (fftw3.h)"
               NO_DEFAULT_PATH
             )
             
##########################################################################################
### Check if everything was found, otherwise give error message
if    ( NOT FFTW_LIBRARIES )
  message ( ERROR " Failed to locate the FFTW3 library. Please install FFTW3 to a standard location or use the cmake -DCUSTOM_FFTW3_LIB_PATH=/path/to/libfftw3.a/so/dylib command line argument to supply the paths to the location of the libfftw3.a/so/dylib library file." )
endif ( NOT FFTW_LIBRARIES )

if    ( NOT FFTW_INCLUDES )
  message ( ERROR " Failed to locate the FFTW3 header. Please install FFTW3 to a standard location or use the cmake -DCUSTOM_FFTW3_INC_PATH=/path/to/fftw3.h command line argument to supply the path to the location of the fftw3.h header file." )
endif ( NOT FFTW_INCLUDES )

##########################################################################################
### Deal with the standard arguments 
include ( FindPackageHandleStandardArgs                                                   )
find_package_handle_standard_args ( FFTW REQUIRED_VARS FFTW_LIBRARIES FFTW_INCLUDES       )

##########################################################################################
### Set the output variables to be accessible by the caller.
mark_as_advanced ( FFTW_LIBRARIES FFTW_INCLUDES FFTW_FOUND FFTW_DEFINITIONS )
