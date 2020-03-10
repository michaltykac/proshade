##########################################################################################
##########################################################################################
################################### ProSHADE cmake file ##################################
##########################################################################################
##########################################################################################

##########################################################################################
################################### Decide version of CMake and use appropriate fn calls
set ( PRODUCE_PYTHON_BINDINGS   “0”                                                       )
set ( USE_MODERN_PATHS          “0”                                                       )

if    ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 6 )
	message ( FATAL_ERROR "Your CMake version is less than 2.6. ProSHADE requires a newer version and cannot be installed using the current version. Please update your CMake and try again." )
endif ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 6 )

if    ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 5 )
	if    ( ${CMAKE_MINOR_VERSION} GREATER 7 )
		message ( "Your CMake version is lower than recommended (3.1). This will result in no Python language bindings and slower compilation. Please consider updating your CMake." )
		set ( PRODUCE_PYTHON_BINDINGS   “0”                                               )
		set ( USE_MODERN_PATHS          “1”                                               )
	else  ( ${CMAKE_MINOR_VERSION} GREATER 7 )
		message ( "Your CMake version is lower than recommended (3.1). This will result in no Python language bindings and slower compilation. Please consider updating your CMake." )
		set ( PRODUCE_PYTHON_BINDINGS   “0”                                               )
		set ( USE_MODERN_PATHS          “0”                                               )
	endif ( ${CMAKE_MINOR_VERSION} GREATER 7 )
endif ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 5 )

if    ( ${CMAKE_MAJOR_VERSION} GREATER 2 )
	if    ( ${CMAKE_MINOR_VERSION} LESS 1 )
		message ( "Your CMake version is lower than recommended (3.1). This will result in no Python language bindings and slower compilation. Please consider updating your CMake." )
		set ( PRODUCE_PYTHON_BINDINGS   “0”                                               )
		set ( USE_MODERN_PATHS          “1”                                               )
	else  ( ${CMAKE_MINOR_VERSION} LESS 1 )
		message ( "Your CMake version is at least the recommended version. Allowing fast build and Python language binding to be created." )
		set ( PRODUCE_PYTHON_BINDINGS   “1”                                               )
		set ( USE_MODERN_PATHS          “1”                                               )
	endif ( ${CMAKE_MINOR_VERSION} LESS 1 )
endif ( ${CMAKE_MAJOR_VERSION} GREATER 2 )

##########################################################################################
################################### Set the correct required version
if    ( ${PRODUCE_PYTHON_BINDINGS} MATCHES 1 )
	cmake_minimum_required  ( VERSION 3.1                                                 )
else  ( ${PRODUCE_PYTHON_BINDINGS} MATCHES 1 )
	cmake_minimum_required  ( VERSION 2.6                                                 )
endif ( ${PRODUCE_PYTHON_BINDINGS} MATCHES 1 )

##########################################################################################
################################### Set project name
project ( proshade VERSION 0.7.1 LANGUAGES CXX                                            )     

##########################################################################################
################################### Disable C++11
set ( CMAKE_CXX_STANDARD 98                                                               )
set ( CMAKE_CXX_STANDARD_REQUIRED ON                                                      )

##########################################################################################
################################### Find the source files
file    ( GLOB SOURCES  "${CMAKE_SOURCE_DIR}/src/proshade/*.cpp"                       )
file    ( GLOB EXEC_SRC "${CMAKE_SOURCE_DIR}/src/bin/*.cpp"                            )

##########################################################################################
################################### Set flags depending on compiler

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  ################################# Using Clang
  set ( CMAKE_CXX_FLAGS "-O3 -g -Wshadow -Wpedantic -Wall -Wno-macro-redefined"           )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  ################################# Using GCC
  set ( CMAKE_CXX_FLAGS "-O3 -g -Wshadow -Wpedantic -Wall -Wno-c++11-compat"                                )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  ################################# Using Intel C++
  ### TO BE COMPLETED
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  ################################# Using Visual Studio C++
  ### TO BE COMPLETED
endif()

##########################################################################################
################################### Add clipper dependency

### Set command line accessible variable for clipper path
set ( CLIPPER_INC    "/Users/mysak/Repositories/CCP4_devtools_clang/install/include"      CACHE STRING "path to clipper/clipper.h" )
set ( CLIPPER_LIB    "/Users/mysak/Repositories/CCP4_devtools_clang/install/lib"          CACHE STRING "path to libclipper-core.*" )
set ( CCP4_CMAP_INC  "/Users/mysak/Repositories/CCP4_devtools_clang/install/include/ccp4" CACHE STRING "path to cmaplib.h"         )

if    ( EXISTS ${CLIPPER_INC} AND EXISTS ${CLIPPER_LIB} AND EXISTS ${CCP4_CMAP_INC} )
	message             ( STATUS "Clipper library path is      : " ${CLIPPER_LIB}         )
	message             ( STATUS "Clipper include path is      : " ${CLIPPER_INC}         )
	message             ( STATUS "CCP4 CMapLib include path is : " ${CCP4_CMAP_INC}       )
	include_directories ( ${CLIPPER_INC}                                                  )
	include_directories ( ${CCP4_CMAP_INC}                                                )
	link_directories    ( ${CLIPPER_LIB}                                                  )            
	set                 ( CLIPPER_LINK ${CLIPPER_LIB} CACHE STRING "RPATH value to clipper libs" )
else  ( EXISTS ${CLIPPER_INC} AND EXISTS ${CLIPPER_LIB} AND EXISTS ${CCP4_CMAP_INC} )
	message             ( FATAL_ERROR "Could not find the Clipper dependency. Please use the following command line arguments to supply the required paths: \n   -DCLIPPER_INC=/path/to/CCP4_devtools/install/include\n   -DCLIPPER_LIB=/path/to/CCP4_devtools/install/lib\n   -DCCP4_CMAP_INC=/path/to/CCP4_devtools/install/include/ccp4" )
endif ( EXISTS ${CLIPPER_INC} AND EXISTS ${CLIPPER_LIB} AND EXISTS ${CCP4_CMAP_INC} )

##########################################################################################
################################### Add fftw3 dependency

### Set command line accessible variable for clipper path
set ( FFTW3_LIB_PATH "/qusr/local/lib"     CACHE STRING "path to libfftw3.a"               )
set ( FFTW3_INC_PATH "/qusr/local/include" CACHE STRING "path to fftw3.h"                  )

if    ( EXISTS ${FFTW3_LIB_PATH} AND EXISTS ${FFTW3_INC_PATH} )
	message                ( STATUS "FFTW3 include path is        : " ${FFTW3_INC_PATH}   )
	include_directories    ( ${FFTW3_INC_PATH}                                            )
	message                ( STATUS "FFTW3 library path is        : " ${FFTW3_LIB_PATH}   )
	link_directories       ( ${FFTW3_LIB_PATH}                                            )
	set                    ( FFTW_LINK ${FFTW3_LIB_PATH} CACHE STRING "RPATH value to FFTW libs" )
else  ( EXISTS ${FFTW3_LIB_PATH} AND EXISTS ${FFTW3_INC_PATH} )
	message             ( FATAL_ERROR "Could not find the FFTW3 dependency. Please use the following command line arguments to supply the required paths: \n   -DFFTW3_LIB_PATH=/path/to/libfftw3.a\n   -DFFTW3_INC_PATH=/path/to/fftw3.h" )
endif ( EXISTS ${FFTW3_LIB_PATH} AND EXISTS ${FFTW3_INC_PATH} )

##########################################################################################
################################### Add SOFT dependency

### Compile Soft-2.0 from supplied folder
add_custom_target (
   soft2.0_lib
   COMMAND make lib FFTWINC=-I${FFTW3_INC_PATH} FFTWLIB1=-L${FFTW3_LIB_PATH} FFTWLIB2=-lfftw3
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/soft-2.0
)

### Report SOFT related info
message                ( STATUS "SOFT library path is         : " ${CMAKE_SOURCE_DIR}/extern/soft-2.0 )
message                ( STATUS "SOFT include path is         : " ${CMAKE_SOURCE_DIR}/extern/soft-2.0/include )
include_directories ( ${CMAKE_SOURCE_DIR}/extern/soft-2.0/include                         )

##########################################################################################
################################### Add LAPACK dependency
 
### Automatic LAPACK detection
find_package            ( lapack REQUIRED                                                 )

### Check for success
if    ( LAPACK_FOUND )    
	link_directories    ( ${LAPACK_LIBRARIES}                                             )
	set                 ( LAPACK_LINK ${LAPACK_LIBRARIES} CACHE STRING "RPATH value to LAPACK libs" )
else  ( LAPACK_FOUND )
    message             ( FATAL_ERROR "Could not find the LAPACK dependency. Please install LAPACK to a standard location or contact the author for help." )
endif ( LAPACK_FOUND )

##########################################################################################
################################### Set installation folders
set ( INSTALL_LOCALLY ON CACHE BOOL "Should installation be done locally?" )
if    ( INSTALL_LOCALLY )
	set ( MY_INSTALL_LOCATION ${CMAKE_SOURCE_DIR}/install CACHE STRING "Installation path" )
	set ( MY_PY2_INSTALL_LOCATION ${CMAKE_SOURCE_DIR}/install/python2 CACHE STRING "Python2 install location" )
	set ( MY_PY3_INSTALL_LOCATION ${CMAKE_SOURCE_DIR}/install/python3 CACHE STRING "Python3 install location" )
else  ( INSTALL_LOCALLY )
	set ( MY_INSTALL_LOCATION ${CMAKE_INSTALL_PREFIX} CACHE STRING "Installation path"    )
	
	execute_process      ( COMMAND python -c "import sys; print ( str ( sys.version_info[0] ) + '.' + str ( sys.version_info[1] ) )" OUTPUT_VARIABLE PY_VERSION )
	string ( REGEX REPLACE "\n$" "" PY_VERSION "${PY_VERSION}"                            )
	execute_process      ( COMMAND python -c "import sys; print ( sys.exec_prefix )" OUTPUT_VARIABLE PY_PREFIX )
	string ( REGEX REPLACE "\n$" "" PY_PREFIX "${PY_PREFIX}"                              )
	string ( CONCAT PYTHON_INSTALL ${PY_PREFIX} "/lib/python"                             )
	string ( REGEX REPLACE "\n$" "" PYTHON_INSTALL "${PYTHON_INSTALL}"                    )
	string ( CONCAT PYTHON_INSTALL ${PYTHON_INSTALL} ${PY_VERSION}                        )
	string ( REGEX REPLACE "\n$" "" PYTHON_INSTALL "${PYTHON_INSTALL}"                    )
	string ( CONCAT PYTHON_INSTALL2 ${PYTHON_INSTALL} "/site-packages"                    )
	string ( REGEX REPLACE "\n$" "" PYTHON_INSTALL2 "${PYTHON_INSTALL2}"                  )
	if    ( EXISTS ${PYTHON_INSTALL2} )
		set ( MY_PY2_INSTALL_LOCATION ${PYTHON_INSTALL2} CACHE PATH "Python2 install location" )
	else  ( EXISTS ${PYTHON_INSTALL2} )
		set ( MY_PY2_INSTALL_LOCATION ${PYTHON_INSTALL} CACHE PATH "Python2 install location" )
	endif ( EXISTS ${PYTHON_INSTALL2} )
	
	execute_process      ( COMMAND python3 -c "import sys; print ( str ( sys.version_info[0] ) + '.' + str ( sys.version_info[1] ) )" OUTPUT_VARIABLE PY3_VERSION )
	string ( REGEX REPLACE "\n$" "" PY3_VERSION "${PY3_VERSION}"                          )
	execute_process      ( COMMAND python3 -c "import sys; print ( sys.exec_prefix )" OUTPUT_VARIABLE PY3_PREFIX )
	string ( REGEX REPLACE "\n$" "" PY3_PREFIX "${PY3_PREFIX}"                            )
	string ( CONCAT PYTHON3_INSTALL ${PY3_PREFIX} "/lib/python"                           )
	string ( REGEX REPLACE "\n$" "" PYTHON3_INSTALL "${PYTHON3_INSTALL}"                  )
	string ( CONCAT PYTHON3_INSTALL ${PYTHON3_INSTALL} ${PY3_VERSION}                     )
	string ( REGEX REPLACE "\n$" "" PYTHON3_INSTALL "${PYTHON3_INSTALL}"                  )
	string ( CONCAT PYTHON3_INSTALL2 ${PYTHON3_INSTALL} "/site-packages"                  )
	string ( REGEX REPLACE "\n$" "" PYTHON3_INSTALL2 "${PYTHON3_INSTALL2}"                )
	if    ( EXISTS ${PYTHON3_INSTALL2} )
		set ( MY_PY3_INSTALL_LOCATION ${PYTHON3_INSTALL2} CACHE PATH "Python3 install location" )
	else  ( EXISTS ${PYTHON3_INSTALL2} )
		set ( MY_PY3_INSTALL_LOCATION ${PYTHON3_INSTALL} CACHE PATH "Python3 install location" )
	endif ( EXISTS ${PYTHON3_INSTALL2} )
	
endif ( INSTALL_LOCALLY )
message                 ( STATUS "Installation is set to       : " ${MY_INSTALL_LOCATION}         )
message                 ( STATUS "Python 2 Installation        : " ${MY_PY2_INSTALL_LOCATION} )
message                 ( STATUS "Python 3 Installation        : " ${MY_PY3_INSTALL_LOCATION} )

##########################################################################################
################################### Build the OBJECTS
if    ( ${USE_MODERN_PATHS} MATCHES 1 )
	add_subdirectory ( ${CMAKE_SOURCE_DIR}/cmake/sources                                  )
	set              ( OBJS $<TARGET_OBJECTS:compileObjs> CACHE INTERNAL "ProSHADE objects" )
endif ( ${USE_MODERN_PATHS} MATCHES 1 )

##########################################################################################
################################### Link executable
if     ( ${USE_MODERN_PATHS} MATCHES 1 )
	add_subdirectory ( ${CMAKE_SOURCE_DIR}/cmake/bin                                      )
elseif ( ${USE_MODERN_PATHS} MATCHES 0 ) 
	include          ( ${CMAKE_SOURCE_DIR}/cmake/bin/bin.cmake                            )
endif  ( ${USE_MODERN_PATHS} MATCHES 1 )

##########################################################################################
################################### Link library
add_subdirectory      ( ${CMAKE_SOURCE_DIR}/cmake/lib                                     )

##########################################################################################
################################### Link python module
#add_subdirectory      ( ${CMAKE_SOURCE_DIR}/cmake/python                                  )

##########################################################################################
################################### Link python3 module
add_subdirectory      ( ${CMAKE_SOURCE_DIR}/cmake/python3                                 )

#####################################################################################
################################### Remove all the ugly files 
set                ( CLEAR_VARIABLE "${CMAKE_SOURCE_DIR}/clearme.txt"                     )
add_custom_command ( 
					 OUTPUT ${CLEAR_VARIABLE}
					 COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --red "Removing all files except for the targets. New cmake run will be required for re-compilation, but all targets remain working."
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/cmake
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/cmake_install.cmake
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/CMakeCache.txt
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/CMakeFiles
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/Makefile
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/python3 
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/install_manifest.txt
				     DEPENDS "${CMAKE_SOURCE_DIR}/CMakeLists.txt"
  					 COMMENT "Removing all CMake created files but keeping the targets untouched."
 					 VERBATIM
                   )
add_custom_target  ( clear DEPENDS ${CLEAR_VARIABLE}                                      )

#####################################################################################
################################### Remove all non-essential files including compiled files 
set                ( DELETE_VARIABLE "${CMAKE_SOURCE_DIR}/delme.txt"                      )
add_custom_command ( 
					 OUTPUT ${DELETE_VARIABLE}
					 COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --red "Removing all files including the targets. New cmake - make run will be required for the targets to be working again."
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/cmake
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/cmake_install.cmake
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/CMakeCache.txt
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/CMakeFiles
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/Makefile
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/python3 
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/install_manifest.txt
					 COMMAND rm -fr ${MY_INSTALL_LOCATION}/bin/proshade
					 COMMAND rm -fr ${MY_INSTALL_LOCATION}/include/ProSHADE.hpp
					 COMMAND rm -fr ${MY_INSTALL_LOCATION}/lib/libproshade.dylib
					 COMMAND rm -fr ${MY_INSTALL_LOCATION}/python3/__pycache__
					 COMMAND rm -fr ${MY_INSTALL_LOCATION}/python3/_proshade.so
					 COMMAND rm -fr ${MY_INSTALL_LOCATION}/python3/proshade.py
				     DEPENDS "${CMAKE_SOURCE_DIR}/CMakeLists.txt"
  					 COMMENT "Removing all CMake created files but keeping the targets untouched."
 					 VERBATIM
                   )
add_custom_target  ( remove DEPENDS ${DELETE_VARIABLE}                                    )
