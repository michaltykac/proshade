##########################################################################################
##########################################################################################
################################### ProSHADE cmake file ##################################
##########################################################################################
##########################################################################################

##########################################################################################
################################### Decide version of CMake and use appropriate fn calls
set ( PRODUCE_PYTHON_BINDINGS   "0"                                                       )
set ( USE_MODERN_PATHS          "0"                                                       )

if    ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 6 )
	message ( FATAL_ERROR "Your CMake version is less than 2.6. ProSHADE requires a newer version and cannot be installed using the current version. Please update your CMake and try again." )
endif ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 6 )

if    ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 5 )
	if    ( ${CMAKE_MINOR_VERSION} GREATER 7 )
		message ( "Your CMake version is lower than recommended (3.1). This will result in no Python language bindings and slower compilation. Please consider updating your CMake." )
		set ( PRODUCE_PYTHON_BINDINGS   "0"                                               )
		set ( USE_MODERN_PATHS          "1"                                               )
	else  ( ${CMAKE_MINOR_VERSION} GREATER 7 )
		message ( "Your CMake version is lower than recommended (3.1). This will result in no Python language bindings and slower compilation. Please consider updating your CMake." )
		set ( PRODUCE_PYTHON_BINDINGS   "0"                                               )
		set ( USE_MODERN_PATHS          "0"                                               )
	endif ( ${CMAKE_MINOR_VERSION} GREATER 7 )
endif ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 5 )

if    ( ${CMAKE_MAJOR_VERSION} GREATER 2 )
	if    ( ${CMAKE_MINOR_VERSION} LESS 1 )
		message ( "Your CMake version is lower than recommended (3.1). This will result in no Python language bindings and slower compilation. Please consider updating your CMake." )
		set ( PRODUCE_PYTHON_BINDINGS   "0"                                               )
		set ( USE_MODERN_PATHS          "1"                                               )
	else  ( ${CMAKE_MINOR_VERSION} LESS 1 )
		message ( "Your CMake version is at least the recommended version. Allowing fast build and Python language binding to be created." )
		set ( PRODUCE_PYTHON_BINDINGS   "1"                                               )
		set ( USE_MODERN_PATHS          "1"                                               )
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
project ( proshade VERSION 0.7.1 LANGUAGES CXX C                                          )     

##########################################################################################
################################### Disable C++11
set ( CMAKE_CXX_STANDARD 98                                                               )
set ( CMAKE_CXX_STANDARD_REQUIRED ON                                                      )

##########################################################################################
################################### Find the source files
file    ( GLOB SOURCES  "${CMAKE_SOURCE_DIR}/src/proshade/*.cpp"                          )
file    ( GLOB EXEC_SRC "${CMAKE_SOURCE_DIR}/src/bin/*.cpp"                               )

##########################################################################################
################################### Set flags depending on compiler

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  ################################# Using Clang
  set ( CMAKE_CXX_FLAGS "-O3 -g -Wshadow -Wpedantic -Wall -Wno-macro-redefined -fPIC"     )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  ################################# Using GCC
  set ( CMAKE_CXX_FLAGS "-O3 -g -Wshadow -Wpedantic -Wall -Wno-c++11-compat -fPIC"        )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  ################################# Using Intel C++
  ### TO BE COMPLETED
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  ################################# Using Visual Studio C++
  ### TO BE COMPLETED
endif()

##########################################################################################
################################### Set targets for clipper installation from supplied sources

### MMDB2
add_custom_command (
   OUTPUT ${PROJECT_SOURCE_DIR}/extern/install/buildMMDB2
   COMMAND ./configure --prefix=${CMAKE_SOURCE_DIR}/extern/install --includedir=${CMAKE_SOURCE_DIR}/extern/install/include --libdir=${CMAKE_SOURCE_DIR}/extern/install/lib --disable-Werror --enable-silent-rules --enable-shared --disable-static > ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.err.txt
   COMMAND make >> ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.err.txt
   COMMAND make install >> ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.err.txt
   COMMAND touch ${PROJECT_SOURCE_DIR}/extern/install/buildMMDB2
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/mmdb2
)

add_custom_target (
   mmdb2_lib
   DEPENDS ${PROJECT_SOURCE_DIR}/extern/install/buildMMDB2
)

### CCP4LIB
add_custom_command (
   OUTPUT ${PROJECT_SOURCE_DIR}/extern/install/buildCCP4LIB
   DEPENDS mmdb2_lib
   COMMAND env LDFLAGS='-L${CMAKE_SOURCE_DIR}/extern/install/lib' CPPFLAGS='-I${CMAKE_SOURCE_DIR}/extern/install/include' ./configure --prefix=${CMAKE_SOURCE_DIR}/extern/install --includedir=${CMAKE_SOURCE_DIR}/extern/install/include --libdir=${CMAKE_SOURCE_DIR}/extern/install/lib --enable-silent-rules --enable-shared --disable-static  > ${CMAKE_SOURCE_DIR}/extern/libccp4/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/libccp4/install.err.txt
   COMMAND make >> ${CMAKE_SOURCE_DIR}/extern/libccp4/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/libccp4/install.err.txt
   COMMAND make install >> ${CMAKE_SOURCE_DIR}/extern/libccp4/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/libccp4/install.err.txt
   COMMAND touch ${PROJECT_SOURCE_DIR}/extern/install/buildCCP4LIB
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/libccp4   
)

add_custom_target (
   ccp4_lib
   DEPENDS ${PROJECT_SOURCE_DIR}/extern/install/buildCCP4LIB
)

### FFTW2
if      ( APPLE )
    add_custom_command (
       OUTPUT ${PROJECT_SOURCE_DIR}/extern/install/buildFFTW2
       COMMAND ./configure --prefix=${CMAKE_SOURCE_DIR}/extern/install --disable-Werror --enable-float --enable-silent-rules --enable-shared --disable-static --includedir=${CMAKE_SOURCE_DIR}/extern/install/include --libdir=${CMAKE_SOURCE_DIR}/extern/install/lib > ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.log.txt 2> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.err.txt
       COMMAND sed -i '' -e 's/SUBDIRS = fftw rfftw tests doc threads mpi/SUBDIRS = fftw rfftw tests threads mpi/' Makefile  
       COMMAND make >> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.err.txt
       COMMAND make install >> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.err.txt
       COMMAND touch ${PROJECT_SOURCE_DIR}/extern/install/buildFFTW2 
       WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5
    )
else    ( APPLE )
    add_custom_command (
       OUTPUT ${PROJECT_SOURCE_DIR}/extern/install/buildFFTW2
       COMMAND ./configure --prefix=${CMAKE_SOURCE_DIR}/extern/install --disable-Werror --enable-float --enable-silent-rules --enable-shared --disable-static --includedir=${CMAKE_SOURCE_DIR}/extern/install/include --libdir=${CMAKE_SOURCE_DIR}/extern/install/lib > ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.log.txt 2> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.err.txt
       COMMAND sed -i 's/SUBDIRS = fftw rfftw tests doc threads mpi/SUBDIRS = fftw rfftw tests threads mpi/' Makefile 
       COMMAND make >> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.err.txt
       COMMAND make install >> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.err.txt
       COMMAND touch ${PROJECT_SOURCE_DIR}/extern/install/buildFFTW2 
       WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5 )
endif   ( APPLE )

add_custom_target (
   fftw2_lib
   DEPENDS ${PROJECT_SOURCE_DIR}/extern/install/buildFFTW2
   )

### CLIPPER
add_custom_command (
   OUTPUT ${PROJECT_SOURCE_DIR}/extern/install/buildCLIPPER
   DEPENDS fftw2_lib
   DEPENDS ccp4_lib
   COMMAND env LDFLAGS='-L${CMAKE_SOURCE_DIR}/extern/install/lib' CPPFLAGS='-I${CMAKE_SOURCE_DIR}/extern/install/include' ./configure --prefix=${CMAKE_SOURCE_DIR}/extern/install -enable-shared --includedir=${CMAKE_SOURCE_DIR}/extern/install/include --libdir=${CMAKE_SOURCE_DIR}/extern/install/lib --with-mmdb=${CMAKE_SOURCE_DIR}/extern/install --with-fftw=${CMAKE_SOURCE_DIR}/extern/install --enable-mmdb --enable-minimol --enable-cif --enable-cns --enable-ccp4 --with-ccp4=${CMAKE_SOURCE_DIR}/extern/install > ${CMAKE_SOURCE_DIR}/extern/clipper/install.log.txt 2> ${CMAKE_SOURCE_DIR}/extern/clipper/install.err.txt 
   COMMAND make >> ${CMAKE_SOURCE_DIR}/extern/clipper/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/clipper/install.err.txt
   COMMAND make install >> ${CMAKE_SOURCE_DIR}/extern/clipper/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/clipper/install.err.txt
   COMMAND touch ${PROJECT_SOURCE_DIR}/extern/install/buildCLIPPER
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/clipper
)

add_custom_target (
   clipper_lib
   DEPENDS ${PROJECT_SOURCE_DIR}/extern/install/buildCLIPPER
)

set               ( CLIPPER_LINK ${CMAKE_SOURCE_DIR}/extern/install/lib CACHE STRING "RPATH value to MMDB2/CMAPLIB/CLIPPER libs" )

### Set command line accessible variable for clipper path
include_directories ( ${CMAKE_SOURCE_DIR}/extern/clipper                                  )
include_directories ( ${CMAKE_SOURCE_DIR}/extern/install/include                          )
link_directories    ( ${CMAKE_SOURCE_DIR}/extern/install/lib                              )            

##########################################################################################
################################### Add fftw3 dependency

### Automatic FFTW3 detection
set ( CUSTOM_FFTW3_LIB_PATH ""     CACHE STRING "path to libfftw3.a"                      )
set ( CUSTOM_FFTW3_INC_PATH ""     CACHE STRING "path to fftw3.h"                         )

if    ( CUSTOM_FFTW3_LIB_PATH OR CUSTOM_FFTW3_INC_PATH )
    message         ( STATUS "Custom FFTW3 path was supplied. Checking for the paths containing the library and the include files." )
    
    file            ( GLOB TESTED_FFTW3_LIB_PATH ${CUSTOM_LAPACK_LIB_PATH}/libfftw3*      )
    if    ( TESTED_FFTW3_LIB_PATH )
        list            ( GET TESTED_FFTW3_LIB_PATH 0 TESTED_FFTW3_LIB_PATH_HEAD          )
        link_directories ( ${TESTED_FFTW3_LIB_PATH_HEAD}                                  )
        set             ( FFTW_LINK ${TESTED_FFTW3_LIB_PATH_HEAD} CACHE STRING "RPATH value to FFTW3 lib" )
        set             ( FFTW3_LIB_PATH ${TESTED_FFTW3_LIB_PATH_HEAD} CACHE STRING "internal fftw3 library path" )
        message         ( STATUS "FFTW3 library path is        : " ${FFTW3_LIB_PATH}      )
    else  ( TESTED_LAPACK_PATH )
        message             ( FATAL_ERROR "Cannot find LAPACK lib in: " ${CUSTOM_LAPACK_LIB_PATH} )
    endif ( TESTED_LAPACK_PATH )
    
    file            ( GLOB TESTED_FFTW3_INC_PATH ${CUSTOM_FFTW3_INC_PATH}/fftw3.h         )
    if    ( TESTED_FFTW3_INC_PATH )
        list            ( GET TESTED_FFTW3_INC_PATH 0 TESTED_FFTW3_INC_PATH_HEAD          )
        include_directories ( ${TESTED_FFTW3_INC_PATH_HEAD}                               )
        set             ( FFTW3_INC_PATH ${TESTED_FFTW3_INC_PATH_HEAD} CACHE STRING "internal fftw3 include path" )
        message         ( STATUS "FFTW3 include path is        : " ${FFTW3_INC_PATH}      )
    else  ( TESTED_FFTW3_INC_PATH )
        message             ( FATAL_ERROR "Cannot find LAPACK inc in: " ${TESTED_FFTW3_INC_PATH} )
    endif ( TESTED_FFTW3_INC_PATH )
else  ( CUSTOM_FFTW3_LIB_PATH OR CUSTOM_FFTW3_INC_PATH )
    set                 ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/findFFTW            )
    find_package        ( FFTW REQUIRED                                                   )
    if    ( FFTW_FOUND )    
        link_directories ( ${FFTW_LIBRARIES}                                              )
        set             ( FFTW_LINK ${FFTW_LIBRARIES} CACHE STRING "RPATH value to FFTW3 lib" )
        include_directories ( ${FFTW_INCLUDES}                                            )
        set             ( FFTW3_LIB_PATH ${FFTW_LIBRARIES} CACHE STRING "internal fftw3 library path" )
        set             ( FFTW3_INC_PATH ${FFTW_INCLUDES} CACHE STRING "internal fftw3 include path" )
        message         ( STATUS "FFTW3 library path is        : " ${FFTW3_LIB_PATH}      )
        message         ( STATUS "FFTW3 include path is        : " ${FFTW3_INC_PATH}      )
    else  ( FFTW_FOUND )    
        message         ( FATAL_ERROR "Could not find the FFTW3 dependency. Please install FFTW3 to a standard location or use the cmake -DCUSTOM_FFTW3_LIB_PATH=/path/to/libfftw3.a/so/dylib and the -DCUSTOM_FFTW3_INC_PATH=/path/to/fftw3.h command line arguments to supply the paths to the locations of the libfftw3.a/so/dylib library and its header file. Alternatively, contact the author for help." )
    endif ( FFTW_FOUND ) 
endif  ( CUSTOM_FFTW3_LIB_PATH OR CUSTOM_FFTW3_INC_PATH )

##########################################################################################
################################### Add SOFT dependency

### Compile Soft-2.0 from supplied folder
add_custom_command (
   OUTPUT ${PROJECT_SOURCE_DIR}/extern/install/buildSOFT2
   COMMAND make lib FFTWINC=-I${FFTW3_INC_PATH} FFTWLIB1=-L${FFTW3_LIB_PATH} FFTWLIB2=-lfftw3 > ${CMAKE_SOURCE_DIR}/extern/soft-2.0/install.log.txt 2> ${CMAKE_SOURCE_DIR}/extern/soft-2.0/install.err.txt
   COMMAND touch ${PROJECT_SOURCE_DIR}/extern/install/buildSOFT2
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/soft-2.0
)

add_custom_target (
   soft2_lib
   DEPENDS ${PROJECT_SOURCE_DIR}/extern/install/buildSOFT2
   )

### Report SOFT related info
message                ( STATUS "SOFT library path is         : " ${CMAKE_SOURCE_DIR}/extern/soft-2.0 )
message                ( STATUS "SOFT include path is         : " ${CMAKE_SOURCE_DIR}/extern/soft-2.0/include )
include_directories    ( ${CMAKE_SOURCE_DIR}/extern/soft-2.0/include                      )

##########################################################################################
################################### Add LAPACK dependency
 
### Automatic LAPACK detection
set                    ( CUSTOM_LAPACK_LIB_PATH "" CACHE STRING "Custom path to LAPACK libs" )

if    ( CUSTOM_LAPACK_LIB_PATH )
    message            ( STATUS "Custom LAPACK path was supplied. Checking for the paths containing the library." )
    file               ( GLOB TESTED_LAPACK_PATH ${CUSTOM_LAPACK_LIB_PATH}/liblapack* )
    if    ( TESTED_LAPACK_PATH )
        list                ( GET TESTED_LAPACK_PATH 0 TESTED_LAPACK_PATH_HEAD            )
        message             ( STATUS "Found LAPACK library in      : " ${TESTED_LAPACK_PATH_HEAD} )
        link_directories    ( ${TESTED_LAPACK_PATH_HEAD}                                  )
        set                 ( LAPACK_LINK ${TESTED_LAPACK_PATH_HEAD} CACHE STRING "RPATH value to LAPACK libs" )
    else  ( TESTED_LAPACK_PATH )
        message             ( FATAL_ERROR "Cannot find LAPACK in   : " ${CUSTOM_LAPACK_LIB_PATH} )
    endif ( TESTED_LAPACK_PATH )
else  ( CUSTOM_LAPACK_LIB_PATH )
    set                 ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/findLapack          )
    find_package        ( lapack REQUIRED                                                 )
    if    ( LAPACK_FOUND )    
        list            ( GET LAPACK_LIBRARIES 0 LAPACK_LIBRARIES_HEAD                    )
        get_filename_component  (LAPACK_LIBRARIES_FOLDER ${LAPACK_LIBRARIES_HEAD} DIRECTORY )
        link_directories ( ${LAPACK_LIBRARIES_FOLDER}                                     )
        message         ( STATUS "Found LAPACK library in      : " ${LAPACK_LIBRARIES_FOLDER} )
        set             ( LAPACK_LINK ${LAPACK_LIBRARIES_FOLDER} CACHE STRING "RPATH value to LAPACK libs" )
    else  ( LAPACK_FOUND )
        message         ( FATAL_ERROR "Could not find the LAPACK dependency. Please install LAPACK to a standard location or use the cmake -DCUSTOM_LAPACK_LIB_PATH=/path/to/liblapack command line argument to supply the path to the location of the liblapack.a/so/dylib library. Alternatively, contact the author for help." )
    endif ( LAPACK_FOUND )
endif  ( CUSTOM_LAPACK_LIB_PATH )

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
  					 COMMENT "Removed all CMake created files but keeping the targets untouched."
 					 VERBATIM
                   )
add_custom_target  ( clear DEPENDS ${CLEAR_VARIABLE}                                      )

#####################################################################################
################################### Remove all files including compiled files and targets
set                ( DELETEMMDB2_VARIABLE "${CMAKE_SOURCE_DIR}/delmeMMDB2.txt"            )
add_custom_command ( 
					 OUTPUT ${DELETEMMDB2_VARIABLE}
					 COMMAND make clean >> ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.err.txt
					 COMMAND make uninstall >> ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.err.txt
					 DEPENDS ${PROJECT_SOURCE_DIR}/CMakeLists.txt
					 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/mmdb2
				   )
				   
set                ( DELETECCP4_VARIABLE "${CMAKE_SOURCE_DIR}/delmeLIBCCP4.txt"           )
add_custom_command ( 
					 OUTPUT ${DELETECCP4_VARIABLE}
					 COMMAND make clean >> ${CMAKE_SOURCE_DIR}/extern/libccp4/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/libccp4/install.err.txt
					 COMMAND make uninstall >> ${CMAKE_SOURCE_DIR}/extern/libccp4/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/libccp4/install.err.txt
					 DEPENDS ${PROJECT_SOURCE_DIR}/CMakeLists.txt
					 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/libccp4
				   )
				   
set                ( DELETEFFTW2_VARIABLE "${CMAKE_SOURCE_DIR}/delmeFFTW2.txt"            )
add_custom_command ( 
					 OUTPUT ${DELETEFFTW2_VARIABLE}
					 COMMAND make clean >> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.err.txt
					 COMMAND make uninstall >> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.err.txt
					 DEPENDS ${PROJECT_SOURCE_DIR}/CMakeLists.txt
					 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5
				   )
				   
set                ( DELETECLIPPER_VARIABLE "${CMAKE_SOURCE_DIR}/delmeCLIPPER.txt"        )
add_custom_command ( 
					 OUTPUT ${DELETECLIPPER_VARIABLE}
					 COMMAND make clean >> ${CMAKE_SOURCE_DIR}/extern/clipper/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/clipper/install.err.txt
					 COMMAND make uninstall >> ${CMAKE_SOURCE_DIR}/extern/clipper/install.log.txt 2>> ${CMAKE_SOURCE_DIR}/extern/clipper/install.err.txt
					 DEPENDS ${PROJECT_SOURCE_DIR}/CMakeLists.txt
					 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/clipper
				   )
				   
set                ( DELETESOFT2_VARIABLE "${CMAKE_SOURCE_DIR}/delmeSOFT2.txt"            )
add_custom_command ( 
					 OUTPUT ${DELETESOFT2_VARIABLE}
					 COMMAND make clean > ${CMAKE_SOURCE_DIR}/extern/soft-2.0/install.log.txt 2> ${CMAKE_SOURCE_DIR}/extern/soft-2.0/install.err.txt
					 DEPENDS ${PROJECT_SOURCE_DIR}/CMakeLists.txt
					 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/soft-2.0
				   )				   

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
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/rotate_so3_fftw.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/rotate_so3_utils.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/s2_cospmls.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/s2_legendreTransforms.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/s2_semi_fly.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/s2_semi_memo.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/so3_correlate_fftw.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/soft_fftw_nt.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/soft_fftw_pc.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/soft_fftw.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/utils_vec_cx.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/wignerTransforms_fftw.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/wrap_s2_rotate_fftw.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/wrap_soft_fftw_cor2.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/lib1/wrap_soft_fftw.o
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/libsoft1.a
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/clipper/install.log.txt
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/clipper/install.err.txt
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.log.txt
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/fftw-2.1.5/install.err.txt
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/libccp4/install.log.txt
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/libccp4/install.err.txt
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.log.txt
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/mmdb2/install.err.txt
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/install.log.txt
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/install.err.txt
					 COMMAND rm -fr ${PROJECT_SOURCE_DIR}/extern/install/buildMMDB2
					 COMMAND rm -fr ${PROJECT_SOURCE_DIR}/extern/install/buildCCP4LIB
					 COMMAND rm -fr ${PROJECT_SOURCE_DIR}/extern/install/buildFFTW2
					 COMMAND rm -fr ${PROJECT_SOURCE_DIR}/extern/install/buildCLIPPER
					 COMMAND rm -fr ${PROJECT_SOURCE_DIR}/extern/install/buildSOFT2
					 COMMAND rm -fr ${PROJECT_SOURCE_DIR}/extern/install/include
					 COMMAND rm -fr ${PROJECT_SOURCE_DIR}/extern/install/lib
					 COMMAND rm -fr ${PROJECT_SOURCE_DIR}/extern/install/lib64
					 COMMAND rm -fr ${PROJECT_SOURCE_DIR}/extern/install/share
					 DEPENDS "${CMAKE_SOURCE_DIR}/delmeMMDB2.txt"
					 DEPENDS "${CMAKE_SOURCE_DIR}/delmeLIBCCP4.txt"
					 DEPENDS "${CMAKE_SOURCE_DIR}/delmeFFTW2.txt"
				     DEPENDS "${CMAKE_SOURCE_DIR}/CMakeLists.txt"
				     DEPENDS "${CMAKE_SOURCE_DIR}/delmeCLIPPER.txt"
				     DEPENDS "${CMAKE_SOURCE_DIR}/delmeSOFT2.txt"
  					 COMMENT "Removed all files including the targets."
 					 VERBATIM
                   )
add_custom_target  ( remove DEPENDS ${DELETE_VARIABLE}                                    )
