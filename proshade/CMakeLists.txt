##########################################################################################
##########################################################################################
################################### ProSHADE cmake file ##################################
##########################################################################################
##########################################################################################

##########################################################################################
################################### Decide version of CMake and use appropriate fn calls
set     ( PRODUCE_PYTHON_BINDINGS   "0"                                                   )
set     ( USE_MODERN_PATHS          "0"                                                   )

if      ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 6 )
	message ( FATAL_ERROR "Your CMake version is less than 2.6. ProSHADE requires a newer version and cannot be installed using the current version. Please update your CMake and try again." )
endif   ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 6 )

if      ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 5 )
	if    ( ${CMAKE_MINOR_VERSION} GREATER 7 )
		message ( "Your CMake version is lower than recommended (3.4). This will result in no Python language bindings. Please consider updating your CMake." )
		set ( PRODUCE_PYTHON_BINDINGS   "0"                                               )
		set ( USE_MODERN_PATHS          "1"                                               )
	else  ( ${CMAKE_MINOR_VERSION} GREATER 7 )
		message ( "Your CMake version is lower than recommended (3.4). This will result in no Python language bindings and slower compilation. Please consider updating your CMake." )
		set ( PRODUCE_PYTHON_BINDINGS   "0"                                               )
		set ( USE_MODERN_PATHS          "0"                                               )
	endif ( ${CMAKE_MINOR_VERSION} GREATER 7 )
endif   ( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 5 )

if      ( ${CMAKE_MAJOR_VERSION} GREATER 2 )
	if    ( ${CMAKE_MINOR_VERSION} LESS 4 )
		message ( "Your CMake version is lower than recommended (3.4). This will result in no Python language bindings. Please consider updating your CMake." )
		set ( PRODUCE_PYTHON_BINDINGS   "0"                                               )
		set ( USE_MODERN_PATHS          "1"                                               )
	else  ( ${CMAKE_MINOR_VERSION} LESS 4 )
		message ( "Your CMake version is at least the recommended version. Allowing fast build and Python language binding to be created." )
		set ( PRODUCE_PYTHON_BINDINGS   "1"                                               )
		set ( USE_MODERN_PATHS          "1"                                               )
	endif ( ${CMAKE_MINOR_VERSION} LESS 4 )
endif   ( ${CMAKE_MAJOR_VERSION} GREATER 2 )

##########################################################################################
################################### Set the correct required version
if      ( ${PRODUCE_PYTHON_BINDINGS} MATCHES 1 )
	cmake_minimum_required  ( VERSION 3.4                                                 )
else    ( ${PRODUCE_PYTHON_BINDINGS} MATCHES 1 )
	cmake_minimum_required  ( VERSION 2.6                                                 )
endif   ( ${PRODUCE_PYTHON_BINDINGS} MATCHES 1 )

##########################################################################################
################################### Set project name
project ( proshade VERSION 0.7.5.5 LANGUAGES CXX C                                        )     

##########################################################################################
################################### Force C++11
set     ( CMAKE_CXX_STANDARD 11                                                           )
set     ( CMAKE_CXX_STANDARD_REQUIRED ON                                                  )

##########################################################################################
################################### Force Position Independent Code (fPIC)
set     ( CMAKE_POSITION_INDEPENDENT_CODE ON                                              )

##########################################################################################
################################### Find the source files
file    ( GLOB SOURCES   "${CMAKE_SOURCE_DIR}/src/proshade/*.cpp"                         )
file    ( GLOB EXEC_SRC  "${CMAKE_SOURCE_DIR}/src/bin/*.cpp"                              )

##########################################################################################
################################### Define the command line options variables and defaults
set     ( INSTALL_LOCALLY        ON      CACHE BOOL   "Should installation be done locally?"   )
set     ( INSTALL_BIN_DIR        bin     CACHE PATH   "Installation directory for executables" )
set     ( INSTALL_LIB_DIR        lib     CACHE PATH   "Installation directory for libraries"   )
set     ( INSTALL_INC_DIR        include CACHE PATH   "Installation directory for headers"     )
set     ( CUSTOM_FFTW3_LIB_PATH  ""      CACHE STRING "Path to libfftw3.a"                     )
set     ( CUSTOM_FFTW3_INC_PATH  ""      CACHE STRING "Path to fftw3.h"                        )
set     ( CUSTOM_LAPACK_LIB_PATH ""      CACHE STRING "Custom path to LAPACK libs"             )
set     ( BUILD_PYTHON           FALSE   CACHE BOOL   "Should python modules be build?"        )
set     ( CMAKE_BUILD_TYPE       Debug                                                         )

##########################################################################################
################################### Warn if non-local installation - sudo may be needed
if    ( NOT INSTALL_LOCALLY )
    message ( WARNING "You have requested installation to system folders. This has two implications:\n 1) You may want to set the -DINSTALL_BIN_DIR=/binary/install/path,\n                        the -DINSTALL_LIB_DIR=/library/install/path\n                    and the -DINSTALL_INC_DIR=/include/install/path\n    to the correct paths on your system, as the defaults are /usr/loca/lib, /usr/local/bin and /usr/local/include.\n 2) You may need to use sudo for 'make install' as well as 'make remove' as this may be required in order to write into the system folders." )
endif ( NOT INSTALL_LOCALLY )

##########################################################################################
################################### Set flags depending on compiler and build type
if         ( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
	if     ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
  		################################# Using Clang
  		set   ( CMAKE_CXX_FLAGS " -O0 -g3 -Wshadow -Wpedantic -Wall -Wextra -Wdouble-promotion -Wformat=2 -Weverything -Wundef -Wconversion -ffunction-sections -fdata-sections -flto -fno-common -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-documentation -Wno-padded -Wno-poison-system-directories -fPIC " )
	elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
  		################################# Using Clang
  		set   ( CMAKE_CXX_FLAGS " -O0 -g3 -Wshadow -Wpedantic -Wall -Wextra -Wdouble-promotion -Wformat=2 -Weverything -Wundef -Wconversion -ffunction-sections -fdata-sections -flto -fno-common -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-documentation -Wno-padded -Wno-poison-system-directories -fPIC " )
	elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"   )
  		################################# Using GCC
  		set   ( CMAKE_CXX_FLAGS " -O0 -g3 -Wshadow -Wpedantic -Wall -Wextra -Wdouble-promotion -Wformat=2 -Wformat-truncation -Wformat-overflow -fstack-usage -Wundef -Wconversion -ffunction-sections -fdata-sections -fno-common -Wno-padded -fPIC"       )
#	elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
#  		################################# Using Intel C++
#  		### TO BE COMPLETED
	elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC"  )
  		################################# Using Visual Studio C++
  		set   ( CMAKE_CXX_FLAGS " /Od /Z7 /DEBUG /W4 "                  )
	endif  ( )   

elseif     ( "${CMAKE_BUILD_TYPE}" STREQUAL "Release" )

	if     ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
  		################################# Using Clang
  		set   ( CMAKE_CXX_FLAGS " -O3 -fPIC "                                             )
  	elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
  		################################# Using Clang
  		set   ( CMAKE_CXX_FLAGS " -O3 -fPIC "                                             ) 	
	elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"   )
  		################################# Using GCC
  		set   ( CMAKE_CXX_FLAGS " -O3 -fPIC "                                             )
#	elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
#  		################################# Using Intel C++
#  		### TO BE COMPLETED
	elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC"  )
  		################################# Using Visual Studio C++
  		set   ( CMAKE_CXX_FLAGS " /O2 /GL /Gw /Qpar /Qspectre /Qspectre-load /Qspectre-load-cf /RTC " )
	endif  ( )   

else       ( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )

	message ( ERROR "Unrecognised build type. Will assume release, but beware that this should not have happened!" )

endif      ( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )  

##########################################################################################
################################### Set flags depending on system
if     ( CYGWIN )
	set              ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__int64=\"long long\""       )
elseif ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  )
	link_directories    ( ${CMAKE_SOURCE_DIR}/winLibs/x64/ZLIB                            )
	include_directories ( ${CMAKE_SOURCE_DIR}/winLibs/x64/ZLIB                            )
	link_directories    ( ${CMAKE_SOURCE_DIR}/winLibs/x64/DLLs                            )
endif  ( CYGWIN )

##########################################################################################
################################### Set targets for dependency installation from supplied sources

### GEMMI
if    ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  )
	add_custom_command (
   		OUTPUT ${PROJECT_SOURCE_DIR}/extern/buildGEMMI
   		COMMAND rmdir "gemmi" /s /q
   		COMMAND git clone https://github.com/project-gemmi/gemmi.git --quiet
   		COMMAND type nul > ${PROJECT_SOURCE_DIR}/extern/buildGEMMI
   		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern
	)
else  ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  )
	add_custom_command (
   		OUTPUT ${PROJECT_SOURCE_DIR}/extern/buildGEMMI
   		COMMAND rm -fr ./gemmi
   		COMMAND git clone https://github.com/project-gemmi/gemmi.git --quiet
   		COMMAND touch ${PROJECT_SOURCE_DIR}/extern/buildGEMMI
   		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern
	)
endif ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  )

add_custom_target (
   gemmi_lib
   DEPENDS ${PROJECT_SOURCE_DIR}/extern/buildGEMMI
)

### GetOpt_port
include_directories ( ${CMAKE_SOURCE_DIR}/extern/                                         ) 

### Set paths to installed dependencies
include_directories ( ${CMAKE_SOURCE_DIR}/extern/gemmi/include                            )  

##########################################################################################
################################### Add fftw3 dependency

### Automatic FFTW3 detection
if    ( CUSTOM_FFTW3_LIB_PATH OR CUSTOM_FFTW3_INC_PATH )
    message         ( STATUS "Custom FFTW3 path was supplied. Checking for the paths containing the library and the include files." )
    
    file            ( GLOB TESTED_FFTW3_LIB_PATH ${CUSTOM_FFTW3_LIB_PATH}/libfftw3*       )
    if    ( TESTED_FFTW3_LIB_PATH )
        list            ( GET TESTED_FFTW3_LIB_PATH 0 TESTED_FFTW3_LIB_PATH_HEAD          )
        get_filename_component ( TESTED_FFTW3_LIB_PATH_HEAD_PARENT ${TESTED_FFTW3_LIB_PATH_HEAD} DIRECTORY )
        link_directories ( ${TESTED_FFTW3_LIB_PATH_HEAD_PARENT}                           )
        set             ( FFTW_LINK ${TESTED_FFTW3_LIB_PATH_HEAD_PARENT} CACHE STRING "RPATH value to FFTW3 lib" )
        set             ( FFTW3_LIB_PATH ${TESTED_FFTW3_LIB_PATH_HEAD_PARENT} CACHE STRING "internal fftw3 library path" )
        message         ( STATUS "FFTW3 library path is        : " ${TESTED_FFTW3_LIB_PATH_HEAD_PARENT} )
    else  ( TESTED_FFTW3_LIB_PATH )
        message             ( FATAL_ERROR "Cannot find FFTW lib in: " ${CUSTOM_FFTW3_LIB_PATH} )
    endif ( TESTED_FFTW3_LIB_PATH )
    
    file            ( GLOB TESTED_FFTW3_INC_PATH ${CUSTOM_FFTW3_INC_PATH}/fftw3.h         )
    if    ( TESTED_FFTW3_INC_PATH )
        list            ( GET TESTED_FFTW3_INC_PATH 0 TESTED_FFTW3_INC_PATH_HEAD          )
        get_filename_component ( TESTED_FFTW3_INC_PATH_HEAD_PARENT ${TESTED_FFTW3_INC_PATH_HEAD} DIRECTORY)
        include_directories ( ${TESTED_FFTW3_INC_PATH_HEAD_PARENT}                        )
        set             ( FFTW3_INC_PATH ${TESTED_FFTW3_INC_PATH_HEAD_PARENT} CACHE STRING "internal fftw3 include path" )
        message         ( STATUS "FFTW3 include path is        : " ${TESTED_FFTW3_INC_PATH_HEAD_PARENT} )
    else  ( TESTED_FFTW3_INC_PATH )
        message             ( FATAL_ERROR "Cannot find FFTW3 inc in: " ${TESTED_FFTW3_INC_PATH} )
    endif ( TESTED_FFTW3_INC_PATH )
elseif( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  )
	link_directories    ( ${CMAKE_SOURCE_DIR}/winLibs/x64/FFTW3                           )
	include_directories ( ${CMAKE_SOURCE_DIR}/winLibs/x64/FFTW3                           )
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
################################### Add LAPACK dependency
 
### Automatic LAPACK detection
if    ( CUSTOM_LAPACK_LIB_PATH )
    message            ( STATUS "Custom LAPACK path was supplied. Checking for the paths containing the library." )
    file               ( GLOB TESTED_LAPACK_PATH ${CUSTOM_LAPACK_LIB_PATH}/liblapack*     )
    if    ( TESTED_LAPACK_PATH )
        list                ( GET TESTED_LAPACK_PATH 0 TESTED_LAPACK_PATH_HEAD            )
        message             ( STATUS "Found LAPACK library in      : " ${TESTED_LAPACK_PATH_HEAD} )
        link_directories    ( ${TESTED_LAPACK_PATH_HEAD}                                  )
        set                 ( LAPACK_LINK ${TESTED_LAPACK_PATH_HEAD} CACHE STRING "RPATH value to LAPACK libs" )
    else  ( TESTED_LAPACK_PATH )
        message             ( FATAL_ERROR "Cannot find LAPACK in   : " ${CUSTOM_LAPACK_LIB_PATH} )
    endif ( TESTED_LAPACK_PATH )
elseif( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  )
	link_directories    ( ${CMAKE_SOURCE_DIR}/winLibs/x64/LAPACK                          )
	include_directories ( ${CMAKE_SOURCE_DIR}/winLibs/x64/LAPACK                          )
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
################################### Add getopt_port dependency
set ( GETOPT_SOURCE  ${CMAKE_SOURCE_DIR}/extern/getopt_port/getopt_port.c        CACHE STRING "Getopt_port source" )
set ( GETOPT_HEADER  ${CMAKE_SOURCE_DIR}/extern/getopt_port/getopt_port.h        CACHE STRING "Getopt_port header" )

##########################################################################################
################################### Add almostEqual dependency
set ( ALMOSTEQUAL_SOURCE  ${CMAKE_SOURCE_DIR}/extern/almostEqual CACHE STRING "almostEqual source" )
include_directories    ( ${ALMOSTEQUAL_SOURCE}                                            )

##########################################################################################
################################### Add SOFT dependency

### Compile Soft-2.0 from supplied folder
if    ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  )
	add_custom_command (
        OUTPUT ${PROJECT_SOURCE_DIR}/extern/buildSOFT2
   		COMMAND type nul > ${PROJECT_SOURCE_DIR}/extern/buildSOFT2
   		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/soft-2.0
   	)
else  ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  ) 
	add_custom_command (
        OUTPUT ${PROJECT_SOURCE_DIR}/extern/buildSOFT2
   		COMMAND make lib FFTWINC=-I${FFTW3_INC_PATH} FFTWLIB1=-L${FFTW3_LIB_PATH} FFTWLIB2=-lfftw3 > ${CMAKE_SOURCE_DIR}/extern/soft-2.0/install.log.txt 2> ${CMAKE_SOURCE_DIR}/extern/soft-2.0/install.err.txt
   		COMMAND touch ${PROJECT_SOURCE_DIR}/extern/buildSOFT2
   		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern/soft-2.0
   	)
   	link_directories    ( ${CMAKE_SOURCE_DIR}/winLibs/x64/SOFT2                           )
endif ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  )

add_custom_target (
   soft2_lib
   DEPENDS ${PROJECT_SOURCE_DIR}/extern/buildSOFT2
   )

### Report SOFT related info
message                ( STATUS "SOFT library path is         : " ${CMAKE_SOURCE_DIR}/extern/soft-2.0 )
message                ( STATUS "SOFT include path is         : " ${CMAKE_SOURCE_DIR}/extern/soft-2.0/include )
include_directories    ( ${CMAKE_SOURCE_DIR}/extern/soft-2.0/include                      )

##########################################################################################
################################### Set installation folders
if    ( INSTALL_LOCALLY )
	set ( MY_INSTALL_LOCATION  ${CMAKE_SOURCE_DIR}/install   CACHE STRING "Installation path" )
else  ( INSTALL_LOCALLY )
	set ( MY_INSTALL_LOCATION  ${CMAKE_INSTALL_PREFIX}       CACHE STRING "Installation path" )
endif ( INSTALL_LOCALLY )

message              ( STATUS "Installation is set to       : " ${MY_INSTALL_LOCATION}    )

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
if     ( ${USE_MODERN_PATHS} MATCHES 1 )
	add_subdirectory ( ${CMAKE_SOURCE_DIR}/cmake/lib                                      )
elseif ( ${USE_MODERN_PATHS} MATCHES 0 ) 
	include          ( ${CMAKE_SOURCE_DIR}/cmake/lib/lib.cmake                            )
endif  ( ${USE_MODERN_PATHS} MATCHES 1 )

##########################################################################################
################################### Build python module using PyBind11
if     ( BUILD_PYTHON )
	if     ( ${PRODUCE_PYTHON_BINDINGS} MATCHES 1 )
		add_subdirectory   ( ${CMAKE_SOURCE_DIR}/cmake/python                             )
	else   ( ${PRODUCE_PYTHON_BINDINGS} MATCHES 1 )
		message            ( ERROR "This CMake version does not allow building the python module - please consider updating your CMake to at least version 3.4 if you want the python module build." )
	endif  ( ${PRODUCE_PYTHON_BINDINGS} MATCHES 1 )
else   ( BUILD_PYTHON )
	message                ( STATUS "Will not attempt to build the python module."        )
endif  ( BUILD_PYTHON )

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
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/python
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/python3 
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/install_manifest.txt
				     DEPENDS "${CMAKE_SOURCE_DIR}/CMakeLists.txt"
  					 COMMENT "Removed all CMake created files but keeping the targets untouched."
 					 VERBATIM
                   )
add_custom_target  ( clear DEPENDS ${CLEAR_VARIABLE}                                      )

#####################################################################################
################################### Remove all files including compiled files and targets
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
					 COMMAND rm -fr ${CMAKE_BINARY_DIR}/install_manifest.txt
					 COMMAND rm -fr ${MY_INSTALL_LOCATION}/bin/proshade
					 COMMAND rm -fr ${MY_INSTALL_LOCATION}/bin/proshade.dSYM
					 COMMAND rm -fr ${MY_INSTALL_LOCATION}/include/ProSHADE.hpp
					 COMMAND rm -fr ${MY_INSTALL_LOCATION}/lib/libproshade.dylib
					 COMMAND rm -fr ${MY_INSTALL_LOCATION}/pythonModule/proshade.*
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
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/install.log.txt
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/soft-2.0/install.err.txt
					 COMMAND rm -fr ${PROJECT_SOURCE_DIR}/extern/buildSOFT2
					 COMMAND rm -fr ${PROJECT_SOURCE_DIR}/extern/buildGEMMI
					 COMMAND rm -fr ${CMAKE_SOURCE_DIR}/extern/gemmi
					 DEPENDS "${CMAKE_SOURCE_DIR}/CMakeLists.txt"
				     DEPENDS "${CMAKE_SOURCE_DIR}/delmeSOFT2.txt"
  					 COMMENT "Removed all files including the targets."
 					 VERBATIM
                   )
add_custom_target  ( remove DEPENDS ${DELETE_VARIABLE}                                    )
