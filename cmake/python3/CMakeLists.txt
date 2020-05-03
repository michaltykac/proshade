##########################################################################################
##########################################################################################
############################### ProSHADE python3 cmake file ##############################
##########################################################################################
##########################################################################################

##########################################################################################
################################### Find the SWIG package 
find_package         ( SWIG REQUIRED                                                      )

##########################################################################################
################################### Set swig module name with the new naming policy
if    ( ${CMAKE_VERSION} GREATER "3.12" )
    cmake_policy     ( SET CMP0078 NEW                                                    )
else  ( ${CMAKE_VERSION} GREATER "3.12" )
    set              ( UseSWIG_TARGET_NAME_PREFERENCE "${PROJECT_NAME}_py3Module"         )
endif ( ${CMAKE_VERSION} GREATER "3.12" )

if    ( ${CMAKE_VERSION} GREATER "3.13" )
    cmake_policy     ( SET CMP0086 NEW                                                    )
endif ( ${CMAKE_VERSION} GREATER "3.13" )

##########################################################################################
################################### Set swig module name with the new naming policy
set ( SWIG_MODULE_NAME ${PROJECT_NAME}_py3Module                                          )

##########################################################################################
################################### Include the swig locations
include              ( ${SWIG_USE_FILE}                                                   )

##########################################################################################
################################### Copy the files to be pythonised into the building folder
execute_process      ( COMMAND bash -c "mkdir python3" )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/src/proshade/ProSHADE_typedefs.hpp ${CMAKE_BINARY_DIR}/python3/ProSHADE_typedefs.hpp" )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/src/proshade/ProSHADE_settings.hpp ${CMAKE_BINARY_DIR}/python3/ProSHADE_settings.hpp" )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/src/proshade/ProSHADE.hpp          ${CMAKE_BINARY_DIR}/python3/ProSHADE.hpp"          )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/src/python/proshade.i              ${CMAKE_BINARY_DIR}/python3/proshade.i"            )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/src/python/numpy.i                 ${CMAKE_BINARY_DIR}/python3/numpy.i"               )

##########################################################################################
################################### Find the python include locations
execute_process      ( COMMAND python3 -c "import sysconfig as sysconfig; import sys; sys.stdout.write (sysconfig.get_path(scheme='posix_prefix',name='include'))" OUTPUT_VARIABLE REAL_PYTHON_INC )
execute_process      ( COMMAND python3 -c "import sys; print ( str ( sys.version_info[0] ) + '.' + str ( sys.version_info[1] ) )" OUTPUT_VARIABLE PY_VERSION )
execute_process      ( COMMAND python3 -c "import numpy; print ( numpy.get_include() )" OUTPUT_VARIABLE PY3_INCLUDE_LOC )

##########################################################################################
################################### Check if python3 is installed as python3, terminate this target otherwise
string               ( REGEX REPLACE "\n$" "" PY_VERSION "${PY_VERSION}"                  )

if    ( "${PY_VERSION}" STREQUAL "" )
    message          ( WARNING "Python3 was not detected (i.e. command \"python3\" was not found). ProSHADE installation will assume no python3 is installed and will therefore not produce python3 module." )
    return           (                                                                    )
endif ( "${PY_VERSION}" STREQUAL "" )

if     ( ${PY_VERSION} LESS "3" )
    message          ( WARNING "Python3 was not found (i.e. command \"python3 --version\" returns version < 3). ProSHADE installation will assume no python3 is installed and will therefore not produce python3 module." )
    return           (                                                                    )
endif  ( ${PY_VERSION} LESS "3" )

##########################################################################################
################################### Include the source code and the python include
include_directories  ( ${CMAKE_SOURCE_DIR}/src/proshade                                   )
include_directories  ( ${REAL_PYTHON_INC}                                                 )
include_directories  ( ${PY3_INCLUDE_LOC}/numpy                                           )

##########################################################################################
################################### Get the libpython library location
# Find the python library containing folder options 
execute_process      ( COMMAND python3 -c "import distutils.sysconfig as sysconfig; fpaths = [sysconfig.get_config_var(pv) for pv in ('LIBDIR', 'LIBPL')]; print ( *fpaths )" OUTPUT_VARIABLE REAL_PYTHON_LIB )

# Split the results into a list
separate_arguments ( REAL_PYTHON_LIB )

# For each possible path
foreach    ( POS_PATH ${REAL_PYTHON_LIB} )
    # Remove the trailing endline
    string               ( REGEX REPLACE "\n$" "" POS_PATH "${POS_PATH}" )

    # Manipulate the path to include the libpython name, including the version
    string               ( CONCAT REAL_PYTHON_LIB_HLP ${POS_PATH} "/libpython" ${PY_VERSION} )
    string               ( REGEX REPLACE "\n$" "" REAL_PYTHON_LIB_HLP "${REAL_PYTHON_LIB_HLP}" )

    # Add the correct extension
    if    ( APPLE )
        string           ( CONCAT REAL_PYTHON_LIB_HLP2 ${REAL_PYTHON_LIB_HLP} ".dylib"    )
    else  ( APPLE )
        string           ( CONCAT REAL_PYTHON_LIB_HLP2 ${REAL_PYTHON_LIB_HLP} ".so"       )
    endif ( APPLE )

    # Remove the trailing endline
    string               ( REGEX REPLACE "\n$" "" REAL_PYTHON_LIB_HLP2 "${REAL_PYTHON_LIB_HLP2}" )

    # Found?
    if    ( EXISTS ${REAL_PYTHON_LIB_HLP2} )
	set              ( PYTHON_LIBRARIES ${REAL_PYTHON_LIB_HLP2}                           )
        break            ( )
    endif ( EXISTS ${REAL_PYTHON_LIB_HLP2} )
    
    # Try the same library name, except with 'm' after the version
	if    ( APPLE )
        string           ( CONCAT REAL_PYTHON_LIB_HLP2 ${REAL_PYTHON_LIB_HLP} "m.dylib"    )
    else  ( APPLE )
        string           ( CONCAT REAL_PYTHON_LIB_HLP2 ${REAL_PYTHON_LIB_HLP} "m.so"       )
    endif ( APPLE )

    # Remove the trailing endline
    string               ( REGEX REPLACE "\n$" "" REAL_PYTHON_LIB_HLP2 "${REAL_PYTHON_LIB_HLP2}" )

    # Found?
    if    ( EXISTS ${REAL_PYTHON_LIB_HLP2} )
	set              ( PYTHON_LIBRARIES ${REAL_PYTHON_LIB_HLP2}                           )
        break            ( )
    endif ( EXISTS ${REAL_PYTHON_LIB_HLP2} )
endforeach ( POS_PATH ${REAL_PYTHON_LIB} )

# Check for fails
if    ( NOT EXISTS ${PYTHON_LIBRARIES} )
	set              ( PYTHON_LIBRARIES ""                                                )
	message          ( FATAL_ERROR "Failed to find the libpython library in " ${REAL_PYTHON_LIB} " Please report this case to the author ..." )
endif ( NOT EXISTS ${PYTHON_LIBRARIES} )

##########################################################################################
################################### Define language and specific swig flags
set_source_files_properties ( ${CMAKE_BINARY_DIR}/python3/proshade.i PROPERTIES CPLUSPLUS ON      )
set_source_files_properties ( ${CMAKE_BINARY_DIR}/python3/proshade.i PROPERTIES SWIG_FLAGS "-py3" )

##########################################################################################
################################### Create target using SWIG
if     ( ${CMAKE_MINOR_VERSION} LESS 8 )
	swig_add_module         ( ${PROJECT_NAME}_py3Module LANGUAGE python SOURCES ${CMAKE_BINARY_DIR}/python3/proshade.i )
else   ( ${CMAKE_MINOR_VERSION} LESS 8 )
	swig_add_library        ( ${PROJECT_NAME}_py3Module LANGUAGE python SOURCES ${CMAKE_BINARY_DIR}/python3/proshade.i )
endif  ( ${CMAKE_MINOR_VERSION} LESS 8 )

##########################################################################################
################################### Link the dynamic library to the target
if    ( ${CMAKE_VERSION} GREATER "3.12" )
    target_link_libraries (  ${SWIG_MODULE_NAME} ${PROJECT_NAME}_LIB ${PYTHON_LIBRARIES}  )
else  ( ${CMAKE_VERSION} GREATER "3.12" )
    target_link_libraries ( _${SWIG_MODULE_NAME} ${PROJECT_NAME}_LIB ${PYTHON_LIBRARIES}  )
endif ( ${CMAKE_VERSION} GREATER "3.12" )

##########################################################################################
################################### Change compiler flags to remove SWIG caused warnings
if     ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" )
  ################################# Using Clang
  set ( CMAKE_CXX_FLAGS "-O3 -g -Wall -Wno-macro-redefined"                               )
elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
  ################################# Using GCC
  set ( CMAKE_CXX_FLAGS "-O3 -g -Wall -Wno-c++11-compat"                                  )
elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
  ################################# Using Intel C++
  ### TO BE COMPLETED
elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" )
  ################################# Using Visual Studio C++
  ### TO BE COMPLETED
endif ( )

##########################################################################################
################################### Deal with MacOS rpath
if (APPLE) 
   # Use rpath for ${PROJECT_NAME}_py3Module target
   set_target_properties ( ${PROJECT_NAME}_py3Module PROPERTIES MACOSX_RPATH TRUE         )

   # append directories in the linker search path and outside the project to the INSTALL_RPATH
   set_target_properties ( ${PROJECT_NAME}_py3Module PROPERTIES CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
endif()

##########################################################################################
################################### Set installation runpath value, settings and output names
set   ( PYTHON_INSTALL_RPATH "${MY_INSTALL_LOCATION}/${INSTALL_LIB_DIR}"  "${CMAKE_SOURCE_DIR}/install/lib" )

if    ( ${CMAKE_VERSION} GREATER "3.12" )
    set_target_properties (  ${SWIG_MODULE_NAME} PROPERTIES INSTALL_RPATH "PYTHON_INSTALL_RPATH" )
    set_target_properties (  ${SWIG_MODULE_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE )
    set_target_properties (  ${SWIG_MODULE_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME}  )
else  ( ${CMAKE_VERSION} GREATER "3.12" )
    set_target_properties ( _${SWIG_MODULE_NAME} PROPERTIES INSTALL_RPATH "PYTHON_INSTALL_RPATH" )
    set_target_properties ( _${SWIG_MODULE_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE )
    set_target_properties ( _${SWIG_MODULE_NAME} PROPERTIES OUTPUT_NAME _${PROJECT_NAME}  )
endif ( ${CMAKE_VERSION} GREATER "3.12" )

##########################################################################################
################################### Install rules
if    ( ${CMAKE_VERSION} GREATER "3.12" )
    install           ( TARGETS  ${SWIG_MODULE_NAME} DESTINATION ${MY_PY3_INSTALL_LOCATION} ) 
else  ( ${CMAKE_VERSION} GREATER "3.12" )
    install           ( TARGETS _${SWIG_MODULE_NAME} DESTINATION ${MY_PY3_INSTALL_LOCATION} ) 
endif ( ${CMAKE_VERSION} GREATER "3.12" )
install               ( FILES   ${CMAKE_BINARY_DIR}/cmake/python3/proshade.py DESTINATION ${MY_PY3_INSTALL_LOCATION} ) 

