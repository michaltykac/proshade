##########################################################################################
##########################################################################################
############################### ProSHADE python3 cmake file ##############################
##########################################################################################
##########################################################################################

### Find the SWIG package with the new naming policy
find_package         ( SWIG REQUIRED                                                      )

### Set swig module name 
if    ( ${CMAKE_VERSION} GREATER "3.12" )
    cmake_policy     ( SET CMP0078 OLD                                                    )
endif ( ${CMAKE_VERSION} GREATER "3.12" )

if    ( ${CMAKE_VERSION} GREATER "3.13" )
    cmake_policy     ( SET CMP0086 NEW                                                    )
    set ( SWIG_MODULE_NAME _${PROJECT_NAME}                                               )
endif ( ${CMAKE_VERSION} GREATER "3.13" )

### Include the swig location 
include              ( ${SWIG_USE_FILE}                                                   )

### Copy the files to be pythonised into the building folder
execute_process      ( COMMAND bash -c "mkdir python3" )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/src/proshade/ProSHADE_typedefs.hpp ${CMAKE_BINARY_DIR}/python3/ProSHADE_typedefs.hpp" )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/src/proshade/ProSHADE_settings.hpp ${CMAKE_BINARY_DIR}/python3/ProSHADE_settings.hpp" )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/src/proshade/ProSHADE.hpp          ${CMAKE_BINARY_DIR}/python3/ProSHADE.hpp" )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/src/python/proshade.i              ${CMAKE_BINARY_DIR}/python3/proshade.i" )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/src/python/numpy.i                 ${CMAKE_BINARY_DIR}/python3/numpy.i" )

### Find the python include location
execute_process      ( COMMAND python3 -c "import sysconfig as sysconfig; import sys; sys.stdout.write (sysconfig.get_path(scheme='posix_prefix',name='include'))" OUTPUT_VARIABLE REAL_PYTHON_INC )
execute_process      ( COMMAND python3 -c "import distutils.sysconfig as sysconfig; import sys; sys.stdout.write(sysconfig.get_config_var('LIBDIR'))" OUTPUT_VARIABLE REAL_PYTHON_LIB )
execute_process      ( COMMAND python3 -c "import sys; print ( str ( sys.version_info[0] ) + '.' + str ( sys.version_info[1] ) )" OUTPUT_VARIABLE PY_VERSION )
execute_process      ( COMMAND python3 -c "import numpy; print ( numpy.get_include() )" OUTPUT_VARIABLE PY3_INCLUDE_LOC )

### Find python site packages location
execute_process      ( COMMAND python3 -c "import sys; print ( sys.exec_prefix )" OUTPUT_VARIABLE PY_PREFIX )
string               ( REGEX REPLACE "\n$" "" PY_PREFIX "${PY_PREFIX}"                    )
string               ( CONCAT PYTHON_INSTALL ${PY_PREFIX} "/lib/python"                   )
string               ( REGEX REPLACE "\n$" "" PYTHON_INSTALL "${PYTHON_INSTALL}"          )
string               ( CONCAT PYTHON_INSTALL ${PYTHON_INSTALL} ${PY_VERSION}              )
string               ( REGEX REPLACE "\n$" "" PYTHON_INSTALL "${PYTHON_INSTALL}"          )
string               ( CONCAT PYTHON_INSTALL2 ${PYTHON_INSTALL} "/site-packages"          )
string               ( REGEX REPLACE "\n$" "" PYTHON_INSTALL2 "${PYTHON_INSTALL2}"        )
if    ( EXISTS ${PYTHON_INSTALL2} )
	set ( PY3_INSTALL_LOCATION ${PYTHON_INSTALL2}                                         )
else  ( EXISTS ${PYTHON_INSTALL2} )
	set ( PY3_INSTALL_LOCATION ${PYTHON_INSTALL}                                          )
endif ( EXISTS ${PYTHON_INSTALL2} )

### Include the source code and the python include
include_directories  ( ${CMAKE_SOURCE_DIR}/src/proshade                                   )
include_directories  ( ${REAL_PYTHON_INC}                                                 )
include_directories  ( ${PY3_INCLUDE_LOC}/numpy                                           )

### Paste all info together to get the proper paths to the python library
if     ( APPLE ) 
	string ( CONCAT REAL_PYTHON_LIB_HLP ${REAL_PYTHON_LIB} "/libpython" ${PY_VERSION}     )
	string ( REGEX REPLACE "\n$" "" REAL_PYTHON_LIB_HLP "${REAL_PYTHON_LIB_HLP}"          )
	string ( CONCAT REAL_PYTHON_LIB_HLP2 ${REAL_PYTHON_LIB_HLP} ".dylib"                  )
	string ( REGEX REPLACE "\n$" "" REAL_PYTHON_LIB_HLP2 "${REAL_PYTHON_LIB_HLP2}"        )
	if    ( EXISTS ${REAL_PYTHON_LIB_HLP2} )
		set          ( PYTHON_LIBRARIES ${REAL_PYTHON_LIB_HLP2}                           )
	else  ( EXISTS ${REAL_PYTHON_LIB_HLP2} )
		set          ( PYTHON_LIBRARIES ""                                                )
	endif ( EXISTS ${REAL_PYTHON_LIB_HLP2} )
else   ( APPLE )
	message ( FATAL_ERROR "Pythonisation of ProSHADE not yet implemented outside MacOS" )
endif  ( APPLE )

### Define language and specific swig flags
set_source_files_properties ( ${CMAKE_BINARY_DIR}/python3/proshade.i PROPERTIES CPLUSPLUS ON                          )
set_source_files_properties ( ${CMAKE_BINARY_DIR}/python3/proshade.i PROPERTIES SWIG_FLAGS "-py3"                     )

### Create target using SWIG
if     ( ${CMAKE_MINOR_VERSION} LESS 8 )
	swig_add_module         ( ${PROJECT_NAME} LANGUAGE python SOURCES ${CMAKE_BINARY_DIR}/python3/proshade.i )
else   ( ${CMAKE_MINOR_VERSION} LESS 8 )
	swig_add_library        ( ${PROJECT_NAME} LANGUAGE python SOURCES ${CMAKE_BINARY_DIR}/python3/proshade.i          )
endif  ( ${CMAKE_MINOR_VERSION} LESS 8 )

### Link the dynamic library to the target
swig_link_libraries ( ${PROJECT_NAME} ${PROJECT_NAME}_LIB ${PYTHON_LIBRARIES}             )

### Change compiler flags to remove SWIG caused warnings
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  ################################# Using Clang
  set ( CMAKE_CXX_FLAGS "-O3 -g -Wall -Wno-macro-redefined"                               )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  ################################# Using GCC
  set ( CMAKE_CXX_FLAGS "-O3 -g -Wall -Wno-c++11-compat"                                  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  ################################# Using Intel C++
  ### TO BE COMPLETED
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  ################################# Using Visual Studio C++
  ### TO BE COMPLETED
endif()

### Deal with MacOS rpath
if (APPLE) 
   # Use rpath for _${PROJECT_NAME}.so
   set_target_properties ( _${PROJECT_NAME} PROPERTIES MACOSX_RPATH TRUE                  )

   # append directories in the linker search path and outside the project to the INSTALL_RPATH
   set_target_properties ( _${PROJECT_NAME} PROPERTIES CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

endif()

### Set installation locations
set_target_properties ( _${PROJECT_NAME} PROPERTIES INSTALL_RPATH "${MY_INSTALL_LOCATION}/lib" )
set_target_properties ( _${PROJECT_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE     )

### Install rules
install               ( TARGETS _${PROJECT_NAME} DESTINATION ${MY_PY3_INSTALL_LOCATION}  ) 
install               ( FILES   ${CMAKE_BINARY_DIR}/cmake/python3/proshade.py DESTINATION ${MY_PY3_INSTALL_LOCATION} ) 

