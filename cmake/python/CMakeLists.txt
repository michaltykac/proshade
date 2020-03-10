##########################################################################################
##########################################################################################
################################ ProSHADE python cmake file ##############################
##########################################################################################
##########################################################################################

### Find the SWIG package
find_package         ( SWIG REQUIRED                                                      )

### Include the swig location 
include              ( ${SWIG_USE_FILE}                                                   )

### Copy the files to be pythonised into the building folder
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/../src/proshade/ProSHADE_typedefs.hpp ${CMAKE_SOURCE_DIR}/python/ProSHADE_typedefs.hpp" )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/../src/proshade/ProSHADE_settings.hpp ${CMAKE_SOURCE_DIR}/python/ProSHADE_settings.hpp" )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/../src/proshade/ProSHADE.hpp          ${CMAKE_SOURCE_DIR}/python/ProSHADE.hpp" )
execute_process      ( COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/../src/python/proshade.i              ${CMAKE_SOURCE_DIR}/python/proshade.i" )

### Find the python include location
execute_process      ( COMMAND python -c "import sysconfig as sysconfig; import sys; sys.stdout.write (sysconfig.get_path(scheme='posix_prefix',name='include'))" OUTPUT_VARIABLE REAL_PYTHON_INC )
execute_process      ( COMMAND python -c "import distutils.sysconfig as sysconfig; import sys; sys.stdout.write(sysconfig.get_config_var('LIBDIR'))" OUTPUT_VARIABLE REAL_PYTHON_LIB )
execute_process      ( COMMAND python -c "import sys; print ( str ( sys.version_info[0] ) + '.' + str ( sys.version_info[1] ) )" OUTPUT_VARIABLE PY_VERSION )

### Include the source code and the python include
include_directories  ( ${PROJECT_SOURCE_DIR}/../src/proshade                              )
include_directories  ( ${REAL_PYTHON_INC}                                                 )

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

### Define language
set_source_files_properties ( proshade.i PROPERTIES CPLUSPLUS ON                          )

### Create target using SWIG
if     ( ${CMAKE_MINOR_VERSION} LESS 8 )
	swig_add_module         ( ${PROJECT_NAME}2 LANGUAGE python SOURCES proshade.i         )
else   ( ${CMAKE_MINOR_VERSION} LESS 8 )
	swig_add_library        ( ${PROJECT_NAME}2 LANGUAGE python SOURCES proshade.i         )
endif  ( ${CMAKE_MINOR_VERSION} LESS 8 )

### Link the dynamic library to the target
swig_link_libraries ( ${PROJECT_NAME}2 ${PROJECT_NAME}_LIB ${PYTHON_LIBRARIES}            )

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
   set_target_properties ( _${PROJECT_NAME}2 PROPERTIES MACOSX_RPATH TRUE                 )

   # append directories in the linker search path and outside the project to the INSTALL_RPATH
   set_target_properties ( _${PROJECT_NAME}2 PROPERTIES CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

endif()

### Set installation locations
set_target_properties ( _${PROJECT_NAME}2 PROPERTIES INSTALL_RPATH "${MY_INSTALL_LOCATION}/lib" )
set_target_properties ( _${PROJECT_NAME}2 PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE     )

### Install rules
install               ( TARGETS _${PROJECT_NAME}2 DESTINATION ${MY_PY2_INSTALL_LOCATION}  ) 
install               ( FILES   ${CMAKE_SOURCE_DIR}/python/proshade.py DESTINATION ${MY_PY2_INSTALL_LOCATION} ) 

