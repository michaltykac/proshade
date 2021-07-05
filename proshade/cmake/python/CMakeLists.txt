##########################################################################################
##########################################################################################
############################### ProSHADE python cmake file ###############################
##########################################################################################
##########################################################################################

##########################################################################################
################################### Add the PyBind11 directory
add_subdirectory        ( ${CMAKE_SOURCE_DIR}/extern/pybind11 build                       )

##########################################################################################
################################### Add the PyBind11 target
pybind11_add_module     ( py${PROJECT_NAME} MODULE ${CMAKE_SOURCE_DIR}/src/python/pyProSHADE_bindings.cpp )

##########################################################################################
################################### Add and link the dependencies
add_dependencies        ( py${PROJECT_NAME} gemmi_lib                                     )
add_dependencies        ( py${PROJECT_NAME} soft2_lib                                     )
add_dependencies        ( py${PROJECT_NAME} ${PROJECT_NAME}_LIB                           )

if     ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  )
	target_link_libraries   ( py${PROJECT_NAME} PRIVATE zlib                              )
	target_link_libraries   ( py${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/winLibs/x64/SOFT2/soft2.lib )
	target_link_libraries   ( py${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/winLibs/x64/FFTW3/libfftw3-3.lib )
	target_link_libraries   ( py${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/winLibs/x64/LAPACK/liblapack.dll.a )
	target_link_libraries   ( py${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/winLibs/x64/LAPACK/libblas.dll.a )
else   ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  )
	target_link_libraries   ( py${PROJECT_NAME} PRIVATE z                                 )
	target_link_libraries   ( py${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/extern/soft-2.0/libsoft1.a )
	target_link_libraries   ( py${PROJECT_NAME} PRIVATE fftw3                             )
	target_link_libraries   ( py${PROJECT_NAME} PRIVATE lapack blas                       )
endif  ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows"  )


##########################################################################################
################################### Include all the required files
include_directories     ( ${CMAKE_SOURCE_DIR}/extern/pybind11/include                     )
include_directories     ( ${CMAKE_SOURCE_DIR}/src/proshade                                )
include_directories     ( ${CMAKE_SOURCE_DIR}/extern/gemmi                                )
include_directories     ( ${CMAKE_SOURCE_DIR}/extern                                      )

##########################################################################################
################################### Install to bin
set_target_properties   ( py${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${PROJECT_NAME}"      )
install                 ( TARGETS py${PROJECT_NAME} DESTINATION ${MY_INSTALL_LOCATION}/pythonModule  )




		
		


