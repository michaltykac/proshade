##########################################################################################
##########################################################################################
############################### ProSHADE library cmake file ##############################
##########################################################################################
##########################################################################################

##########################################################################################
################################### Add the library target
add_library             ( ${PROJECT_NAME}_LIB SHARED ${OBJS}                              )

##########################################################################################
### Add and link the dependencies
add_dependencies        ( ${PROJECT_NAME} gemmi_lib                                       )
add_dependencies        ( ${PROJECT_NAME} soft2_lib                                       )
target_link_libraries   ( ${PROJECT_NAME}_LIB z                                           )
target_link_libraries   ( ${PROJECT_NAME}_LIB ${CMAKE_SOURCE_DIR}/extern/soft-2.0/libsoft1.a )
target_link_libraries   ( ${PROJECT_NAME}_LIB fftw3                                       )
target_link_libraries   ( ${PROJECT_NAME}_LIB lapack blas                                 )

##########################################################################################
################################### Set RPATH for the installed library
set_property            (  TARGET ${PROJECT_NAME}_LIB
                           PROPERTY INSTALL_RPATH
                                 "${FFTW_LINK}"
                                 "${SOFT_LINK}"
                                 "${LAPACK_LINK}" )

##########################################################################################
################################### Deal with MacOS RPATH
if (APPLE) 
   set_target_properties ( ${PROJECT_NAME}_LIB PROPERTIES MACOSX_RPATH TRUE               )
   set_target_properties ( ${PROJECT_NAME}_LIB PROPERTIES CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

##########################################################################################
################################### Set paths
set_target_properties ( ${PROJECT_NAME}_LIB PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE   )

##########################################################################################
################################### Set installation
set_target_properties ( ${PROJECT_NAME}_LIB PROPERTIES OUTPUT_NAME "${PROJECT_NAME}"      )
install               ( TARGETS ${PROJECT_NAME}_LIB DESTINATION ${MY_INSTALL_LOCATION}/${INSTALL_LIB_DIR} ) 
install               ( FILES   ${CMAKE_SOURCE_DIR}/src/proshade/ProSHADE.hpp DESTINATION ${MY_INSTALL_LOCATION}/${INSTALL_INC_DIR} ) 