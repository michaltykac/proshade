##########################################################################################
##########################################################################################
################################# ProSHADE bin cmake file ################################
##########################################################################################
##########################################################################################

##########################################################################################
################################### Add the executable target
add_executable          ( ${PROJECT_NAME} ${EXEC_SRC} ${OBJS}                             )

##########################################################################################
### Add and link the dependencies
add_dependencies        ( ${PROJECT_NAME} gemmi_lib                                       )
add_dependencies        ( ${PROJECT_NAME} soft2_lib                                       )
add_dependencies        ( ${PROJECT_NAME} mmdb2_lib                                       )
add_dependencies        ( ${PROJECT_NAME} ccp4_lib                                        )
target_link_libraries   ( ${PROJECT_NAME} ccp4c mmdb2 z )
target_link_libraries   ( ${PROJECT_NAME} fftw3                                           )
target_link_libraries   ( ${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/extern/soft-2.0/libsoft1.a  )
target_link_libraries   ( ${PROJECT_NAME} lapack blas                                     )

##########################################################################################
################################### Set RPATH for the installed executable
set_property            (
              				TARGET ${PROJECT_NAME}
                			PROPERTY INSTALL_RPATH
                			        "${CMAKE_SOURCE_DIR}/install"
                			        "${FFTW_LINK}"
                			        "${SOFT_LINK}"
                			        "${LAPACK_LINK}" )

##########################################################################################
################################### Install to bin
install                 ( TARGETS ${PROJECT_NAME} DESTINATION ${MY_INSTALL_LOCATION}/${INSTALL_BIN_DIR}  )