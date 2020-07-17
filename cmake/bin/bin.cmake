##########################################################################################
##########################################################################################
################################# ProSHADE bin cmake file ################################
##########################################################################################
##########################################################################################

##########################################################################################
################################### Find the sources
file                    ( GLOB SOURCES  "${CMAKE_SOURCE_DIR}/src/proshade/*.cpp"          )
file                    ( GLOB EXEC_SRC "${CMAKE_SOURCE_DIR}/src/bin/*.cpp"	              )

##########################################################################################
################################### Link the executable
add_executable          ( ${PROJECT_NAME} ${EXEC_SRC} ${SOURCES}                          )
add_dependencies        ( ${PROJECT_NAME} soft2.0_lib                                     )
target_link_libraries   ( ${PROJECT_NAME} ccp4c mmdb2 )
target_link_libraries   ( ${PROJECT_NAME} fftw3                                           )
target_link_libraries   ( ${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/extern/soft-2.0/libsoft1.a  )
target_link_libraries   ( ${PROJECT_NAME} lapack blas                                     )

##########################################################################################
################################### Set RPATH for the installed executable
set_property            (
              				TARGET ${PROJECT_NAME}
                			PROPERTY INSTALL_RPATH
                			        "${FFTW_LINK}"
                			        "${SOFT_LINK}"
                			        "${LAPACK_LINK}" )

##########################################################################################
################################### Install to bin
#install                 ( TARGETS ${PROJECT_NAME} DESTINATION ${MY_INSTALL_LOCATION}/bin  )
#install                 ( FILES ${MY_INSTALL_LOCATION}/bin/${PROJECT_NAME} DESTINATION ${CMAKE_SOURCE_DIR}/bin PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
