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
################################### Add the executable target
add_executable          ( ${PROJECT_NAME} ${EXEC_SRC} ${SOURCES}                          )

##########################################################################################
################################### Link the executable
add_dependencies        ( ${PROJECT_NAME} gemmi_lib                                       )
add_dependencies        ( ${PROJECT_NAME} soft2_lib                                       )
target_link_libraries   ( ${PROJECT_NAME} z                                               )
target_link_libraries   ( ${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/extern/soft-2.0/libsoft1.a  )
target_link_libraries   ( ${PROJECT_NAME} fftw3                                           )
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
