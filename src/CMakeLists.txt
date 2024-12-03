# Set source files for the library
set(SOURCES
    aidisp.c  
    balance.c  
    clapack.c  
    disp.c  
    efp.c  
    elec.c  
    electerms.c  
    int.c  
    log.c  
    parse.c  
    pol.c  
    poldirect.c  
    stream.c  
    swf.c  
    util.c  
    xr.c
)

# Create the library target
add_library(efp ${SOURCES})

# Set include directories
target_include_directories(efp PUBLIC 
			${CMAKE_CURRENT_SOURCE_DIR}
			${CMAKE_SOURCE_DIR}/efpmd/src
			${CMAKE_SOURCE_DIR}/efpmd/libff
			${CMAKE_SOURCE_DIR}/efpmd/libopt
			)

# Link libraries if needed
target_link_libraries(efp PRIVATE lapack blas)

# Optionally, install the library
#install(TARGETS efp DESTINATION lib)
#install(FILES efp.h DESTINATION include)

install(TARGETS efp
        EXPORT "${PN}Targets"
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})


install(FILES efp.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_package_config_file(${CMAKE_SOURCE_DIR}/cmake/${PN}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
                                 VERSION ${${PN}_VERSION}
                                 COMPATIBILITY AnyNewerVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindTargetLAPACK.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

#install(EXPORT "${PN}Targets"
#        NAMESPACE "${PN}::"
#        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
