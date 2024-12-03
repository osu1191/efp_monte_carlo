cmake_minimum_required(VERSION 3.10)
project(libefp VERSION 1.7.2 LANGUAGES C)

# Add subdirectories
#find_package(TargetLAPACK REQUIRED)

set(libefp_AUTHORS "Ilya A. Kaliman, Lori A. Burns, Dmitry Morozov, Carlos H. Borca, Yen (Terri) Bui, Yongbin Kim, Lyudmila V. Slipchenko")
set(libefp_DESCRIPTION "Parallel implementation of the Effective Fragment Potential method")
set(libefp_URL "https://github.com/libefp2/libefp")
set(libefp_LICENSE "BSD 2-clause")
set(CMAKE_INSTALL_PREFIX "/depot/lslipche/data/skp/efp_monte_carlo")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(psi4OptionsTools)
option_with_default(CMAKE_BUILD_TYPE "Build type" Release)
option_with_print(BUILD_SHARED_LIBS "Build final library as shared, not static" ON)

include(autocmake_safeguards)
include(autocmake_omp)
include(custom_static_library)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(PN ${PROJECT_NAME})
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")

add_subdirectory(src)
add_subdirectory(efpmd)





 
set(FRAGLIB_DATADIRS "")
file(GLOB_RECURSE _dotefps RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "fraglib/*.efp")
foreach(_dotefp ${_dotefps})
    get_filename_component(_efpdir ${_dotefp} DIRECTORY)
    get_filename_component(_efpfile ${_dotefp} NAME)
    file(READ ${_dotefp} _strefpL)
    string(REPLACE "_L" "" _strefp ${_strefpL})

    if(FRAGLIB_DEEP)
        set(_destdir ${_efpdir})
    else()
        set(_destdir fraglib)
    endif()
    if(FRAGLIB_UNDERSCORE_L)
        set(_destcontents ${_strefpL})
    else()
        set(_destcontents ${_strefp})
    endif()

    list(APPEND FRAGLIB_DATADIRS ${_destdir})
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${_destdir}/${_efpfile} ${_destcontents})
endforeach()
list(REMOVE_DUPLICATES FRAGLIB_DATADIRS)

install(FILES fraglib/makefp.inp
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${PN})
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fraglib
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${PN})


