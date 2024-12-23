# libefpConfig.cmake
# ------------------
#
# LIBEFP cmake module.
# This module sets the following variables in your project::
#
#   libefp_FOUND - true if libefp and all required components found on the system
#   libefp_VERSION - libefp version in format Major.Minor.Release
#   libefp_INCLUDE_DIRS - Directory where libefp header is located.
#   libefp_INCLUDE_DIR - same as DIRS
#   libefp_DEFINITIONS - Definitions necessary to use libefp, namely USING_libefp.
#   libefp_LIBRARIES - libefp library to link against.
#   libefp_LIBRARY - same as LIBRARIES
#   libefp_FRAGLIB_DIRS - Directories (list) where EFP fragments are located
#
#
# Available components: shared static ::
#
#   shared - search for only shared library
#   static - search for only static library
#   shallow - search for only fragment library where directory structure has been collapsed
#
#
# Exported targets::
#
# If libefp is found, this module defines the following :prop_tgt:`IMPORTED`
# target. Target is shared _or_ static, so, for both, use separate, not
# overlapping, installations. ::
#
#   libefp::efp - the main libefp library with header & defs attached.
#
#
# Suggested usage::
#
#   find_package(libefp)
#   find_package(libefp 1.5.0 EXACT CONFIG REQUIRED COMPONENTS shared)
#
#
# The following variables can be set to guide the search for this package::
#
#   libefp_DIR - CMake variable, set to directory containing this Config file
#   CMAKE_PREFIX_PATH - CMake variable, set to root directory of this package
#   PATH - environment variable, set to bin directory of this package
#   CMAKE_DISABLE_FIND_PACKAGE_libefp - CMake variable, disables
#     find_package(libefp) when not REQUIRED, perhaps to force internal build


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was libefpConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

set(PN libefp)
set (_valid_components
    static
    shared
    shallow
)

# find includes
unset(_temp_h CACHE)
find_path(_temp_h
          NAMES efp.h
          PATHS ${PACKAGE_PREFIX_DIR}/include
          NO_DEFAULT_PATH)
if(_temp_h)
    set(${PN}_INCLUDE_DIR "${_temp_h}")
    set(${PN}_INCLUDE_DIRS ${${PN}_INCLUDE_DIR})
else()
    set(${PN}_FOUND 0)
    if(NOT CMAKE_REQUIRED_QUIET)
        message(STATUS "${PN}Config missing component: header (${PN}: ${_temp_h})")
    endif()
endif()

# find library: shared, static, or whichever
set(_hold_library_suffixes ${CMAKE_FIND_LIBRARY_SUFFIXES})
list(FIND ${PN}_FIND_COMPONENTS "shared" _seek_shared)
list(FIND ${PN}_FIND_COMPONENTS "static" _seek_static)
if(_seek_shared GREATER -1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX})
elseif(_seek_static GREATER -1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()
unset(_temp CACHE)
find_library(_temp
             NAMES efp
             PATHS ${PACKAGE_PREFIX_DIR}/lib64
             NO_DEFAULT_PATH)
if(_temp)
    set(${PN}_LIBRARY "${_temp}")
    if(_seek_shared GREATER -1)
        set(${PN}_shared_FOUND 1)
    elseif(_seek_static GREATER -1)
        set(${PN}_static_FOUND 1)
    endif()
else()
    if(_seek_shared GREATER -1)
        if(NOT CMAKE_REQUIRED_QUIET)
            message(STATUS "${PN}Config missing component: shared library (${PN}: ${_temp})")
        endif()
    elseif(_seek_static GREATER -1)
        if(NOT CMAKE_REQUIRED_QUIET)
            message(STATUS "${PN}Config missing component: static library (${PN}: ${_temp})")
        endif()
    else()
        set(${PN}_FOUND 0)
        if(NOT CMAKE_REQUIRED_QUIET)
            message(STATUS "${PN}Config missing component: library (${PN}: ${_temp})")
        endif()
    endif()
endif()
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_hold_library_suffixes})
set(${PN}_LIBRARIES ${${PN}_LIBRARY})
set(${PN}_DEFINITIONS USING_${PN})

# find fraglibs
list(FIND ${PN}_FIND_COMPONENTS "shallow" _seek_shallow)
string(REGEX REPLACE "([^;]+)" "${PACKAGE_PREFIX_DIR}/share/${PN}/\\1" ${PN}_FRAGLIB_DIRS "fraglib;fraglib/databases")
if(_seek_shallow GREATER -1)
    list(LENGTH ${PN}_FRAGLIB_DIRS _temp_len)
    if(_temp_len EQUAL 1)
        set(${PN}_shallow_FOUND 1)
    else()
        if(NOT CMAKE_REQUIRED_QUIET)
            message(STATUS "${PN}Config missing component: shallow fraglib (${PN}: ${${PN}_FRAGLIB_DIRS})")
        endif()
    endif()
endif()

check_required_components(${PN})

# make detectable the FindTarget*.cmake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already built this as a subproject
#-----------------------------------------------------------------------------
if(NOT TARGET ${PN}::efp)
    include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets.cmake")

    include(CMakeFindDependencyMacro)
    if(NOT TARGET tgt::lapack)
        find_dependency(TargetLAPACK)
    endif()
endif()

