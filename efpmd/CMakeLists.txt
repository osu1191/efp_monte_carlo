# Set source files for the executable

set(SOURCES
    src/cfg.c
    src/common.c
    src/efield.c
    src/energy.c
    src/grad.c
    src/gtest.c
    src/hess.c
    src/md.c
    src/mc.c
    src/msg.c
    src/opt.c
    src/parse.c
    src/rand.c
    src/main.c
    src/sp.c
)

# Create the executable target

add_executable(efpmd ${SOURCES})

#add_library(efpmd ${SOURCES})
  
add_subdirectory(libopt)
add_subdirectory(libff)

# Set include directories
target_include_directories(efpmd PRIVATE 
		${CMAKE_CURRENT_SOURCE_DIR}/src
		${CMAKE_CURRENT_SOURCE_DIR}/libopt
		${CMAKE_CURRENT_SOURCE_DIR}/libff
		)

# Link libraries
target_link_libraries(efpmd PRIVATE efp ff opt m gfortran)
#target_link_libraries(efpmd PRIVATE efp ${CMAKE_CURRENT_SOURCE_DIR}/libff/libff.a ${CMAKE_CURRENT_SOURCE_DIR}/libopt/libopt.a m)
