# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /apps/spack/bell/apps/cmake/3.18.2-gcc-10.2.0-z3hxs65/bin/cmake

# The command to remove a file.
RM = /apps/spack/bell/apps/cmake/3.18.2-gcc-10.2.0-z3hxs65/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff/build

# Include any dependencies generated for this target.
include CMakeFiles/ff.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ff.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ff.dir/flags.make

CMakeFiles/ff.dir/ff.c.o: CMakeFiles/ff.dir/flags.make
CMakeFiles/ff.dir/ff.c.o: ../ff.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/ff.dir/ff.c.o"
	/apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/ff.dir/ff.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff/ff.c

CMakeFiles/ff.dir/ff.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ff.dir/ff.c.i"
	/apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff/ff.c > CMakeFiles/ff.dir/ff.c.i

CMakeFiles/ff.dir/ff.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ff.dir/ff.c.s"
	/apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff/ff.c -o CMakeFiles/ff.dir/ff.c.s

# Object files for target ff
ff_OBJECTS = \
"CMakeFiles/ff.dir/ff.c.o"

# External object files for target ff
ff_EXTERNAL_OBJECTS =

libff.a: CMakeFiles/ff.dir/ff.c.o
libff.a: CMakeFiles/ff.dir/build.make
libff.a: CMakeFiles/ff.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library libff.a"
	$(CMAKE_COMMAND) -P CMakeFiles/ff.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ff.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ff.dir/build: libff.a

.PHONY : CMakeFiles/ff.dir/build

CMakeFiles/ff.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ff.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ff.dir/clean

CMakeFiles/ff.dir/depend:
	cd /depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff /depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff /depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff/build /depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff/build /depot/lslipche/data/skp/efp_monte_carlo/efpmd/libff/build/CMakeFiles/ff.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ff.dir/depend

