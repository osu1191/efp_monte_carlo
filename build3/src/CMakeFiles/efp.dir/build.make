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
CMAKE_SOURCE_DIR = /depot/lslipche/data/skp/efp_monte_carlo

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /depot/lslipche/data/skp/efp_monte_carlo/build3

# Include any dependencies generated for this target.
include src/CMakeFiles/efp.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/efp.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/efp.dir/flags.make

src/CMakeFiles/efp.dir/aidisp.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/aidisp.c.o: ../src/aidisp.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/CMakeFiles/efp.dir/aidisp.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/aidisp.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/aidisp.c

src/CMakeFiles/efp.dir/aidisp.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/aidisp.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/aidisp.c > CMakeFiles/efp.dir/aidisp.c.i

src/CMakeFiles/efp.dir/aidisp.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/aidisp.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/aidisp.c -o CMakeFiles/efp.dir/aidisp.c.s

src/CMakeFiles/efp.dir/balance.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/balance.c.o: ../src/balance.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/CMakeFiles/efp.dir/balance.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/balance.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/balance.c

src/CMakeFiles/efp.dir/balance.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/balance.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/balance.c > CMakeFiles/efp.dir/balance.c.i

src/CMakeFiles/efp.dir/balance.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/balance.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/balance.c -o CMakeFiles/efp.dir/balance.c.s

src/CMakeFiles/efp.dir/clapack.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/clapack.c.o: ../src/clapack.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/CMakeFiles/efp.dir/clapack.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/clapack.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/clapack.c

src/CMakeFiles/efp.dir/clapack.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/clapack.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/clapack.c > CMakeFiles/efp.dir/clapack.c.i

src/CMakeFiles/efp.dir/clapack.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/clapack.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/clapack.c -o CMakeFiles/efp.dir/clapack.c.s

src/CMakeFiles/efp.dir/disp.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/disp.c.o: ../src/disp.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object src/CMakeFiles/efp.dir/disp.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/disp.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/disp.c

src/CMakeFiles/efp.dir/disp.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/disp.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/disp.c > CMakeFiles/efp.dir/disp.c.i

src/CMakeFiles/efp.dir/disp.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/disp.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/disp.c -o CMakeFiles/efp.dir/disp.c.s

src/CMakeFiles/efp.dir/efp.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/efp.c.o: ../src/efp.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object src/CMakeFiles/efp.dir/efp.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/efp.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/efp.c

src/CMakeFiles/efp.dir/efp.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/efp.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/efp.c > CMakeFiles/efp.dir/efp.c.i

src/CMakeFiles/efp.dir/efp.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/efp.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/efp.c -o CMakeFiles/efp.dir/efp.c.s

src/CMakeFiles/efp.dir/elec.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/elec.c.o: ../src/elec.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object src/CMakeFiles/efp.dir/elec.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/elec.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/elec.c

src/CMakeFiles/efp.dir/elec.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/elec.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/elec.c > CMakeFiles/efp.dir/elec.c.i

src/CMakeFiles/efp.dir/elec.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/elec.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/elec.c -o CMakeFiles/efp.dir/elec.c.s

src/CMakeFiles/efp.dir/electerms.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/electerms.c.o: ../src/electerms.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object src/CMakeFiles/efp.dir/electerms.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/electerms.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/electerms.c

src/CMakeFiles/efp.dir/electerms.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/electerms.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/electerms.c > CMakeFiles/efp.dir/electerms.c.i

src/CMakeFiles/efp.dir/electerms.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/electerms.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/electerms.c -o CMakeFiles/efp.dir/electerms.c.s

src/CMakeFiles/efp.dir/int.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/int.c.o: ../src/int.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object src/CMakeFiles/efp.dir/int.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/int.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/int.c

src/CMakeFiles/efp.dir/int.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/int.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/int.c > CMakeFiles/efp.dir/int.c.i

src/CMakeFiles/efp.dir/int.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/int.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/int.c -o CMakeFiles/efp.dir/int.c.s

src/CMakeFiles/efp.dir/log.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/log.c.o: ../src/log.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object src/CMakeFiles/efp.dir/log.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/log.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/log.c

src/CMakeFiles/efp.dir/log.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/log.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/log.c > CMakeFiles/efp.dir/log.c.i

src/CMakeFiles/efp.dir/log.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/log.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/log.c -o CMakeFiles/efp.dir/log.c.s

src/CMakeFiles/efp.dir/parse.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/parse.c.o: ../src/parse.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object src/CMakeFiles/efp.dir/parse.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/parse.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/parse.c

src/CMakeFiles/efp.dir/parse.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/parse.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/parse.c > CMakeFiles/efp.dir/parse.c.i

src/CMakeFiles/efp.dir/parse.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/parse.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/parse.c -o CMakeFiles/efp.dir/parse.c.s

src/CMakeFiles/efp.dir/pol.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/pol.c.o: ../src/pol.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object src/CMakeFiles/efp.dir/pol.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/pol.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/pol.c

src/CMakeFiles/efp.dir/pol.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/pol.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/pol.c > CMakeFiles/efp.dir/pol.c.i

src/CMakeFiles/efp.dir/pol.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/pol.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/pol.c -o CMakeFiles/efp.dir/pol.c.s

src/CMakeFiles/efp.dir/poldirect.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/poldirect.c.o: ../src/poldirect.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building C object src/CMakeFiles/efp.dir/poldirect.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/poldirect.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/poldirect.c

src/CMakeFiles/efp.dir/poldirect.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/poldirect.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/poldirect.c > CMakeFiles/efp.dir/poldirect.c.i

src/CMakeFiles/efp.dir/poldirect.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/poldirect.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/poldirect.c -o CMakeFiles/efp.dir/poldirect.c.s

src/CMakeFiles/efp.dir/stream.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/stream.c.o: ../src/stream.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building C object src/CMakeFiles/efp.dir/stream.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/stream.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/stream.c

src/CMakeFiles/efp.dir/stream.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/stream.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/stream.c > CMakeFiles/efp.dir/stream.c.i

src/CMakeFiles/efp.dir/stream.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/stream.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/stream.c -o CMakeFiles/efp.dir/stream.c.s

src/CMakeFiles/efp.dir/swf.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/swf.c.o: ../src/swf.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building C object src/CMakeFiles/efp.dir/swf.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/swf.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/swf.c

src/CMakeFiles/efp.dir/swf.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/swf.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/swf.c > CMakeFiles/efp.dir/swf.c.i

src/CMakeFiles/efp.dir/swf.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/swf.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/swf.c -o CMakeFiles/efp.dir/swf.c.s

src/CMakeFiles/efp.dir/util.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/util.c.o: ../src/util.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building C object src/CMakeFiles/efp.dir/util.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/util.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/util.c

src/CMakeFiles/efp.dir/util.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/util.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/util.c > CMakeFiles/efp.dir/util.c.i

src/CMakeFiles/efp.dir/util.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/util.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/util.c -o CMakeFiles/efp.dir/util.c.s

src/CMakeFiles/efp.dir/xr.c.o: src/CMakeFiles/efp.dir/flags.make
src/CMakeFiles/efp.dir/xr.c.o: ../src/xr.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building C object src/CMakeFiles/efp.dir/xr.c.o"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/efp.dir/xr.c.o -c /depot/lslipche/data/skp/efp_monte_carlo/src/xr.c

src/CMakeFiles/efp.dir/xr.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/efp.dir/xr.c.i"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /depot/lslipche/data/skp/efp_monte_carlo/src/xr.c > CMakeFiles/efp.dir/xr.c.i

src/CMakeFiles/efp.dir/xr.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/efp.dir/xr.c.s"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && /apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /depot/lslipche/data/skp/efp_monte_carlo/src/xr.c -o CMakeFiles/efp.dir/xr.c.s

# Object files for target efp
efp_OBJECTS = \
"CMakeFiles/efp.dir/aidisp.c.o" \
"CMakeFiles/efp.dir/balance.c.o" \
"CMakeFiles/efp.dir/clapack.c.o" \
"CMakeFiles/efp.dir/disp.c.o" \
"CMakeFiles/efp.dir/efp.c.o" \
"CMakeFiles/efp.dir/elec.c.o" \
"CMakeFiles/efp.dir/electerms.c.o" \
"CMakeFiles/efp.dir/int.c.o" \
"CMakeFiles/efp.dir/log.c.o" \
"CMakeFiles/efp.dir/parse.c.o" \
"CMakeFiles/efp.dir/pol.c.o" \
"CMakeFiles/efp.dir/poldirect.c.o" \
"CMakeFiles/efp.dir/stream.c.o" \
"CMakeFiles/efp.dir/swf.c.o" \
"CMakeFiles/efp.dir/util.c.o" \
"CMakeFiles/efp.dir/xr.c.o"

# External object files for target efp
efp_EXTERNAL_OBJECTS =

src/libefp.a: src/CMakeFiles/efp.dir/aidisp.c.o
src/libefp.a: src/CMakeFiles/efp.dir/balance.c.o
src/libefp.a: src/CMakeFiles/efp.dir/clapack.c.o
src/libefp.a: src/CMakeFiles/efp.dir/disp.c.o
src/libefp.a: src/CMakeFiles/efp.dir/efp.c.o
src/libefp.a: src/CMakeFiles/efp.dir/elec.c.o
src/libefp.a: src/CMakeFiles/efp.dir/electerms.c.o
src/libefp.a: src/CMakeFiles/efp.dir/int.c.o
src/libefp.a: src/CMakeFiles/efp.dir/log.c.o
src/libefp.a: src/CMakeFiles/efp.dir/parse.c.o
src/libefp.a: src/CMakeFiles/efp.dir/pol.c.o
src/libefp.a: src/CMakeFiles/efp.dir/poldirect.c.o
src/libefp.a: src/CMakeFiles/efp.dir/stream.c.o
src/libefp.a: src/CMakeFiles/efp.dir/swf.c.o
src/libefp.a: src/CMakeFiles/efp.dir/util.c.o
src/libefp.a: src/CMakeFiles/efp.dir/xr.c.o
src/libefp.a: src/CMakeFiles/efp.dir/build.make
src/libefp.a: src/CMakeFiles/efp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/depot/lslipche/data/skp/efp_monte_carlo/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Linking C static library libefp.a"
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && $(CMAKE_COMMAND) -P CMakeFiles/efp.dir/cmake_clean_target.cmake
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/efp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/efp.dir/build: src/libefp.a

.PHONY : src/CMakeFiles/efp.dir/build

src/CMakeFiles/efp.dir/clean:
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3/src && $(CMAKE_COMMAND) -P CMakeFiles/efp.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/efp.dir/clean

src/CMakeFiles/efp.dir/depend:
	cd /depot/lslipche/data/skp/efp_monte_carlo/build3 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /depot/lslipche/data/skp/efp_monte_carlo /depot/lslipche/data/skp/efp_monte_carlo/src /depot/lslipche/data/skp/efp_monte_carlo/build3 /depot/lslipche/data/skp/efp_monte_carlo/build3/src /depot/lslipche/data/skp/efp_monte_carlo/build3/src/CMakeFiles/efp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/efp.dir/depend

