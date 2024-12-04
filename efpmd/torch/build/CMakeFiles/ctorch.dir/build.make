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
CMAKE_SOURCE_DIR = /depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch/build

# Include any dependencies generated for this target.
include CMakeFiles/ctorch.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ctorch.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ctorch.dir/flags.make

CMakeFiles/ctorch.dir/c_libtorch.cc.o: CMakeFiles/ctorch.dir/flags.make
CMakeFiles/ctorch.dir/c_libtorch.cc.o: ../c_libtorch.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ctorch.dir/c_libtorch.cc.o"
	/apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ctorch.dir/c_libtorch.cc.o -c /depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch/c_libtorch.cc

CMakeFiles/ctorch.dir/c_libtorch.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ctorch.dir/c_libtorch.cc.i"
	/apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch/c_libtorch.cc > CMakeFiles/ctorch.dir/c_libtorch.cc.i

CMakeFiles/ctorch.dir/c_libtorch.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ctorch.dir/c_libtorch.cc.s"
	/apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch/c_libtorch.cc -o CMakeFiles/ctorch.dir/c_libtorch.cc.s

# Object files for target ctorch
ctorch_OBJECTS = \
"CMakeFiles/ctorch.dir/c_libtorch.cc.o"

# External object files for target ctorch
ctorch_EXTERNAL_OBJECTS =

libctorch.so: CMakeFiles/ctorch.dir/c_libtorch.cc.o
libctorch.so: CMakeFiles/ctorch.dir/build.make
libctorch.so: /home/paulsk/libtorch/lib/libtorch.so
libctorch.so: /home/paulsk/libtorch/lib/libc10.so
libctorch.so: /home/paulsk/libtorch/lib/libkineto.a
libctorch.so: /home/paulsk/libtorch/lib/libc10.so
libctorch.so: CMakeFiles/ctorch.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library libctorch.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ctorch.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ctorch.dir/build: libctorch.so

.PHONY : CMakeFiles/ctorch.dir/build

CMakeFiles/ctorch.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ctorch.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ctorch.dir/clean

CMakeFiles/ctorch.dir/depend:
	cd /depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch /depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch /depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch/build /depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch/build /depot/lslipche/data/skp/tstgit/efp_monte_carlo/efpmd/torch/build/CMakeFiles/ctorch.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ctorch.dir/depend
