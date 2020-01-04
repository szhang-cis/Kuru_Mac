# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/skoll/FEM/vulcan

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/skoll/FEM/vulcan/build

# Include any dependencies generated for this target.
include src/genren/trunk/CMakeFiles/malla3d.dir/depend.make

# Include the progress variables for this target.
include src/genren/trunk/CMakeFiles/malla3d.dir/progress.make

# Include the compile flags for this target's objects.
include src/genren/trunk/CMakeFiles/malla3d.dir/flags.make

src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o: src/genren/trunk/CMakeFiles/malla3d.dir/flags.make
src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o: ../src/genren/trunk/malla3d.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/skoll/FEM/vulcan/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o"
	cd /media/skoll/FEM/vulcan/build/src/genren/trunk && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /media/skoll/FEM/vulcan/src/genren/trunk/malla3d.f -o CMakeFiles/malla3d.dir/malla3d.f.o

src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/malla3d.dir/malla3d.f.i"
	cd /media/skoll/FEM/vulcan/build/src/genren/trunk && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /media/skoll/FEM/vulcan/src/genren/trunk/malla3d.f > CMakeFiles/malla3d.dir/malla3d.f.i

src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/malla3d.dir/malla3d.f.s"
	cd /media/skoll/FEM/vulcan/build/src/genren/trunk && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /media/skoll/FEM/vulcan/src/genren/trunk/malla3d.f -o CMakeFiles/malla3d.dir/malla3d.f.s

src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o.requires:

.PHONY : src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o.requires

src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o.provides: src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o.requires
	$(MAKE) -f src/genren/trunk/CMakeFiles/malla3d.dir/build.make src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o.provides.build
.PHONY : src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o.provides

src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o.provides.build: src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o


# Object files for target malla3d
malla3d_OBJECTS = \
"CMakeFiles/malla3d.dir/malla3d.f.o"

# External object files for target malla3d
malla3d_EXTERNAL_OBJECTS =

src/genren/trunk/malla3d: src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o
src/genren/trunk/malla3d: src/genren/trunk/CMakeFiles/malla3d.dir/build.make
src/genren/trunk/malla3d: src/genren/trunk/CMakeFiles/malla3d.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/skoll/FEM/vulcan/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable malla3d"
	cd /media/skoll/FEM/vulcan/build/src/genren/trunk && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/malla3d.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/genren/trunk/CMakeFiles/malla3d.dir/build: src/genren/trunk/malla3d

.PHONY : src/genren/trunk/CMakeFiles/malla3d.dir/build

src/genren/trunk/CMakeFiles/malla3d.dir/requires: src/genren/trunk/CMakeFiles/malla3d.dir/malla3d.f.o.requires

.PHONY : src/genren/trunk/CMakeFiles/malla3d.dir/requires

src/genren/trunk/CMakeFiles/malla3d.dir/clean:
	cd /media/skoll/FEM/vulcan/build/src/genren/trunk && $(CMAKE_COMMAND) -P CMakeFiles/malla3d.dir/cmake_clean.cmake
.PHONY : src/genren/trunk/CMakeFiles/malla3d.dir/clean

src/genren/trunk/CMakeFiles/malla3d.dir/depend:
	cd /media/skoll/FEM/vulcan/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/skoll/FEM/vulcan /media/skoll/FEM/vulcan/src/genren/trunk /media/skoll/FEM/vulcan/build /media/skoll/FEM/vulcan/build/src/genren/trunk /media/skoll/FEM/vulcan/build/src/genren/trunk/CMakeFiles/malla3d.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/genren/trunk/CMakeFiles/malla3d.dir/depend
