# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /home/mehdi/Libs/cmake-3.15.0-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /home/mehdi/Libs/cmake-3.15.0-Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/src/main.f90.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/main.f90.o: ../src/main.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/main.dir/src/main.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/main.f90 -o CMakeFiles/main.dir/src/main.f90.o

CMakeFiles/main.dir/src/main.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.dir/src/main.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/main.f90 > CMakeFiles/main.dir/src/main.f90.i

CMakeFiles/main.dir/src/main.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.dir/src/main.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/main.f90 -o CMakeFiles/main.dir/src/main.f90.s

CMakeFiles/main.dir/src/m_interpolations.f90.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/m_interpolations.f90.o: ../src/m_interpolations.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/main.dir/src/m_interpolations.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_interpolations.f90 -o CMakeFiles/main.dir/src/m_interpolations.f90.o

CMakeFiles/main.dir/src/m_interpolations.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.dir/src/m_interpolations.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_interpolations.f90 > CMakeFiles/main.dir/src/m_interpolations.f90.i

CMakeFiles/main.dir/src/m_interpolations.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.dir/src/m_interpolations.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_interpolations.f90 -o CMakeFiles/main.dir/src/m_interpolations.f90.s

CMakeFiles/main.dir/src/m_least_square.f90.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/m_least_square.f90.o: ../src/m_least_square.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object CMakeFiles/main.dir/src/m_least_square.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_least_square.f90 -o CMakeFiles/main.dir/src/m_least_square.f90.o

CMakeFiles/main.dir/src/m_least_square.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.dir/src/m_least_square.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_least_square.f90 > CMakeFiles/main.dir/src/m_least_square.f90.i

CMakeFiles/main.dir/src/m_least_square.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.dir/src/m_least_square.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_least_square.f90 -o CMakeFiles/main.dir/src/m_least_square.f90.s

CMakeFiles/main.dir/src/m_math.f90.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/m_math.f90.o: ../src/m_math.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object CMakeFiles/main.dir/src/m_math.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_math.f90 -o CMakeFiles/main.dir/src/m_math.f90.o

CMakeFiles/main.dir/src/m_math.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.dir/src/m_math.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_math.f90 > CMakeFiles/main.dir/src/m_math.f90.i

CMakeFiles/main.dir/src/m_math.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.dir/src/m_math.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_math.f90 -o CMakeFiles/main.dir/src/m_math.f90.s

CMakeFiles/main.dir/src/mod_test.f90.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/mod_test.f90.o: ../src/mod_test.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building Fortran object CMakeFiles/main.dir/src/mod_test.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/mod_test.f90 -o CMakeFiles/main.dir/src/mod_test.f90.o

CMakeFiles/main.dir/src/mod_test.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.dir/src/mod_test.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/mod_test.f90 > CMakeFiles/main.dir/src/mod_test.f90.i

CMakeFiles/main.dir/src/mod_test.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.dir/src/mod_test.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/mod_test.f90 -o CMakeFiles/main.dir/src/mod_test.f90.s

CMakeFiles/main.dir/src/mod_write.f90.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/mod_write.f90.o: ../src/mod_write.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building Fortran object CMakeFiles/main.dir/src/mod_write.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/mod_write.f90 -o CMakeFiles/main.dir/src/mod_write.f90.o

CMakeFiles/main.dir/src/mod_write.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.dir/src/mod_write.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/mod_write.f90 > CMakeFiles/main.dir/src/mod_write.f90.i

CMakeFiles/main.dir/src/mod_write.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.dir/src/mod_write.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/mod_write.f90 -o CMakeFiles/main.dir/src/mod_write.f90.s

CMakeFiles/main.dir/src/m_tests.f90.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/m_tests.f90.o: ../src/m_tests.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building Fortran object CMakeFiles/main.dir/src/m_tests.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_tests.f90 -o CMakeFiles/main.dir/src/m_tests.f90.o

CMakeFiles/main.dir/src/m_tests.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.dir/src/m_tests.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_tests.f90 > CMakeFiles/main.dir/src/m_tests.f90.i

CMakeFiles/main.dir/src/m_tests.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.dir/src/m_tests.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/src/m_tests.f90 -o CMakeFiles/main.dir/src/m_tests.f90.s

# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/src/main.f90.o" \
"CMakeFiles/main.dir/src/m_interpolations.f90.o" \
"CMakeFiles/main.dir/src/m_least_square.f90.o" \
"CMakeFiles/main.dir/src/m_math.f90.o" \
"CMakeFiles/main.dir/src/mod_test.f90.o" \
"CMakeFiles/main.dir/src/mod_write.f90.o" \
"CMakeFiles/main.dir/src/m_tests.f90.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

main: CMakeFiles/main.dir/src/main.f90.o
main: CMakeFiles/main.dir/src/m_interpolations.f90.o
main: CMakeFiles/main.dir/src/m_least_square.f90.o
main: CMakeFiles/main.dir/src/m_math.f90.o
main: CMakeFiles/main.dir/src/mod_test.f90.o
main: CMakeFiles/main.dir/src/mod_write.f90.o
main: CMakeFiles/main.dir/src/m_tests.f90.o
main: CMakeFiles/main.dir/build.make
main: ../third_party/liblapack.a
main: ../third_party/libblas.a
main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking Fortran executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: main

.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2 /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2 /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/build/CMakeFiles/main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend
