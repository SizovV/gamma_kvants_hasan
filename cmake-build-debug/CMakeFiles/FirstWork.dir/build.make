# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.15

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2019.2.4\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2019.2.4\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\University\8Semestr\CalculationPhys\FirstWork

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\University\8Semestr\CalculationPhys\FirstWork\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/FirstWork.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/FirstWork.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/FirstWork.dir/flags.make

CMakeFiles/FirstWork.dir/main.cpp.obj: CMakeFiles/FirstWork.dir/flags.make
CMakeFiles/FirstWork.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\University\8Semestr\CalculationPhys\FirstWork\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/FirstWork.dir/main.cpp.obj"
	C:\PROGRA~1\mingw-w64\x86_64-8.1.0-win32-seh-rt_v6-rev0\mingw64\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\FirstWork.dir\main.cpp.obj -c D:\University\8Semestr\CalculationPhys\FirstWork\main.cpp

CMakeFiles/FirstWork.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FirstWork.dir/main.cpp.i"
	C:\PROGRA~1\mingw-w64\x86_64-8.1.0-win32-seh-rt_v6-rev0\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\University\8Semestr\CalculationPhys\FirstWork\main.cpp > CMakeFiles\FirstWork.dir\main.cpp.i

CMakeFiles/FirstWork.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FirstWork.dir/main.cpp.s"
	C:\PROGRA~1\mingw-w64\x86_64-8.1.0-win32-seh-rt_v6-rev0\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\University\8Semestr\CalculationPhys\FirstWork\main.cpp -o CMakeFiles\FirstWork.dir\main.cpp.s

# Object files for target FirstWork
FirstWork_OBJECTS = \
"CMakeFiles/FirstWork.dir/main.cpp.obj"

# External object files for target FirstWork
FirstWork_EXTERNAL_OBJECTS =

FirstWork.exe: CMakeFiles/FirstWork.dir/main.cpp.obj
FirstWork.exe: CMakeFiles/FirstWork.dir/build.make
FirstWork.exe: CMakeFiles/FirstWork.dir/linklibs.rsp
FirstWork.exe: CMakeFiles/FirstWork.dir/objects1.rsp
FirstWork.exe: CMakeFiles/FirstWork.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\University\8Semestr\CalculationPhys\FirstWork\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable FirstWork.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\FirstWork.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/FirstWork.dir/build: FirstWork.exe

.PHONY : CMakeFiles/FirstWork.dir/build

CMakeFiles/FirstWork.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\FirstWork.dir\cmake_clean.cmake
.PHONY : CMakeFiles/FirstWork.dir/clean

CMakeFiles/FirstWork.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\University\8Semestr\CalculationPhys\FirstWork D:\University\8Semestr\CalculationPhys\FirstWork D:\University\8Semestr\CalculationPhys\FirstWork\cmake-build-debug D:\University\8Semestr\CalculationPhys\FirstWork\cmake-build-debug D:\University\8Semestr\CalculationPhys\FirstWork\cmake-build-debug\CMakeFiles\FirstWork.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/FirstWork.dir/depend

