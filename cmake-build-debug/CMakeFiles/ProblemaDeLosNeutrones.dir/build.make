# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.28

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "E:\Clion\CLion 2024.1.1\bin\cmake\win\x64\bin\cmake.exe"

# The command to remove a file.
RM = "E:\Clion\CLion 2024.1.1\bin\cmake\win\x64\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/ProblemaDeLosNeutrones.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/ProblemaDeLosNeutrones.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/ProblemaDeLosNeutrones.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ProblemaDeLosNeutrones.dir/flags.make

CMakeFiles/ProblemaDeLosNeutrones.dir/main.c.obj: CMakeFiles/ProblemaDeLosNeutrones.dir/flags.make
CMakeFiles/ProblemaDeLosNeutrones.dir/main.c.obj: E:/AA\ TEC/5\ Semestre\ 1\ 2024/Lenguajes/Proyecto/ProblemaDeLosNeutrones/main.c
CMakeFiles/ProblemaDeLosNeutrones.dir/main.c.obj: CMakeFiles/ProblemaDeLosNeutrones.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/ProblemaDeLosNeutrones.dir/main.c.obj"
	"E:\Clion\CLion 2024.1.1\bin\mingw\bin\gcc.exe" $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/ProblemaDeLosNeutrones.dir/main.c.obj -MF CMakeFiles\ProblemaDeLosNeutrones.dir\main.c.obj.d -o CMakeFiles\ProblemaDeLosNeutrones.dir\main.c.obj -c "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\main.c"

CMakeFiles/ProblemaDeLosNeutrones.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/ProblemaDeLosNeutrones.dir/main.c.i"
	"E:\Clion\CLion 2024.1.1\bin\mingw\bin\gcc.exe" $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\main.c" > CMakeFiles\ProblemaDeLosNeutrones.dir\main.c.i

CMakeFiles/ProblemaDeLosNeutrones.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/ProblemaDeLosNeutrones.dir/main.c.s"
	"E:\Clion\CLion 2024.1.1\bin\mingw\bin\gcc.exe" $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\main.c" -o CMakeFiles\ProblemaDeLosNeutrones.dir\main.c.s

CMakeFiles/ProblemaDeLosNeutrones.dir/cJSON.c.obj: CMakeFiles/ProblemaDeLosNeutrones.dir/flags.make
CMakeFiles/ProblemaDeLosNeutrones.dir/cJSON.c.obj: E:/AA\ TEC/5\ Semestre\ 1\ 2024/Lenguajes/Proyecto/ProblemaDeLosNeutrones/cJSON.c
CMakeFiles/ProblemaDeLosNeutrones.dir/cJSON.c.obj: CMakeFiles/ProblemaDeLosNeutrones.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/ProblemaDeLosNeutrones.dir/cJSON.c.obj"
	"E:\Clion\CLion 2024.1.1\bin\mingw\bin\gcc.exe" $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/ProblemaDeLosNeutrones.dir/cJSON.c.obj -MF CMakeFiles\ProblemaDeLosNeutrones.dir\cJSON.c.obj.d -o CMakeFiles\ProblemaDeLosNeutrones.dir\cJSON.c.obj -c "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\cJSON.c"

CMakeFiles/ProblemaDeLosNeutrones.dir/cJSON.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/ProblemaDeLosNeutrones.dir/cJSON.c.i"
	"E:\Clion\CLion 2024.1.1\bin\mingw\bin\gcc.exe" $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\cJSON.c" > CMakeFiles\ProblemaDeLosNeutrones.dir\cJSON.c.i

CMakeFiles/ProblemaDeLosNeutrones.dir/cJSON.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/ProblemaDeLosNeutrones.dir/cJSON.c.s"
	"E:\Clion\CLion 2024.1.1\bin\mingw\bin\gcc.exe" $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\cJSON.c" -o CMakeFiles\ProblemaDeLosNeutrones.dir\cJSON.c.s

# Object files for target ProblemaDeLosNeutrones
ProblemaDeLosNeutrones_OBJECTS = \
"CMakeFiles/ProblemaDeLosNeutrones.dir/main.c.obj" \
"CMakeFiles/ProblemaDeLosNeutrones.dir/cJSON.c.obj"

# External object files for target ProblemaDeLosNeutrones
ProblemaDeLosNeutrones_EXTERNAL_OBJECTS =

ProblemaDeLosNeutrones.exe: CMakeFiles/ProblemaDeLosNeutrones.dir/main.c.obj
ProblemaDeLosNeutrones.exe: CMakeFiles/ProblemaDeLosNeutrones.dir/cJSON.c.obj
ProblemaDeLosNeutrones.exe: CMakeFiles/ProblemaDeLosNeutrones.dir/build.make
ProblemaDeLosNeutrones.exe: CMakeFiles/ProblemaDeLosNeutrones.dir/linkLibs.rsp
ProblemaDeLosNeutrones.exe: CMakeFiles/ProblemaDeLosNeutrones.dir/objects1.rsp
ProblemaDeLosNeutrones.exe: CMakeFiles/ProblemaDeLosNeutrones.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable ProblemaDeLosNeutrones.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\ProblemaDeLosNeutrones.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ProblemaDeLosNeutrones.dir/build: ProblemaDeLosNeutrones.exe
.PHONY : CMakeFiles/ProblemaDeLosNeutrones.dir/build

CMakeFiles/ProblemaDeLosNeutrones.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\ProblemaDeLosNeutrones.dir\cmake_clean.cmake
.PHONY : CMakeFiles/ProblemaDeLosNeutrones.dir/clean

CMakeFiles/ProblemaDeLosNeutrones.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones" "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones" "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\cmake-build-debug" "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\cmake-build-debug" "E:\AA TEC\5 Semestre 1 2024\Lenguajes\Proyecto\ProblemaDeLosNeutrones\cmake-build-debug\CMakeFiles\ProblemaDeLosNeutrones.dir\DependInfo.cmake" "--color=$(COLOR)"
.PHONY : CMakeFiles/ProblemaDeLosNeutrones.dir/depend

