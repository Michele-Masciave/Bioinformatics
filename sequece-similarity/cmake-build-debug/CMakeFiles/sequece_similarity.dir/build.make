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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/michele/Desktop/sequece-similarity

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/michele/Desktop/sequece-similarity/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/sequece_similarity.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/sequece_similarity.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/sequece_similarity.dir/flags.make

CMakeFiles/sequece_similarity.dir/main.cpp.o: CMakeFiles/sequece_similarity.dir/flags.make
CMakeFiles/sequece_similarity.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/michele/Desktop/sequece-similarity/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/sequece_similarity.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/sequece_similarity.dir/main.cpp.o -c /Users/michele/Desktop/sequece-similarity/main.cpp

CMakeFiles/sequece_similarity.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sequece_similarity.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/michele/Desktop/sequece-similarity/main.cpp > CMakeFiles/sequece_similarity.dir/main.cpp.i

CMakeFiles/sequece_similarity.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sequece_similarity.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/michele/Desktop/sequece-similarity/main.cpp -o CMakeFiles/sequece_similarity.dir/main.cpp.s

# Object files for target sequece_similarity
sequece_similarity_OBJECTS = \
"CMakeFiles/sequece_similarity.dir/main.cpp.o"

# External object files for target sequece_similarity
sequece_similarity_EXTERNAL_OBJECTS =

sequece_similarity: CMakeFiles/sequece_similarity.dir/main.cpp.o
sequece_similarity: CMakeFiles/sequece_similarity.dir/build.make
sequece_similarity: CMakeFiles/sequece_similarity.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/michele/Desktop/sequece-similarity/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable sequece_similarity"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sequece_similarity.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/sequece_similarity.dir/build: sequece_similarity

.PHONY : CMakeFiles/sequece_similarity.dir/build

CMakeFiles/sequece_similarity.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/sequece_similarity.dir/cmake_clean.cmake
.PHONY : CMakeFiles/sequece_similarity.dir/clean

CMakeFiles/sequece_similarity.dir/depend:
	cd /Users/michele/Desktop/sequece-similarity/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/michele/Desktop/sequece-similarity /Users/michele/Desktop/sequece-similarity /Users/michele/Desktop/sequece-similarity/cmake-build-debug /Users/michele/Desktop/sequece-similarity/cmake-build-debug /Users/michele/Desktop/sequece-similarity/cmake-build-debug/CMakeFiles/sequece_similarity.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/sequece_similarity.dir/depend
