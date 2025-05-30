cmake_minimum_required(VERSION 3.16)
project(shasta)

# This is the top level cmake file for a Shasta build. 

# Before starting a build,
# run shasta/scripts/InstallPrerequisites-Ubuntu.sh

# Use "cmake -D" to modify the options defined below 
# in each "option" line. For example, to turn off
# building the dynamic executable use
# cmake -DBUILD_DYNAMIC_EXECUTABLE=OFF



# Figure out if we are on Intel (as opposed to ARM).
set(X86_64 OFF)
if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
    set(X86_64 ON)
endif()



# Summarize the environment.
message(STATUS "Processor architecture is " ${CMAKE_HOST_SYSTEM_PROCESSOR})
message(STATUS "CMAKE_SYSTEM_NAME is " ${CMAKE_SYSTEM_NAME})



# Decide what we want to build.
option(BUILD_STATIC_LIBRARY "Build the static library." ON)
option(BUILD_STATIC_EXECUTABLE "Build the static executable." ON)
option(BUILD_DYNAMIC_LIBRARY "Build the shared library." ON)
option(BUILD_DYNAMIC_EXECUTABLE "Build the dynamic executable." OFF)



# To build the static executable, we also need the static library.
if(BUILD_STATIC_EXECUTABLE AND NOT BUILD_STATIC_LIBRARY)
    set(BUILD_STATIC_LIBRARY ON)
    message(STATUS "Turned on BUILD_STATIC_LIBRARY because it is needed for BUILD_STATIC_EXECUTABLE.")
endif(BUILD_STATIC_EXECUTABLE AND NOT BUILD_STATIC_LIBRARY)



# To build the dynamic executable, we also need the dynamic library.
if(BUILD_DYNAMIC_EXECUTABLE AND NOT BUILD_DYNAMIC_LIBRARY)
    set(BUILD_DYNAMIC_LIBRARY ON)
    message(STATUS "Turned on BUILD_DYNAMIC_LIBRARY because it is needed for BUILD_DYNAMIC_EXECUTABLE.")
endif(BUILD_DYNAMIC_EXECUTABLE AND NOT BUILD_DYNAMIC_LIBRARY)



# Write out what we will actually build.
message(STATUS "BUILD_STATIC_LIBRARY is " ${BUILD_STATIC_LIBRARY})
message(STATUS "BUILD_STATIC_EXECUTABLE is " ${BUILD_STATIC_EXECUTABLE})
message(STATUS "BUILD_DYNAMIC_LIBRARY is " ${BUILD_DYNAMIC_LIBRARY})
message(STATUS "BUILD_DYNAMIC_EXECUTABLE is " ${BUILD_DYNAMIC_EXECUTABLE})



# Option to build with -march=native.
option(BUILD_NATIVE "Build with -march=native." OFF)
message(STATUS "BUILD_NATIVE is " ${BUILD_NATIVE})



# Option to request a debug build.
option(BUILD_DEBUG "Make a debuggable build." OFF)
message(STATUS "BUILD_DEBUG is " ${BUILD_DEBUG})



# The BUILD_ID can be specified to identify the build
# This is normally used only when building a new GitHub release,
# in which case we use the following option when running Cmake:
# -DBUILD_ID="Shasta Release X.Y.Z"
# BUILD_ID should not contain a comma, otherwise the C preprocessor does not handle it correctly.
if(NOT DEFINED BUILD_ID)
    set(BUILD_ID "Shasta development build. This is not a released version.")
endif(NOT DEFINED BUILD_ID)
message(STATUS "BUILD_ID is: " ${BUILD_ID})


# Option to request a build with long markers.
option(SHASTA_LONG_MARKERS "Make a build with SHASTA_LONG_MARKERS set." OFF)
message(STATUS "SHASTA_LONG_MARKERS is " ${SHASTA_LONG_MARKERS})



# Add the subdirectories we need.
if(BUILD_STATIC_LIBRARY)
    add_subdirectory(staticLibrary)
endif(BUILD_STATIC_LIBRARY)

if(BUILD_STATIC_EXECUTABLE)
    add_subdirectory(staticExecutable)
endif(BUILD_STATIC_EXECUTABLE)

if(BUILD_DYNAMIC_LIBRARY)
    add_subdirectory(dynamicLibrary)
endif(BUILD_DYNAMIC_LIBRARY)

if(BUILD_DYNAMIC_EXECUTABLE)
    add_subdirectory(dynamicExecutable)
endif(BUILD_DYNAMIC_EXECUTABLE)



# Install to the shasta-install directory.
set(CMAKE_INSTALL_PREFIX .)

# Install the scripts.
file(GLOB SCRIPTS scripts/*.py scripts/*.sh)
install(PROGRAMS ${SCRIPTS} DESTINATION shasta-install/bin)

# Install the configuration files.
install(DIRECTORY conf DESTINATION shasta-install USE_SOURCE_PERMISSIONS)

# Install the docs directory.
install(DIRECTORY docs DESTINATION shasta-install)

# The targets built in each subdirectory are
# installed by the cmake file of each subdirectory.




