cmake_minimum_required(VERSION 3.16)
project(shastaStaticExecutable)

# C++ dialect.
add_definitions(-std=c++20)

# Compilation warnings.
add_definitions(-Wall -Wconversion -Wno-unused-result -Wno-trigraphs -Wno-psabi -Wunused-parameter)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O0)
else(BUILD_DEBUG)
    add_definitions(-g0)
    add_definitions(-O3)
    # NDEBUG is required to turn off SeqAn debug code.
    add_definitions(-DNDEBUG)
endif(BUILD_DEBUG)

# 16-byte compare and swap.
# This is recommended for dset64.hpp/dset64-gccAtomic.hpp".
# It's available only on x86 architectures.
if(X86_64)
    add_definitions(-mcx16)
endif(X86_64)

# Native build.
if(BUILD_NATIVE)
    add_definitions(-march=native)
endif(BUILD_NATIVE)

# Build id.
add_definitions(-DBUILD_ID=${BUILD_ID})

# SHASTA_LONG_MARKERS.
if(SHASTA_LONG_MARKERS)
    add_definitions(-DSHASTA_LONG_MARKERS)
endif(SHASTA_LONG_MARKERS)

# Definitions needed to eliminate dependency on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# This is needed to avoid some Boost warnings about deprecated messages.
# We canot fix this in Shasta as it is a Boost problem.
# It can be removed if Boost is fixed.
add_definitions(-DBOOST_ALLOW_DEPRECATED_HEADERS)

# Source files
file(GLOB SOURCES ../srcMain/*.cpp)

# Include directory.
include_directories(../src)

# Define our executable.
add_executable(shastaStaticExecutable ${SOURCES})
set_target_properties(shastaStaticExecutable PROPERTIES OUTPUT_NAME "shasta")

# Request a static executable.
set_target_properties(shastaStaticExecutable PROPERTIES LINK_FLAGS "-static" )

# Libraries to link with.
# For arcane reasons, statically linking with the pthread
# library on Linux requires "--whole-archive".
if(X86_64)
    target_link_libraries(
        shastaStaticExecutable
        shastaStaticLibrary
        atomic boost_system boost_program_options boost_chrono boost_serialization  spoa png z
        lapack blas gfortran quadmath
        -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
else(X86_64)
    target_link_libraries(
        shastaStaticExecutable
        shastaStaticLibrary
        atomic boost_system boost_program_options boost_chrono boost_serialization spoa png z
        lapack blas gfortran
        -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
endif(X86_64) 
  
# The static executable goes to the bin directory.
install(TARGETS shastaStaticExecutable DESTINATION shasta-install/bin)


