cmake_minimum_required(VERSION 3.15)
project(PDN_Model VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Core library sources
set(CORE_SOURCES
    src/PixelModel.cpp
    src/PDNNetwork.cpp
)

# Example executable
add_executable(pdn_example src/example.cpp ${CORE_SOURCES})

# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pdn_example PRIVATE -Wall -Wextra)
endif()

# Installation (optional)
install(TARGETS pdn_example DESTINATION bin)

