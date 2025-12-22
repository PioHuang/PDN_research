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
    src/VoltSpotVirtualGrid.cpp
    src/PDNLoadMapper.cpp
)

# Example executable
add_executable(pdn_example src/examples/example.cpp ${CORE_SOURCES})

# VoltSpot virtual grid extraction example
add_executable(voltspot_grid_example src/examples/voltspot_grid_example.cpp ${CORE_SOURCES})

# PDN load mapper example
add_executable(pdn_load_example src/examples/pdn_load_example.cpp ${CORE_SOURCES})

# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pdn_example PRIVATE -Wall -Wextra)
    target_compile_options(voltspot_grid_example PRIVATE -Wall -Wextra)
    target_compile_options(pdn_load_example PRIVATE -Wall -Wextra)
endif()

# Installation (optional)
install(TARGETS pdn_example DESTINATION bin)
install(TARGETS voltspot_grid_example DESTINATION bin)
install(TARGETS pdn_load_example DESTINATION bin)
