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

# Find Eigen3
find_package(Eigen3 QUIET)

if(Eigen3_FOUND AND EIGEN3_INCLUDE_DIR)
    include_directories(${EIGEN3_INCLUDE_DIR})
    message(STATUS "Found Eigen3 via find_package: ${EIGEN3_INCLUDE_DIR}")
else()
    set(EIGEN3_POSSIBLE_PATHS
        /opt/homebrew/opt/eigen/include/eigen3
        /opt/homebrew/include/eigen3
        /usr/local/opt/eigen/include/eigen3
        /usr/local/include/eigen3
        /usr/include/eigen3
        ${CMAKE_SOURCE_DIR}/third_party/eigen
    )
    
    find_path(EIGEN3_INCLUDE_DIR 
        NAMES Eigen/Dense
        PATHS ${EIGEN3_POSSIBLE_PATHS}
    )
    
    if(EIGEN3_INCLUDE_DIR)
        include_directories(${EIGEN3_INCLUDE_DIR})
        message(STATUS "Found Eigen3 at: ${EIGEN3_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR "Eigen3 not found. Please install Eigen3:")
        message(STATUS "  macOS: brew install eigen")
        message(STATUS "  Ubuntu: sudo apt-get install libeigen3-dev")
        message(STATUS "  Or set EIGEN3_INCLUDE_DIR manually")
        message(STATUS "  Searched in: ${EIGEN3_POSSIBLE_PATHS}")
    endif()
endif()

# Core library sources
set(CORE_SOURCES
    src/PixelModel.cpp
    src/PDNNetwork.cpp
    src/VoltSpotVirtualGrid.cpp
    src/PDNLoadMapper.cpp
    src/PDNSolver.cpp
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
