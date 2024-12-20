cmake_minimum_required(VERSION 3.14)
project(pyrat)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(CHECK_START "Fetching Eigen3")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

include(FetchContent)
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
# note: To disable eigen tests,
# you should put this code in a add_subdirectory to avoid to change
# BUILD_TESTING for your own project too since variables are directory
# scoped
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_MPL2_ONLY ON)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_BUILD_DOC OFF)
FetchContent_MakeAvailable(Eigen)

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "fetched")

# Collect all source files from src directory
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# Add the -fPIC flag to the compilaBtion options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# Add the pyrat_lib library (core functionality)
add_library(pyrat_lib ${SOURCES})

# Include headers directories
target_include_directories(pyrat_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link libraries if needed (e.g., Python)
target_link_libraries(pyrat_lib 
  PRIVATE 
  # ENDFtk
  Python3::Python 
  Eigen3::Eigen)

# Find and include pybind11
find_package(pybind11 REQUIRED)

# Add the Pybind11 bindings as a separate module
add_library(pyrat MODULE pybind_wrapper.cpp)

# Ensure correct naming for Python import
set_target_properties(pyrat PROPERTIES
    PREFIX ""
    SUFFIX ".so"  # Correct extension for Python modules
)

# Link the bindings with the core library
target_link_libraries(pyrat PRIVATE pyrat_lib pybind11::module Eigen3::Eigen)

# Ensure that the module is placed in the correct output directory
set_target_properties(pyrat PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python  # Output to a 'python' folder within your build directory
)