cmake_minimum_required(VERSION 3.20)
project(NuclearDataModels LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
message(STATUS "Using Python headers from: ${Python3_INCLUDE_DIRS}")
message(STATUS "Using Python library from: ${Python3_LIBRARIES}")

# Fetch Eigen (if needed)
include(FetchContent)
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
  GIT_SHALLOW TRUE
)
set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_MPL2_ONLY ON)
FetchContent_MakeAvailable(Eigen)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.13.6
)
FetchContent_MakeAvailable(pybind11)


# Common include directory
include_directories(${CMAKE_SOURCE_DIR}/include)

# Add subdirectories for each model library:
add_subdirectory(src/RMatrix)
# add_subdirectory(src/Thermal)
# add_subdirectory(src/Optical)

# Add Python bindings subdirectories:
add_subdirectory(bindings/python/RMatrix)
# add_subdirectory(bindings/python/Thermal)
# add_subdirectory(bindings/python/Optical)

# Custom clean command that preserves _deps directory
add_custom_target(clean-except-deps
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaning build directory except _deps/"
    COMMAND find ${CMAKE_BINARY_DIR} -mindepth 1 -maxdepth 1 -not -name "_deps" -type d -exec ${CMAKE_COMMAND} -E remove_directory {} +
    COMMAND find ${CMAKE_BINARY_DIR} -mindepth 1 -maxdepth 1 -not -name "_deps" -type f -exec ${CMAKE_COMMAND} -E remove {} +
    COMMENT "Cleaning build directory while preserving _deps/"
)