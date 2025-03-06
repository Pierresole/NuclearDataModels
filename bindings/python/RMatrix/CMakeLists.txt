cmake_minimum_required(VERSION 3.14)
project(pyRMatrix LANGUAGES CXX)

find_package(pybind11 REQUIRED)

# Create the Python module from your binding source (bindings.cpp)
add_library(pyRMatrix MODULE bindings.cpp)

# Set properties so that the module is importable from Python:
set_target_properties(pyRMatrix PROPERTIES
    PREFIX ""         # No 'lib' prefix
    SUFFIX ".so"      # Use .so (or .pyd on Windows)
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
)

# Link the module with the Thermal model library and pybind11
target_link_libraries(pyRMatrix PRIVATE RMatrix_model pybind11::module Eigen3::Eigen)

# Optionally include the Thermal headers (if needed)
target_include_directories(pyRMatrix PRIVATE ${CMAKE_SOURCE_DIR}/include/RMatrix ${eigen_SOURCE_DIR})

# Fix site-packages detection
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c "import site; print(site.getsitepackages()[0])"
  OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Python site-packages directory: ${PYTHON_SITE_PACKAGES}")

# Install target to the site-packages directory
install(TARGETS pyRMatrix 
        LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}"
        COMPONENT python)