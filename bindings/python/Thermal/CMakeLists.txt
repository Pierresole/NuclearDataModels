cmake_minimum_required(VERSION 3.14)
project(pyThermal LANGUAGES CXX)

find_package(pybind11 REQUIRED)

# Create the Python module from your binding source (bindings.cpp)
add_library(pyThermal MODULE bindings.cpp)

# Set properties so that the module is importable from Python:
set_target_properties(pyThermal PROPERTIES
    PREFIX ""         # No 'lib' prefix
    SUFFIX ".so"      # Use .so (or .pyd on Windows)
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
)

# Link the module with the Thermal model library and pybind11
target_link_libraries(pyThermal PRIVATE Thermal_model pybind11::module)

# Optionally include the Thermal headers (if needed)
target_include_directories(pyThermal PRIVATE ${CMAKE_SOURCE_DIR}/include/Thermal)
