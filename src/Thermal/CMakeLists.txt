cmake_minimum_required(VERSION 3.14)
project(ThermalModel LANGUAGES CXX)

# Gather sources (e.g., thermal.cpp and any related files)
file(GLOB THERMAL_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Build the Thermal model library (can be STATIC or SHARED)
add_library(Thermal_model ${THERMAL_SOURCES})

# Set include directories for Thermal model
target_include_directories(Thermal_model PUBLIC ${CMAKE_SOURCE_DIR}/include/Thermal)

# Link common dependencies (e.g., Eigen)
target_link_libraries(Thermal_model PRIVATE Eigen3::Eigen)
