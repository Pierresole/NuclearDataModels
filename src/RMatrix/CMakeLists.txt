cmake_minimum_required(VERSION 3.14)
project(RMatrixModel LANGUAGES CXX)

# Add the -fPIC flag to the compilaBtion options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# Gather sources (e.g., RMatrix.cpp and any related files)
file(GLOB RMATRIX_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Build the RMatrix model library (can be STATIC or SHARED)
add_library(RMatrix_model ${RMATRIX_SOURCES})

# Set include directories for RMatrix model
target_include_directories(RMatrix_model 
    PUBLIC 
        ${CMAKE_SOURCE_DIR}/include/RMatrix
    PRIVATE
        ${eigen_SOURCE_DIR}
)

# Link Eigen properly using modern CMake
target_link_libraries(RMatrix_model PUBLIC Eigen3::Eigen)
