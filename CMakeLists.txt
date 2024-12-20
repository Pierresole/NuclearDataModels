# Top-level CMakeLists.txt
cmake_minimum_required(VERSION 3.14)
project(pyrat LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # Specify Python 3.10 paths explicitly
# set(Python3_EXECUTABLE /usr/local/bin/python3.10)
# set(Python3_INCLUDE_DIR /usr/local/include/python3.10)
# set(Python3_LIBRARY /usr/local/lib/libpython3.10.so)
# set(PYTHON_EXECUTABLE /usr/local/bin/python3.10 CACHE FILEPATH "Python interpreter")
# set(PYTHON_LIBRARY /usr/local/lib/libpython3.10.so CACHE FILEPATH "Python library")
# set(PYTHON_INCLUDE_DIR /usr/local/include/python3.10 CACHE FILEPATH "Python include directory")

# # Prevent CMake from finding other Python versions
# unset(Python_ADDITIONAL_VERSIONS)
# set(Python_ADDITIONAL_VERSIONS "3.10")

find_package(PythonLibs 3.10 REQUIRED)
find_package(PythonInterp 3.10 REQUIRED)

find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter Development)

find_package(pybind11 REQUIRED)

add_subdirectory(cpp)
