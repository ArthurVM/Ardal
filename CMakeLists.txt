cmake_minimum_required(VERSION 3.15)
project(Ardal LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Since pybind11 headers are vendored in the 'include' directory,
# we don't use find_package(pybind11). Instead, we find Python
# and add the include directory manually.
find_package(PythonLibs 3.8 REQUIRED)
find_package(OpenMP REQUIRED)

# Create the Python module using the standard CMake add_library command.
# This replaces the pybind11_add_module helper, which is not available
# when not using find_package(pybind11).
add_library(_ardal MODULE
    src/bindings.cpp
    src/BitMatrix.cpp
    src/RoaringMatrix.cpp
    src/HybridMatrix.cpp
)

# Add all necessary include directories to the target
target_include_directories(_ardal PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"  # For <pybind11/pybind11.h>
    "${CMAKE_CURRENT_SOURCE_DIR}/src"       # For your project's headers
    ${PYTHON_INCLUDE_DIRS}                  # For <Python.h>
)

# Add compile options from your setup.py
target_compile_options(_ardal PRIVATE
    -O3
    -march=native
    -ffast-math
)

target_link_libraries(_ardal PRIVATE ${PYTHON_LIBRARIES} OpenMP::OpenMP_CXX)

# Set properties to create a valid Python module file (e.g., _ardal.cpython-38-x86_64-linux-gnu.so)
set_target_properties(_ardal PROPERTIES PREFIX "" SUFFIX "${PYTHON_EXTENSION_SUFFIX}")

# Set the output directory to match your package structure for easy testing
set_target_properties(_ardal PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/ardal")