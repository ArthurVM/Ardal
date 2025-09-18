cmake_minimum_required(VERSION 3.15)
project(Ardal LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ARDAL_USE_OPENMP "Enable OpenMP" ON)
option(ARDAL_USE_SIMD  "Enable SIMD (AVX2/AVX-512 if available)" ON)

find_package(Python REQUIRED COMPONENTS Interpreter Development)

set(ARDAL_CORE_SRCS
  src/bindings/bindings_core.cpp
  src/utils/Logger.cpp
  src/core/BitMatrix.cpp
  src/core/RoaringMatrix.cpp
  src/core/HybridMatrix.cpp
  src/roaring/roaring.c
)
set(ARDAL_PACK_SRCS src/bindings/bindings_pack.cpp)

## expose vendored pybind11 headers
add_library(pybind11_headers INTERFACE)
target_include_directories(pybind11_headers INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(_ardal MODULE ${ARDAL_CORE_SRCS})
add_library(_ardal_pack MODULE ${ARDAL_PACK_SRCS})

foreach(tgt IN ITEMS _ardal _ardal_pack)
  target_include_directories(${tgt} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/roaring
  )
  target_link_libraries(${tgt} PRIVATE pybind11_headers Python::Python)
  target_compile_features(${tgt} PRIVATE cxx_std_17)
endforeach()

execute_process(
  COMMAND ${Python_EXECUTABLE} - << "PYCODE"
import sysconfig, sys
print(sysconfig.get_config_var("EXT_SUFFIX") or "")
PYCODE
  OUTPUT_VARIABLE PY_EXT_SUFFIX OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT PY_EXT_SUFFIX)
  ## fallback: .so on Unix, .pyd on Windows
  if(WIN32)
    set(PY_EXT_SUFFIX ".pyd")
  else()
    set(PY_EXT_SUFFIX ".so")
  endif()
endif()
foreach(tgt IN ITEMS _ardal _ardal_pack)
  set_target_properties(${tgt} PROPERTIES
    PREFIX ""                # drop 'lib'
    SUFFIX "${PY_EXT_SUFFIX}"
  )
endforeach()

## OpenMP / SIMD same as above
if(ARDAL_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(_ardal PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(_ardal_pack PRIVATE OpenMP::OpenMP_CXX)
endif()
include(CheckCXXCompilerFlag)
if(ARDAL_USE_SIMD AND NOT MSVC)
  check_cxx_compiler_flag("-mavx2" HAS_AVX2)
  if(HAS_AVX2)
    target_compile_options(_ardal PRIVATE -mavx2 -ffast-math)
    target_compile_options(_ardal_pack PRIVATE -mavx2 -ffast-math)
  endif()
elseif(ARDAL_USE_SIMD AND MSVC)
  target_compile_options(_ardal PRIVATE /arch:AVX2)
  target_compile_options(_ardal_pack PRIVATE /arch:AVX2)
endif()

install(TARGETS _ardal _ardal_pack LIBRARY DESTINATION ardal)
