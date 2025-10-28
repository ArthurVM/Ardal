cmake_minimum_required(VERSION 3.15)
project(Ardal LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ARDAL_USE_OPENMP "Enable OpenMP" ON)
option(ARDAL_USE_SIMD  "Enable SIMD (x86-64-v3/AVX2 if available)" ON)
option(ARDAL_MARCH_NATIVE "Use -march=native (developer builds only)" OFF)
option(ARDAL_ENABLE_LTO "Enable link-time optimization" OFF)

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

## vendored pybind11 headers
add_library(pybind11_headers INTERFACE)
target_include_directories(pybind11_headers INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(ardal MODULE ${ARDAL_CORE_SRCS})
add_library(ardal_pack MODULE ${ARDAL_PACK_SRCS})

foreach(tgt IN ITEMS ardal ardal_pack)
  target_include_directories(${tgt} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/roaring
    ${CMAKE_CURRENT_SOURCE_DIR}/include
  )
  target_link_libraries(${tgt} PRIVATE pybind11_headers Python::Python)
  target_compile_features(${tgt} PRIVATE cxx_std_17)
endforeach()

## compute Python extension suffix
execute_process(
  COMMAND ${Python_EXECUTABLE} - << "PYCODE"
import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX") or "")
PYCODE
  OUTPUT_VARIABLE PY_EXT_SUFFIX OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT PY_EXT_SUFFIX)
  if(WIN32)
    set(PY_EXT_SUFFIX ".pyd")
  else()
    set(PY_EXT_SUFFIX ".so")
  endif()
endif()
foreach(tgt IN ITEMS ardal ardal_pack)
  set_target_properties(${tgt} PROPERTIES PREFIX "" SUFFIX "${PY_EXT_SUFFIX}")
endforeach()

# ------------------------------
# OpenMP
# ------------------------------
if(ARDAL_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  foreach(tgt IN ITEMS ardal ardal_pack)
    target_link_libraries(${tgt} PRIVATE OpenMP::OpenMP_CXX OpenMP::OpenMP_C)
  endforeach()
endif()

# ------------------------------
# SIMD / CPU features
# ------------------------------
include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)

## optional AVX-512 VPOPCNTDQ (off by default)
option(ARDAL_ENABLE_AVX512_VPOPCNTDQ "Enable AVX-512 VPOPCNTDQ kernels" OFF)

if(ARDAL_USE_SIMD)
  if(MSVC)
    foreach(tgt IN ITEMS ardal ardal_pack)
      target_compile_options(${tgt} PRIVATE /arch:AVX2)
      target_compile_definitions(${tgt} PUBLIC ARDAL_HAVE_AVX2_KERNELS=1)
      target_compile_definitions(${tgt} PUBLIC ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS=0)
    endforeach()
  else()
    ## prefer x86-64-v3 if supported
    ## else fall back to -mavx2 flags
    check_cxx_compiler_flag("-march=x86-64-v3" HAS_X64V3_CXX)
    check_c_compiler_flag("-march=x86-64-v3"  HAS_X64V3_C)

    if(HAS_X64V3_CXX AND HAS_X64V3_C)
      foreach(tgt IN ITEMS ardal ardal_pack)
        target_compile_options(${tgt} PRIVATE
          $<$<COMPILE_LANGUAGE:CXX>:-march=x86-64-v3 -ffast-math>
          $<$<COMPILE_LANGUAGE:C>:-march=x86-64-v3>
        )
        target_compile_definitions(${tgt} PUBLIC ARDAL_HAVE_AVX2_KERNELS=1)
      endforeach()
    else()
      foreach(tgt IN ITEMS ardal ardal_pack)
        target_compile_options(${tgt} PRIVATE
          $<$<COMPILE_LANGUAGE:CXX>:-mavx2 -mbmi -mbmi2 -mpopcnt -mfma -mlzcnt -ffast-math>
          $<$<COMPILE_LANGUAGE:C>:-mavx2 -mbmi -mbmi2 -mpopcnt -mfma -mlzcnt>
        )
        target_compile_definitions(${tgt} PUBLIC ARDAL_HAVE_AVX2_KERNELS=1)
      endforeach()
    endif()

    ## optional AVX-512 VPOPCNTDQ (explicit opt-in)
    if(ARDAL_ENABLE_AVX512_VPOPCNTDQ)
      foreach(tgt IN ITEMS ardal ardal_pack)
        target_compile_options(${tgt} PRIVATE
          $<$<COMPILE_LANGUAGE:CXX>:-mavx512f -mavx512bw -mavx512vl -mavx512vpopcntdq>
          $<$<COMPILE_LANGUAGE:C>:-mavx512f -mavx512bw -mavx512vl -mavx512vpopcntdq>
        )
        target_compile_definitions(${tgt} PUBLIC ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS=1)
      endforeach()
    else()
      foreach(tgt IN ITEMS ardal ardal_pack)
        target_compile_definitions(${tgt} PUBLIC ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS=0)
      endforeach()
    endif()
  endif()

  ## optional: builds with -march=native
  if(ARDAL_MARCH_NATIVE AND NOT MSVC)
    foreach(tgt IN ITEMS ardal ardal_pack)
      target_compile_options(${tgt} PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-march=native>
        $<$<COMPILE_LANGUAGE:C>:-march=native>
      )
    endforeach()
  endif()

else()
  ## SIMD disabled: macros explicitly false
  foreach(tgt IN ITEMS ardal ardal_pack)
    target_compile_definitions(${tgt} PUBLIC
      ARDAL_HAVE_AVX2_KERNELS=0
      ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS=0
    )
  endforeach()
endif()

# ------------------------------
# LTO
# ------------------------------
if(ARDAL_ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT has_ipo OUTPUT ipo_msg)
  if(has_ipo)
    set_property(TARGET ardal ardal_pack PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(WARNING "IPO/LTO not supported: ${ipo_msg}")
  endif()
endif()

install(TARGETS ardal ardal_pack LIBRARY DESTINATION ardal)