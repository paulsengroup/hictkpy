# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...3.30)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

project(
  hictkpy
  LANGUAGES
    C
    CXX
  HOMEPAGE_URL https://github.com/paulsengroup/hictkpy
  DESCRIPTION "Python bindings for hictk."
)

include(FetchContent)
FetchContent_Declare(
  _hictkpy_project_options
  URL
    "${CMAKE_CURRENT_SOURCE_DIR}/external/project_options-v0.36.6.tar.xz"
  URL_HASH "SHA256=b4df2a4107847248d4711ee48477cdf93c1a1c7d51a931937759eb9ffba67af8"
  EXCLUDE_FROM_ALL
  SYSTEM
)

FetchContent_Declare(
  hictk
  URL
    "${CMAKE_CURRENT_SOURCE_DIR}/external/hictk-v2.1.2.tar.xz"
  URL_HASH "SHA256=82cb76d595c0760eb1dffbcd09b98141797bd231aa269aef76e1f9b63f8e0ee5"
  EXCLUDE_FROM_ALL
  OVERRIDE_FIND_PACKAGE
  SYSTEM
)

FetchContent_Declare(
  nanobind
  URL
    "${CMAKE_CURRENT_SOURCE_DIR}/external/nanobind-v2.9.2.tar.xz"
  URL_HASH "SHA256=1310c9e6c8801ba18dd9dbc9699e3ac823d4f36a9dcb8e06e0ab6855fbad0d81"
  EXCLUDE_FROM_ALL
  OVERRIDE_FIND_PACKAGE
  SYSTEM
)

# Tweak hictk
set(HICTK_BUILD_BENCHMARKS OFF)
set(HICTK_BUILD_EXAMPLES OFF)
set(HICTK_BUILD_TOOLS OFF)
set(HICTK_ENABLE_FUZZY_TESTING OFF)
set(HICTK_ENABLE_TESTING OFF)
set(HICTK_INSTALL OFF)
set(HICTK_WITH_EIGEN ON)
set(HICTK_ENABLE_GIT_VERSION_TRACKING OFF)
find_package(Eigen3 5 REQUIRED QUIET) # Versions <5 are not supported (and will likely result in segfaults)

# We are finding hictk here to avoid name collisions with the ENABLE_* variables defined in the coming lines
find_package(hictk REQUIRED)

get_property(BUILDING_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILDING_MULTI_CONFIG)
  if(NOT CMAKE_BUILD_TYPE)
    # Make sure that all supported configuration types have their associated conan packages available. You can reduce
    # this list to only the configuration types you use, but only if one is not forced-set on the command line for VS
    message(TRACE "Setting up multi-config build types")
    set(
      CMAKE_CONFIGURATION_TYPES
      Debug
      Release
      RelWithDebInfo
      CACHE STRING
      "Enabled build types"
      FORCE
    )
  else()
    message(TRACE "User chose a specific build type, so we are using that")
    set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} CACHE STRING "Enabled build types" FORCE)
  endif()
endif()

FetchContent_MakeAvailable(_hictkpy_project_options)

include("${_hictkpy_project_options_SOURCE_DIR}/src/DynamicProjectOptions.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerWarnings.cmake")

set(ENABLE_CACHE_DEFAULT ON)
set(ENABLE_COMPILE_COMMANDS_SYMLINK_DEFAULT OFF)
set(ENABLE_CONAN_DEFAULT OFF)
set(ENABLE_CPPCHECK_DEFAULT OFF)
set(ENABLE_DOXYGEN_DEFAULT OFF)
set(ENABLE_INTERPROCEDURAL_OPTIMIZATION_DEFAULT ON)
set(ENABLE_NATIVE_OPTIMIZATION_DEFAULT OFF)
set(ENABLE_PCH_DEFAULT OFF)

set(ENABLE_SANITIZER_ADDRESS_DEFAULT OFF)
set(ENABLE_SANITIZER_LEAK_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_COMPARE_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_SUBTRACT_DEFAULT OFF)
set(ENABLE_SANITIZER_UNDEFINED_DEFAULT OFF)

dynamic_project_options(
        PREFIX
        hictkpy
        # CLANG_WARNINGS # Override the defaults for the CLANG warnings GCC_WARNINGS   # Override the defaults for the GCC
        CPPCHECK_OPTIONS
        --enable=performance,portability,style,warning
        --inline-suppr
        # We cannot act on a bug/missing feature of cppcheck
        --suppress=internalAstError
        # if a file does not have an internalAstError, we get an unmatchedSuppression error
        --suppress=unmatchedSuppression
        --suppress=passedByValue
        --inconclusive
        MSVC_WARNINGS
        "${MSVC_WARNINGS}"
        CLANG_WARNINGS
        "${CLANG_WARNINGS}"
        GCC_WARNINGS
        "${GCC_WARNINGS}"
        CUDA_WARNINGS
        "${CUDA_WARNINGS}"
)

target_compile_features(hictkpy_project_options INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")

option(HICTKPY_BUILD_PORTABLE_WHEELS_WIN "Build portable Windows wheels" OFF)

# Tweak fmt
target_compile_definitions(
  hictkpy_project_options
  INTERFACE
    FMT_HEADER_ONLY
    FMT_ENFORCE_COMPILE_STRING
)
# Tweak spdlog
target_compile_definitions(hictkpy_project_options INTERFACE SPDLOG_FMT_EXTERNAL)

# Tweak hictk
target_compile_definitions(
  hictkpy_project_options
  INTERFACE
    HICTK_WITH_ARROW
    HICTK_WITH_EIGEN
)

if(WIN32)
  target_compile_definitions(
    hictkpy_project_options
    INTERFACE
      NOMINMAX
      _CRT_SECURE_NO_WARNINGS
      _DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR # https://github.com/gabime/spdlog/issues/3212
  )

  if(HICTKPY_BUILD_PORTABLE_WHEELS_WIN)
    target_compile_definitions(
      hictkpy_project_options
      INTERFACE
        _WIN32_WINNT=_WIN32_WINNT_WIN10
        /FI<sdkddkver.h> # required to use _WIN32_WINNT_WIN10
    )
  endif()
endif()

add_subdirectory(src)
