# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...4.1)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

project(
  hictkpy
  LANGUAGES
    C
    CXX
  HOMEPAGE_URL https://github.com/paulsengroup/hictkpy
  DESCRIPTION "Python bindings for hictk."
)

if(NOT DEFINED CMAKE_FIND_PACKAGE_PREFER_CONFIG)
  set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)
endif()

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
    "${CMAKE_CURRENT_SOURCE_DIR}/external/hictk-v2.2.0.tar.xz"
  URL_HASH "SHA256=8341a22a6f6115863549170da72cfe5282c0bcc206cbd872831e60925e3b1eb6"
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

if(DEFINED CMAKE_BUILD_TYPE)
  set(CONAN_INSTALL_ARGS "${CONAN_INSTALL_ARGS};--settings=build_type=${CMAKE_BUILD_TYPE}")
endif()

if(DEFINED CMAKE_CXX_STANDARD)
  set(CONAN_INSTALL_ARGS "${CONAN_INSTALL_ARGS};--settings:all=compiler.cppstd=${CMAKE_CXX_STANDARD}")
else()
  set(CONAN_INSTALL_ARGS "${CONAN_INSTALL_ARGS};--settings:all=compiler.cppstd=17")
endif()

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

# Try to look for ArrowCompute explicitly to find the packages when using e.g., conda.
# See https://github.com/paulsengroup/hictk/pull/477
find_package(Arrow QUIET)
if(NOT TARGET ArrowCompute::arrow_compute_shared AND NOT TARGET ArrowCompute::arrow_compute_static)
  find_package(ArrowCompute REQUIRED QUIET)
endif()
# We are finding hictk here to avoid name collisions with the ENABLE_* variables defined in the coming lines
find_package(hictk REQUIRED QUIET)

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

# gitignore build folders
file(WRITE "${CMAKE_BINARY_DIR}/.gitignore" "*")

target_compile_features(
  hictkpy_project_options
  INTERFACE
    "cxx_std_$<IF:$<BOOL:${CMAKE_CXX_STANDARD}>,17,${CMAKE_CXX_STANDARD}>"
)

option(HICTKPY_BUILD_PORTABLE_WHEELS_WIN "Build portable Windows wheels" OFF)
option(HICTKPY_ENABLE_TESTING "Build native C++ tests" OFF)

# Tweak fmt
target_compile_definitions(
  hictkpy_project_options
  INTERFACE
    # Tweak fmt
    FMT_HEADER_ONLY
    # starting with fmt v12, clang-tidy breaks when FMT_ENFORCE_COMPILE_STRING is defined
    # https://github.com/fmtlib/fmt/commit/619b3a5aa031a25e8d06fec1e03992e5e727c5a5#diff-bdc6f79e8e9f5b4331d66fb785636a87d29f55cf729865e13925b4209424c878R4223
    $<$<NOT:$<BOOL:OPT_ENABLE_CLANG_TIDY>>:FMT_ENFORCE_COMPILE_STRING>
    FMT_USE_FULL_CACHE_DRAGONBOX
    HICTK_WITH_ARROW
    HICTK_WITH_EIGEN
    # Tweak spdlog
    SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG
    SPDLOG_CLOCK_COARSE
    SPDLOG_FMT_EXTERNAL
    SPDLOG_NO_THREAD_ID
    # Windows-specific tweaks
    $<$<BOOL:${WIN32}>:NOMINMAX>
    $<$<BOOL:${WIN32}>:_CRT_SECURE_NO_WARNINGS>
    $<$<BOOL:${WIN32}>:_DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR> # https://github.com/gabime/spdlog/issues/3212
)

if(WIN32 AND HICTKPY_BUILD_PORTABLE_WHEELS_WIN)
  target_compile_definitions(
    hictkpy_project_options
    INTERFACE
      _WIN32_WINNT=_WIN32_WINNT_WIN10
      /FI<sdkddkver.h> # required to use _WIN32_WINNT_WIN10
  )
endif()

add_subdirectory(src)
