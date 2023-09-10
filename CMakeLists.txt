# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...3.27)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

# Not ideal to use this global variable, but necessary to make sure that tooling and projects use the same version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)

# strongly encouraged to enable this globally to avoid conflicts between -Wpedantic being enabled and -std=c++20 and
# -std=gnu++20 for example when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

project(
  hictkpy
  LANGUAGES C CXX
  HOMEPAGE_URL https://github.com/paulsengroup/hictkpy
  DESCRIPTION "Python bindings for hictk.")

include(FetchContent)
FetchContent_Declare(
  hictk
  URL       "${CMAKE_CURRENT_SOURCE_DIR}/external/hictk-v0.0.2.tar.xz"
  URL_HASH  "SHA256=e2df4eb708e51017f595ac7e54cb5fa5c72aadb36766fbaabb50af7cf67d3656"
  SYSTEM)

set(BUILD_SHARED_LIBS ON)
set(HICTK_ENABLE_TESTING OFF)
set(HICTK_BUILD_EXAMPLES OFF)
set(HICTK_BUILD_BENCHMARKS OFF)
set(HICTK_BUILD_PYTHON_BINDINGS OFF)
set(HICTK_WITH_EIGEN OFF)
set(HICTK_BUILD_TOOLS OFF)
set(HICTK_INSTALL OFF)
FetchContent_MakeAvailable(hictk)

add_library(hictkpy_project_options INTERFACE)
target_compile_features(hictkpy_project_options INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerWarnings.cmake")
if(WIN32)
  target_compile_definitions(hictkpy_project_options INTERFACE NOMINMAX _CRT_SECURE_NO_WARNINGS)
endif()

add_library(hictkpy_project_warnings INTERFACE)
target_compile_options(
        hictkpy_project_warnings
        INTERFACE # C++ warnings
        $<$<COMPILE_LANGUAGE:CXX>:${HICTKPY_PROJECT_WARNINGS_CXX}>
        # C warnings
        $<$<COMPILE_LANGUAGE:C>:${HICTKPY_PROJECT_WARNINGS_C}>
        # Cuda warnings
        $<$<COMPILE_LANGUAGE:CUDA>:${HICTKPY_PROJECT_WARNINGS_CUDA}>)

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif()

add_subdirectory(src)
