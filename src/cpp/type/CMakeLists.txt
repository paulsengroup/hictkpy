# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Python 3.10 COMPONENTS Development.Module REQUIRED QUIET)

find_package(FMT REQUIRED QUIET)

add_library(hictkpy_type STATIC)

target_sources(hictkpy_type PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/type.cpp")

target_include_directories(hictkpy_type PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_type
  PRIVATE
    hictkpy_locking
    hictkpy_project_options
    hictkpy_project_warnings
    fmt::fmt-header-only
  PUBLIC
    hictkpy_common
    hictkpy_nanobind
    Python::Module
)
