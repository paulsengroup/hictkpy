# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Python 3.10 COMPONENTS Development.Module REQUIRED QUIET)

find_package(FMT REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictkpy_locking STATIC)

target_sources(hictkpy_locking PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/locking.cpp")

target_include_directories(hictkpy_locking PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_locking
  PRIVATE
    hictkpy_common
    hictkpy_project_options
    hictkpy_project_warnings
    fmt::fmt-header-only
    spdlog::spdlog_header_only
  PUBLIC
    Python::Module
)
