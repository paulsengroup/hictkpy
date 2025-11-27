# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Python 3.10 COMPONENTS Development.Module REQUIRED QUIET)

find_package(FMT REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictkpy_locking INTERFACE)

target_include_directories(hictkpy_locking INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_locking
  INTERFACE
    hictkpy_common
    hictkpy_project_options
    hictkpy_project_warnings
    fmt::fmt-header-only
    Python::Module
    spdlog::spdlog_header_only
)
