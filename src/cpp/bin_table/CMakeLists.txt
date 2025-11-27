# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(spdlog REQUIRED QUIET)

add_library(hictkpy_bin_table STATIC)

target_sources(hictkpy_bin_table PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/bin_table.cpp")

target_include_directories(hictkpy_bin_table PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_bin_table
  PRIVATE
    hictkpy_arrow
    hictkpy_locking
    hictkpy_project_options
    hictkpy_project_warnings
    hictkpy_table
    hictk::genomic_interval
    fmt::fmt-header-only
  PUBLIC
    hictkpy_nanobind
    hictkpy_reference
    hictk::bin_table
    hictk::reference
)
