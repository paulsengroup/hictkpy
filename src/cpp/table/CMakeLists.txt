# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictkpy_table STATIC)

target_sources(
  hictkpy_table
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/table.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/table_impl.hpp"
)

target_include_directories(hictkpy_table PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/table")
target_link_libraries(
  hictkpy_table
  PRIVATE
    hictkpy_locking
    hictkpy_project_options
    hictkpy_project_warnings
    hictkpy_type
    fmt::fmt-header-only
  PUBLIC
    hictkpy_arrow
    hictkpy_common
    hictkpy_nanobind
)

add_library(hictkpy_pixel_table STATIC)

target_sources(
  hictkpy_pixel_table
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_table.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_table_coo.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_table_bg2.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_table_impl.hpp"
)

target_include_directories(hictkpy_pixel_table PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/pixel_table")
target_link_libraries(
  hictkpy_pixel_table
  PRIVATE
    hictkpy_project_options
    hictkpy_project_warnings
    hictk::reference
  PUBLIC
    hictkpy_arrow
    hictkpy_arrow_compute
    hictkpy_common
    hictkpy_table
    hictkpy_type
    fmt::fmt-header-only
    hictk::bin_table
    hictk::pixel
    spdlog::spdlog_header_only
)
