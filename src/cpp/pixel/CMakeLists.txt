# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED QUIET)

add_library(hictkpy_pixel STATIC)

target_sources(
  hictkpy_pixel
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_impl.hpp"
)

target_include_directories(hictkpy_pixel PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_pixel
  PRIVATE
    hictkpy_common
    hictkpy_project_options
    hictkpy_project_warnings
    hictkpy_type
    fmt::fmt-header-only
  PUBLIC
    hictkpy_nanobind
    hictk::bin_table
    hictk::pixel
)
