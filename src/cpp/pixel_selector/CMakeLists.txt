# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED QUIET)
find_package(HDF5 REQUIRED QUIET COMPONENTS C)
find_package(phmap REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_library(hictkpy_pixel_selector STATIC)

target_sources(hictkpy_pixel_selector PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/pixel_selector.cpp")

target_include_directories(hictkpy_pixel_selector PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_pixel_selector
  PRIVATE
    hictkpy_pixel
    hictkpy_pixel_aggregator
    hictkpy_project_options
    hictkpy_project_warnings
    hictkpy_table
    hictkpy_type
    fmt::fmt-header-only
    hictk::format
    phmap
  PUBLIC
    hictkpy_common
    hictkpy_locking
    hictkpy_nanobind
    hictk::bin_table
    hictk::cooler
    hictk::hic
    hictk::pixel
    hictk::transformers
)
