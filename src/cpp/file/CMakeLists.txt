# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED QUIET)
find_package(HDF5 REQUIRED QUIET COMPONENTS C)
find_package(phmap REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_library(hictkpy_file STATIC)

target_sources(
  hictkpy_file
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/file.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/multires_file.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/singlecell_file.cpp"
)

target_include_directories(hictkpy_file PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_file
  PRIVATE
    hictkpy_arrow
    hictkpy_bin_table
    hictkpy_common
    hictkpy_numpy
    hictkpy_pixel_selector
    hictkpy_project_options
    hictkpy_project_warnings
    hictkpy_reference
    hictkpy_table
    hictkpy_type
    HDF5::HDF5
    fmt::fmt-header-only
    hictk::balancing
    hictk::bin_table
    hictk::cooler
    hictk::file
    hictk::genomic_interval
    hictk::hic
    phmap
  PUBLIC
    hictkpy_locking
    hictkpy_nanobind
    hictk::cooler
    hictk::hic
    hictk::file
    std::filesystem
)
