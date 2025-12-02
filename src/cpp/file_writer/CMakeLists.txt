# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_library(hictkpy_file_writer_helper STATIC)

target_sources(hictkpy_file_writer_helper PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/file_writer_helpers.cpp")

target_include_directories(hictkpy_file_writer_helper PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/helpers")
target_link_libraries(
  hictkpy_file_writer_helper
  PRIVATE
    hictkpy_common
    hictkpy_project_options
    hictkpy_project_warnings
    fmt::fmt-header-only
  PUBLIC
    hictkpy_arrow
    hictkpy_nanobind
    hictkpy_table
)

add_library(hictkpy_file_writer_cooler STATIC)

target_sources(hictkpy_file_writer_cooler PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/cooler_file_writer.cpp")

target_include_directories(hictkpy_file_writer_cooler PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/cooler")
target_link_libraries(
  hictkpy_file_writer_cooler
  PRIVATE
    hictkpy_arrow
    hictkpy_common
    hictkpy_file_writer_helper
    hictkpy_locking
    hictkpy_pixel_table
    hictkpy_project_options
    hictkpy_project_warnings
    hictkpy_type
    fmt::fmt-header-only
    hictk::file
    spdlog::spdlog_header_only
  PUBLIC
    hictkpy_bin_table
    hictkpy_file
    hictkpy_nanobind
    hictkpy_reference
    hictk::bin_table
    hictk::cooler
    hictk::reference
    hictk::tmpdir
    std::filesystem
)

add_library(hictkpy_file_writer_hic STATIC)

target_sources(hictkpy_file_writer_hic PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/hic_file_writer.cpp")

target_include_directories(hictkpy_file_writer_hic PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/hic")
target_link_libraries(
  hictkpy_file_writer_hic
  PRIVATE
    hictkpy_common
    hictkpy_file_writer_helper
    hictkpy_locking
    hictkpy_numpy
    hictkpy_pixel_table
    hictkpy_project_options
    hictkpy_project_warnings
    fmt::fmt-header-only
    hictk::bin_table
    hictk::file
    spdlog::spdlog_header_only
  PUBLIC
    hictkpy_bin_table
    hictkpy_file
    hictkpy_nanobind
    hictkpy_reference
    hictk::hic
    hictk::reference
    hictk::tmpdir
)
