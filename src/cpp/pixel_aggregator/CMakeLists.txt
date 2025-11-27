# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED QUIET)
find_package(phmap REQUIRED QUIET)

add_library(hictkpy_pixel_aggregator INTERFACE)

target_include_directories(hictkpy_pixel_aggregator INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_pixel_aggregator
  INTERFACE
    hictkpy_common
    hictkpy_project_options
    hictkpy_project_warnings
    phmap
)
