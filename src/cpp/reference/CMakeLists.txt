# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictkpy_reference STATIC)

target_sources(hictkpy_reference PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/reference.cpp")

target_include_directories(hictkpy_reference PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_reference
  PRIVATE
    hictkpy_locking
    hictkpy_project_options
    hictkpy_project_warnings
  PUBLIC
    hictkpy_nanobind
    hictk::reference
)
