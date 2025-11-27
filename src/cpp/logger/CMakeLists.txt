# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictkpy_logger STATIC)

target_sources(hictkpy_logger PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/logger.cpp")

target_include_directories(hictkpy_logger PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_logger
  PRIVATE
    hictkpy_common
    hictkpy_locking
    hictkpy_project_options
    hictkpy_project_warnings
    fmt::fmt-header-only
  PUBLIC
    hictkpy_nanobind
    spdlog::spdlog_header_only
)
