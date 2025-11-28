# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictkpy_numpy INTERFACE)

target_sources(hictkpy_numpy INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/to_numpy_impl.hpp")

target_include_directories(hictkpy_numpy INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_numpy
  INTERFACE
    hictkpy_locking
    hictkpy_nanobind
    hictkpy_project_options
    hictkpy_project_warnings
)
