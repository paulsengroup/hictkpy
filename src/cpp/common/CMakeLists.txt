# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictkpy_common INTERFACE)

target_include_directories(hictkpy_common INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_common
  INTERFACE
    hictkpy_project_options
    hictkpy_project_warnings
    hictk::common
    hictk::type_traits
)
