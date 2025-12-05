# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Catch2 REQUIRED QUIET)
include(CTest)
include(Catch)

add_executable(hictkpy_unit_tests)

target_sources(hictkpy_unit_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/pixel_table_test.cpp")

target_link_libraries(hictkpy_unit_tests PRIVATE hictkpy_pixel_table)

target_link_system_libraries(
  hictkpy_unit_tests
  PRIVATE
    Catch2::Catch2WithMain
)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Testing/")

catch_discover_tests(
  hictkpy_unit_tests
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/Testing/"
  EXTRA_ARGS
  --skip-benchmarks
)
