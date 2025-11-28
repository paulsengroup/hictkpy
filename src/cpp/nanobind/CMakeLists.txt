# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(
  Python
  3.10
  COMPONENTS
    Interpreter
    Development.Module
    NumPy
  REQUIRED
)

find_package(FMT REQUIRED QUIET)
find_package(nanobind REQUIRED QUIET)

if(NOT TARGET nanobind-static)
  set(HICTKPY_NANOBIND_DUMMY_FILE_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/nanobind/")
  set(HICTKPY_NANOBIND_DUMMY_FILE "${HICTKPY_NANOBIND_DUMMY_FILE_DIR}/dummy.cpp")
  file(MAKE_DIRECTORY "${HICTKPY_NANOBIND_DUMMY_FILE_DIR}")
  file(TOUCH "${HICTKPY_NANOBIND_DUMMY_FILE}")

  if(CMAKE_BUILD_TYPE STREQUAL Release)
    set(
      HICTKPY_MODULE_ARGS
      LTO
      NOMINSIZE
    )
  else()
    set(HICTKPY_MODULE_ARGS NOSTRIP)
  endif()

  nanobind_add_module(
          hictkpy_nanobind_dummy
          NB_DOMAIN hictkpy
          NB_STATIC
          ${HICTKPY_MODULE_ARGS}
          NB_SUPPRESS_WARNINGS
          MODULE
          "${HICTKPY_NANOBIND_DUMMY_FILE}"
  )
endif()

add_library(hictkpy_nanobind STATIC)

target_sources(
  hictkpy_nanobind
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/nanobind.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/nanobind_impl.hpp"
)

target_include_directories(hictkpy_nanobind PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_nanobind
  PRIVATE
    hictk::numeric
    hictkpy_project_options
    hictkpy_project_warnings
  PUBLIC
    hictkpy_locking
    fmt::fmt-header-only
    nanobind-static
    Python::Module
    Python::NumPy
)
