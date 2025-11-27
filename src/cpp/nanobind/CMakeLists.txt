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
find_package(spdlog REQUIRED QUIET)

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

add_library(hictkpy_nanobind INTERFACE)

target_include_directories(hictkpy_nanobind INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictkpy_nanobind
  INTERFACE
    hictkpy_locking
    hictkpy_project_options
    hictkpy_project_warnings
    fmt::fmt-header-only
    nanobind-static
    Python::Module
    Python::NumPy
    spdlog::spdlog_header_only
)
