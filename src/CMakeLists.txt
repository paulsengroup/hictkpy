# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(
  Python 3.9
  COMPONENTS Interpreter Development.Module
  REQUIRED)

nanobind_add_module(
  _hictkpy
  NB_STATIC
  LTO
  MODULE
  hictkpy.cpp
  hictkpy_file.cpp
  hictkpy_file_creation.cpp
  hictkpy_multires_file.cpp
  hictkpy_pixel_selector.cpp
  hictkpy_singlecell_file.cpp)

target_include_directories(_hictkpy PRIVATE include)
target_link_libraries(
  _hictkpy
  PRIVATE hictkpy_project_options
          hictkpy_project_warnings
          hictk::cooler
          hictk::file
          hictk::hic)

install(TARGETS _hictkpy LIBRARY DESTINATION hictkpy)
