# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(
  Python 3.8
  COMPONENTS Interpreter Development.Module
  REQUIRED)

set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(
  _hictkpy
  MODULE
  hictkpy.cpp
  hictkpy_file.cpp
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
