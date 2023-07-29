# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(
  Python 3.6
  COMPONENTS Interpreter Development.Module
  REQUIRED)

find_package(Filesystem REQUIRED)

find_package(pybind11 CONFIG REQUIRED)

# TODO: rename module after moving hictkpy out of hictk
pybind11_add_module(hictkpy_ MODULE hictkpy.cpp)

target_link_libraries(
  hictkpy_
  PRIVATE hictkpy_project_options hictkpy_project_warnings
  PUBLIC hictk::cooler hictk::file hictk::hic std::filesystem)
