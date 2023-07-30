# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(
  Python 3.6
  COMPONENTS Interpreter Development.Module
  REQUIRED)

# For some reason linking to std::filesystem breaks cibw builds for Apple Silicon
# find_package(Filesystem REQUIRED)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(hictkpy MODULE hictkpy.cpp)

target_link_libraries(
  hictkpy
  PRIVATE hictkpy_project_options hictkpy_project_warnings
  PUBLIC hictk::cooler
         hictk::file
         hictk::hic)
