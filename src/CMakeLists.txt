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

# Automatic stub generation for hictkpy does not work, instead run:
# venv/bin/python -m nanobind.stubgen -m hictkpy -o src/hictkpy/hictkpy.pyi
# Remember to remove the version annotation from the stub

nanobind_add_stub(
  hictkpy._hictkpy_stub
  INSTALL_TIME
  MODULE
  hictkpy._hictkpy
  MARKER_FILE
  "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/py.typed"
  OUTPUT
  "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/_hictkpy.pyi"
  VERBOSE)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/hictkpy/hictkpy.pyi" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/_hictkpy.pyi" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/py.typed" DESTINATION hictkpy)
