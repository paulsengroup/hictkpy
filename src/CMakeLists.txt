# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Pyarrow REQUIRED)

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

target_include_directories(_hictkpy PRIVATE include "${PYARROW_INCLUDE_DIR}" "${NUMPY_INCLUDE_DIR}")
target_link_directories(_hictkpy PRIVATE "${PYARROW_LIB_DIRS}")

target_link_libraries(
  _hictkpy
  PRIVATE hictkpy_project_options
          hictkpy_project_warnings
          hictk::cooler
          hictk::file
          hictk::hic
          hictk::transformers
          Arrow::python
)

install(TARGETS _hictkpy LIBRARY DESTINATION hictkpy)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/hictkpy")
file(TOUCH "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/__init__.pyi")

configure_file("${PROJECT_SOURCE_DIR}/cmake/NanobindStubgen.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/NanobindStubgen.cmake" @ONLY)
install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/NanobindStubgen.cmake")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/py.typed" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/__init__.pyi" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/hictkpy.pyi" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/cooler.pyi" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/hic.pyi" DESTINATION hictkpy)
