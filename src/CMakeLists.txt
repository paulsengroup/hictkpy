# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(
  Python 3.9
  COMPONENTS Interpreter Development.Module
  REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import pyarrow; print(pyarrow.get_include(), end='')"
  RESULT_VARIABLE STATUS
  OUTPUT_VARIABLE PYARROW_INCLUDE_DIR)
if(STATUS EQUAL 0)
  message(STATUS "Found pyarrow include directory: ${PYARROW_INCLUDE_DIR}")
else()
  message(FATAL_ERROR "Unable to find pyarrow include directory")
endif()

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import pyarrow; print(' '.join(pyarrow.get_library_dirs()), end='')"
  RESULT_VARIABLE STATUS
  OUTPUT_VARIABLE PYARROW_LIB_DIRS)
if(STATUS EQUAL 0)
  message(STATUS "Found pyarrow link directory: ${PYARROW_LIB_DIRS}")
else()
  message(FATAL_ERROR "Unable to find pyarrow link directory directory")
endif()

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c
          "import pyarrow; print(' '.join(lib for lib in pyarrow.get_libraries() if lib != 'arrow'), end='')"
  RESULT_VARIABLE STATUS
  OUTPUT_VARIABLE PYARROW_LIBS)
if(STATUS EQUAL 0)
  message(STATUS "Found pyarrow libraries: ${PYARROW_LIBS}")
else()
  message(FATAL_ERROR "Unable to find pyarrow libraries")
endif()

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import numpy; print(numpy.get_include(), end='')"
  RESULT_VARIABLE STATUS
  OUTPUT_VARIABLE NUMPY_INCLUDE_DIR)
if(STATUS EQUAL 0)
  message(STATUS "Found numpy include directory: ${NUMPY_INCLUDE_DIR}")
else()
  message(FATAL_ERROR "Unable to find numpy include directory")
endif()

execute_process(COMMAND "${Python_EXECUTABLE}" -c "import pyarrow; pyarrow.create_library_symlinks()"
                RESULT_VARIABLE STATUS)
if(NOT
   STATUS
   EQUAL
   0)
  message(FATAL_ERROR "Unable to create symlink to arrow libraries")
endif()

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
  PRIVATE hictkpy_project_options hictkpy_project_warnings
  PUBLIC hictk::cooler
         hictk::file
         hictk::hic
         ${PYARROW_LIBS})

install(TARGETS _hictkpy LIBRARY DESTINATION hictkpy)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/hictkpy")
file(TOUCH "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/__init__.pyi")

include("${PROJECT_SOURCE_DIR}/cmake/NBStubgen.cmake")

_nanobind_add_stub(
  hictkpy._hictkpy_stub
  INSTALL_TIME
  MODULE
  hictkpy._hictkpy
  MARKER_FILE
  "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/py.typed"
  OUTPUT
  "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/hictkpy.pyi"
  LD_LIBRARY_PATH
  "${PYARROW_LIB_DIRS}"
  VERBOSE)

_nanobind_add_stub(
  hictkpy._hictkpy.cooler_stub
  INSTALL_TIME
  MODULE
  hictkpy._hictkpy.cooler
  OUTPUT
  "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/cooler.pyi"
  LD_LIBRARY_PATH
  "${PYARROW_LIB_DIRS}"
  VERBOSE)

_nanobind_add_stub(
  hictkpy._hictkpy.hic_stub
  INSTALL_TIME
  MODULE
  hictkpy._hictkpy.hic
  OUTPUT
  "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/hic.pyi"
  LD_LIBRARY_PATH
  "${PYARROW_LIB_DIRS}"
  VERBOSE)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/py.typed" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/__init__.pyi" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/hictkpy.pyi" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/cooler.pyi" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/hic.pyi" DESTINATION hictkpy)
