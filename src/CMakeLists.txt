# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
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

find_package(Arrow REQUIRED QUIET)
if(NOT TARGET ArrowCompute::arrow_compute_shared AND NOT TARGET ArrowCompute::arrow_compute_static)
  find_package(ArrowCompute REQUIRED QUIET)
endif()
find_package(FMT REQUIRED QUIET)
find_package(HDF5 REQUIRED QUIET COMPONENTS C)
find_package(nanobind REQUIRED QUIET)
find_package(phmap REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

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
        _hictkpy
        NB_DOMAIN hictkpy
        NB_STATIC
        ${HICTKPY_MODULE_ARGS}
        NB_SUPPRESS_WARNINGS
        MODULE
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/bin_table/bin_table.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/file/file.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/file/multires_file.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/file/singlecell_file.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/file_writer/cooler_file_writer.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/file_writer/file_writer_helpers.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/file_writer/hic_file_writer.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/locking/locking.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/logger/logger.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/nanobind/nanobind.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/nanobind/nanobind_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/numpy/to_numpy_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/pixel/pixel.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/pixel/pixel_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/pixel_aggregator/pixel_aggregator_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/pixel_selector/pixel_selector.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/reference/reference.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/table/pixel_table.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/table/pixel_table_bg2.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/table/pixel_table_coo.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/table/pixel_table_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/table/table.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/table/table_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/type/type.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/hictkpy.cpp"
)

target_include_directories(
  _hictkpy
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/bin_table/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/common/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/file/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/file_writer/include/cooler"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/file_writer/include/helpers"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/file_writer/include/hic"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/locking/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/logger/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/nanobind/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/numpy/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/pixel/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/pixel_aggregator/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/pixel_selector/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/reference/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/table/include/pixel_table"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/table/include/table"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/type/include"
)

target_link_libraries(
  _hictkpy
  PRIVATE
    hictkpy_project_options
    hictkpy_project_warnings
    Arrow::arrow_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>
    ArrowCompute::arrow_compute_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>
    fmt::fmt-header-only
    HDF5::HDF5
    hictk::balancing
    hictk::bin_table
    hictk::common
    hictk::cooler
    hictk::file
    hictk::format
    hictk::genomic_interval
    hictk::hic
    hictk::pixel
    hictk::reference
    hictk::tmpdir
    hictk::transformers
    hictk::type_traits
    Python::Module
    Python::NumPy
    spdlog::spdlog_header_only
    std::filesystem
)

install(TARGETS _hictkpy LIBRARY DESTINATION hictkpy)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/hictkpy")

if(CMAKE_BUILD_TYPE STREQUAL Release)
  # This properly escapes Windows paths
  string(REPLACE "\\" "/" HICTKPY_PROJECT_SOURCE_DIR ${PROJECT_SOURCE_DIR})
  string(REPLACE "\\" "/" HICTKPY_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
  string(REPLACE "\\" "/" HICTKPY_CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
  string(REPLACE "\\" "/" HICTKPY_Python_EXECUTABLE ${Python_EXECUTABLE})

  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/NanobindStubgen.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/NanobindStubgen.cmake"
    @ONLY
  )
  install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/NanobindStubgen.cmake")

  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/py.typed" DESTINATION hictkpy)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/__init__.pyi" DESTINATION hictkpy)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/cooler.pyi" DESTINATION hictkpy)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/hic.pyi" DESTINATION hictkpy)

  # Disable clang-tidy on nanobind
  # This seems to be the only reliable way to do so...
  configure_file("${PROJECT_SOURCE_DIR}/cmake/.clang-tidy.in" "${nanobind_SOURCE_DIR}/.clang-tidy" COPYONLY)
endif()

if(HICTKPY_ENABLE_TESTING)
  add_subdirectory(cpp)
endif()
