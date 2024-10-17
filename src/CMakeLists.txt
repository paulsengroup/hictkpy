# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED)
find_package(nanobind REQUIRED)
find_package(Pyarrow REQUIRED)
find_package(spdlog REQUIRED)
find_package(Filesystem REQUIRED)

nanobind_add_module(
        _hictkpy
        NB_STATIC
        NB_DOMAIN hictkpy
        LTO
        MODULE
        "${CMAKE_CURRENT_SOURCE_DIR}/common.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cooler_file_writer.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/file.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/hic_file_writer.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/hictkpy.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/logger.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/multires_file.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/pixel_selector.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/reference.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/singlecell_file.cpp"
)

target_include_directories(
  _hictkpy
  PRIVATE
    include
    "${PYARROW_INCLUDE_DIR}"
    "${NUMPY_INCLUDE_DIR}"
)
target_link_directories(_hictkpy PRIVATE "${PYARROW_LIB_DIRS}")

target_link_libraries(
  _hictkpy
  PRIVATE
    hictkpy_project_options
    hictkpy_project_warnings
    hictk::cooler
    hictk::file
    hictk::hic
    hictk::transformers
    Arrow::python
    fmt::fmt-header-only
    spdlog::spdlog_header_only
    std::filesystem
)

# TODO fixme
# There is something very odd going on when trying to call most spdlog functions from within
# Python bindings on recent versions of Windows.
# Possibly related to https://github.com/gabime/spdlog/issues/3212
target_compile_definitions(_hictkpy PRIVATE $<$<PLATFORM_ID:Windows>:SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_OFF>)

install(TARGETS _hictkpy LIBRARY DESTINATION hictkpy)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/hictkpy")
file(TOUCH "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/__init__.pyi")

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
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/hictkpy.pyi" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/cooler.pyi" DESTINATION hictkpy)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkpy/hic.pyi" DESTINATION hictkpy)

# Disable clang-tidy on nanobind
# This seems to be the only reliable way to do so...
configure_file("${PROJECT_SOURCE_DIR}/cmake/.clang-tidy.in" "${nanobind_SOURCE_DIR}/.clang-tidy" COPYONLY)
