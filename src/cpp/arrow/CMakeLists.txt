# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Arrow REQUIRED QUIET)
if(NOT TARGET ArrowCompute::arrow_compute_shared AND NOT TARGET ArrowCompute::arrow_compute_static)
  find_package(ArrowCompute REQUIRED QUIET)
endif()

set(IS_WINDOWS "$<BOOL:${WIN32}>")
set(STATIC_BUILD "$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>")

# gersemi: off
add_library(hictkpy_arrow INTERFACE)
target_link_libraries(
  hictkpy_arrow
  INTERFACE
    "Arrow::arrow_$<IF:${STATIC_BUILD},static,shared>"
)
# TODO ARROW_STATIC should've been defined by Conan
target_compile_definitions(
  hictkpy_arrow
  INTERFACE
    "$<$<AND:${IS_WINDOWS},${STATIC_BUILD}>:ARROW_STATIC>"
)

add_library(hictkpy_arrow_compute INTERFACE)
target_link_libraries(
  hictkpy_arrow_compute
  INTERFACE
    "ArrowCompute::arrow_compute_$<IF:${STATIC_BUILD},static,shared>"
)
# TODO ARROW_COMPUTE_STATIC should've been defined by Conan
target_compile_definitions(
  hictkpy_arrow_compute
  INTERFACE
    "$<$<AND:${IS_WINDOWS},${STATIC_BUILD}>:ARROW_COMPUTE_STATIC>"
)
# gersemi: on
