# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Arrow REQUIRED QUIET)
if(NOT TARGET ArrowCompute::arrow_compute_shared AND NOT TARGET ArrowCompute::arrow_compute_static)
  find_package(ArrowCompute REQUIRED QUIET)
endif()

add_library(hictkpy_arrow INTERFACE)
target_link_libraries(hictkpy_arrow INTERFACE Arrow::arrow_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>)

add_library(hictkpy_arrow_compute INTERFACE)
target_link_libraries(
  hictkpy_arrow_compute
  INTERFACE
    ArrowCompute::arrow_compute_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>
)
