# Copyright (c) 2022 Wenzel Jakob <wenzel.jakob@epfl.ch>, All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
# following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
#   disclaimer.
#
# 1. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
#   disclaimer in the documentation and/or other materials provided with the distribution.
#
# 1. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
#   products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# This file contains a modified version of the nanobind_add_stub CMake function Source:
# https://github.com/wjakob/nanobind/blob/master/cmake/nanobind-config.cmake The function was modified to accept a
# LD_LIBRARY_PATH parameter to allow adding custom paths to the list of paths checked by the lib loader.

function(_nanobind_add_stub name)
  cmake_parse_arguments(
    PARSE_ARGV
    1
    ARG
    "VERBOSE;INCLUDE_PRIVATE;EXCLUDE_DOCSTRINGS;INSTALL_TIME;EXCLUDE_FROM_ALL"
    "MODULE;OUTPUT;MARKER_FILE;COMPONENT;PATTERN_FILE;LD_LIBRARY_PATH"
    "PYTHON_PATH;DEPENDS")

  if(EXISTS ${NB_DIR}/src/stubgen.py)
    set(NB_STUBGEN "${NB_DIR}/src/stubgen.py")
  elseif(EXISTS ${NB_DIR}/stubgen.py)
    set(NB_STUBGEN "${NB_DIR}/stubgen.py")
  else()
    message(FATAL_ERROR "nanobind_add_stub(): could not locate 'stubgen.py'!")
  endif()

  if(NOT ARG_VERBOSE)
    list(APPEND NB_STUBGEN_ARGS -q)
  else()
    set(NB_STUBGEN_EXTRA USES_TERMINAL)
  endif()

  if(ARG_INCLUDE_PRIVATE)
    list(APPEND NB_STUBGEN_ARGS -P)
  endif()

  if(ARG_EXCLUDE_DOCSTRINGS)
    list(APPEND NB_STUBGEN_ARGS -D)
  endif()

  foreach(TMP IN LISTS ARG_PYTHON_PATH)
    list(
      APPEND
      NB_STUBGEN_ARGS
      -i
      "${TMP}")
  endforeach()

  if(ARG_PATTERN_FILE)
    list(
      APPEND
      NB_STUBGEN_ARGS
      -p
      "${ARG_PATTERN_FILE}")
  endif()

  if(ARG_MARKER_FILE)
    list(
      APPEND
      NB_STUBGEN_ARGS
      -M
      "${ARG_MARKER_FILE}")
    list(APPEND NB_STUBGEN_OUTPUTS "${ARG_MARKER_FILE}")
  endif()

  if(NOT ARG_MODULE)
    message(FATAL_ERROR "nanobind_add_stub(): a 'MODULE' argument must be specified!")
  else()
    list(
      APPEND
      NB_STUBGEN_ARGS
      -m
      "${ARG_MODULE}")
  endif()

  if(NOT ARG_OUTPUT)
    message(FATAL_ERROR "nanobind_add_stub(): an 'OUTPUT' argument must be specified!")
  else()
    list(
      APPEND
      NB_STUBGEN_ARGS
      -o
      "${ARG_OUTPUT}")
    list(APPEND NB_STUBGEN_OUTPUTS "${ARG_OUTPUT}")
  endif()

  file(TO_CMAKE_PATH ${Python_EXECUTABLE} NB_Python_EXECUTABLE)

  set(NB_STUBGEN_CMD "${NB_Python_EXECUTABLE}" "${NB_STUBGEN}" ${NB_STUBGEN_ARGS})

  if(NOT ARG_INSTALL_TIME)
    add_custom_command(
      OUTPUT ${NB_STUBGEN_OUTPUTS}
      COMMAND ${NB_STUBGEN_CMD}
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      DEPENDS ${ARG_DEPENDS}
              "${NB_STUBGEN}"
              "${ARG_PATTERN_FILE}"
              ${NB_STUBGEN_EXTRA})
    add_custom_target(${name} ALL DEPENDS ${NB_STUBGEN_OUTPUTS})
  else()
    set(NB_STUBGEN_EXTRA "")
    if(ARG_COMPONENT)
      list(
        APPEND
        NB_STUBGEN_EXTRA
        COMPONENT
        ${ARG_COMPONENT})
    endif()
    if(ARG_EXCLUDE_FROM_ALL)
      list(APPEND NB_STUBGEN_EXTRA EXCLUDE_FROM_ALL)
    endif()
    # \${CMAKE_INSTALL_PREFIX} has same effect as $<INSTALL_PREFIX> This is for compatibility with CMake < 3.27. For
    # more info: https://github.com/wjakob/nanobind/issues/420#issuecomment-1971353531
    if(WIN32)
      install(
        CODE "set(CMD \"${NB_STUBGEN_CMD}\")\nset(ENV{PATH} \"$ENV{PATH};${ARG_LD_LIBRARY_PATH}\")\nexecute_process(\n COMMAND \$\{CMD\}\n WORKING_DIRECTORY \"\${CMAKE_INSTALL_PREFIX}\"\n)"
        ${NB_STUBGEN_EXTRA})
    else()
      install(
        CODE "set(CMD \"${NB_STUBGEN_CMD}\")\nset(ENV{LD_LIBRARY_PATH} \"$ENV{LD_LIBRARY_PATH}:${ARG_LD_LIBRARY_PATH}\")\nexecute_process(\n COMMAND \$\{CMD\}\n WORKING_DIRECTORY \"\${CMAKE_INSTALL_PREFIX}\"\n)"
        ${NB_STUBGEN_EXTRA})
    endif()

  endif()
endfunction()
