# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

#[=======================================================================[.rst:
FindPyarrow
---------

Finds pyarrow library that is shipped as part of the pyarrow wheels.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``Arrow::python
  The pyarrow library


Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``Pyarrow_FOUND``
  True if the system has the Pyarrow library.
``Pyarrow_VERSION``
  The version of the Pyarrow library which was found.
``Pyarrow_INCLUDE_DIRS``
  Include directories needed to use Pyarrow.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``Pyarrow_LIBRARY``
  The path to the Arrow::python library.

#]=======================================================================]

find_package(
  Python
  3.9
  COMPONENTS
    Interpreter
    NumPy
  REQUIRED
)

if(TARGET Arrow::python)
  message(DEBUG "Arrow::python has already been defined")
  return()
endif()

# Try to import pyarrow
execute_process(
  COMMAND
    "${Python_EXECUTABLE}" -c "import pyarrow"
  RESULT_VARIABLE STATUS
)
if(NOT STATUS EQUAL 0)
  message(WARNING "Failed to import pyarrow")
  set(Pyarrow_FOUND FALSE)
  return()
endif()

# Get pyarrow version
execute_process(
  COMMAND
    "${Python_EXECUTABLE}" -c "from importlib.metadata import version; print(version('pyarrow'), end='')"
  RESULT_VARIABLE STATUS
  OUTPUT_VARIABLE Pyarrow_VERSION
)
if(NOT STATUS EQUAL 0)
  message(WARNING "Unable to detect pyarrow version")
  set(Pyarrow_FOUND FALSE)
  unset(Pyarrow_VERSION)
  return()
endif()

if(Pyarrow_VERSION VERSION_LESS 15.0.0)
  message(WARNING "pyarrow version ${Pyarrow_VERSION} is too old. Minimum version required: 15.0.0")
  set(Pyarrow_FOUND FALSE)
  unset(Pyarrow_VERSION)
  return()
endif()

set(Pyarrow_VERSION_STRING "${Pyarrow_VERSION}")

# Get include dirs
execute_process(
  COMMAND
    "${Python_EXECUTABLE}" -c "import pyarrow; print(pyarrow.get_include(), end='')"
  RESULT_VARIABLE STATUS
  OUTPUT_VARIABLE Pyarrow_INCLUDE_DIRS
)
if(NOT STATUS EQUAL 0 OR Pyarrow_INCLUDE_DIRS STREQUAL "")
  message(WARNING "Unable to detect pyarrow include directories")
  set(Pyarrow_FOUND FALSE)
  unset(Pyarrow_VERSION)
  unset(Pyarrow_INCLUDE_DIRS)
  return()
endif()

# Get library dirs
execute_process(
  COMMAND
    "${Python_EXECUTABLE}" -c "import pyarrow; print(' '.join(pyarrow.get_library_dirs()), end='')"
  RESULT_VARIABLE STATUS
  OUTPUT_VARIABLE Pyarrow_LIBRARY_DIRS
)
if(NOT STATUS EQUAL 0 OR Pyarrow_LIBRARY_DIRS STREQUAL "")
  message(WARNING "Unable to detect pyarrow library directories")
  set(Pyarrow_FOUND FALSE)
  unset(Pyarrow_VERSION)
  unset(Pyarrow_INCLUDE_DIRS)
  unset(Pyarrow_LIBRARY_DIRS)
  return()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Pyarrow
  FOUND_VAR Pyarrow_FOUND
  REQUIRED_VARS
    Pyarrow_INCLUDE_DIRS
    Pyarrow_LIBRARY_DIRS
  VERSION_VAR Pyarrow_VERSION
)

find_library(Pyarrow_LIBRARY arrow_python REQUIRED PATHS ${Pyarrow_LIBRARY_DIRS} NO_DEFAULT_PATH)

if(Pyarrow_FOUND AND NOT TARGET Arrow::python)
  if(WIN32)
    add_library(Arrow::python UNKNOWN IMPORTED)
  else()
    add_library(Arrow::python SHARED IMPORTED)
  endif()
  set_target_properties(
    Arrow::python
    PROPERTIES
      IMPORTED_LOCATION
        "${Pyarrow_LIBRARY}"
      IMPORTED_CONFIGURATION
        Release
  )
  target_include_directories(
    Arrow::python
    INTERFACE
      ${Pyarrow_INCLUDE_DIRS}
      ${Python_NumPy_INCLUDE_DIR}
  )
  target_link_directories(Arrow::python INTERFACE ${Pyarrow_LIBRARY_DIRS})
endif()

mark_as_advanced(
  Pyarrow_VERSION
  Pyarrow_INCLUDE_DIRS
  Pyarrow_LIBRARY_DIRS
)
