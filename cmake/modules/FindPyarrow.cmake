# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

#[=======================================================================[.rst:
FindPyarrow
-----------

Finds the pyarrow library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``arrow_python``
  The arrow_python library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``Foo_FOUND``
  True if the system has the Foo library.
``Foo_VERSION``
  The version of the Foo library which was found.
``Foo_INCLUDE_DIRS``
  Include directories needed to use Foo.
``Foo_LIBRARIES``
  Libraries needed to link to Foo.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``Foo_INCLUDE_DIR``
  The directory containing ``foo.h``.
``Foo_LIBRARY``
  The path to the Foo library.

#]=======================================================================]

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
