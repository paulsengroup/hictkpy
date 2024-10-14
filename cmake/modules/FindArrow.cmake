# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

#[=======================================================================[.rst:
FindArrow
---------

Finds the Arrow library that is shipped as part of the pyarrow wheels.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``Arrow::arrow_shared``
  The Arrow shared library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``Arrow_FOUND``
  True if the system has the Arrow library.
``Arrow_VERSION``
  The version of the Arrow library which was found.
``Arrow_INCLUDE_DIRS``
  Include directories needed to use Arrow.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``Arrow_LIBRARY``
  The path to the Arrow library.

#]=======================================================================]

find_package(Pyarrow)
