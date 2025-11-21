..
   Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Logging API
===========

.. py:module:: hictkpy.logging
.. py:currentmodule:: hictkpy.logging

.. autofunction:: setLevel

Using this method instead should be preferred to tweak the log level with :code:`logging.getLogger("hictkpy").setLevel("INFO")`,
as :py:meth:`hictkpy.logging.setLevel` sets the log level also on the underlying the C++ logger, thus avoiding needlessly
producing log messages on the C++ side to then discard them once they reach the Python logger.
