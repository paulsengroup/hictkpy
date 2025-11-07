..
   Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Python API Reference
####################

This section provides a detailed reference for the ``hictkpy`` Python API.

The API is organized in three categories:

* :doc:`Generic API <generic>` -- Documents classes and functions such as :py:class:`hictkpy.File` and :py:class:`hictkpy.PixelSelector` that are used to open, inspect, and fetch interactions from ``.[m]cool`` and ``.hic`` files.
  The documentation also covers several general-purpose classes and data structures such as :py:class:`hictkpy.BinTable`, :py:class:`hictkpy.Bin`, and :py:class:`hictkpy.Pixel` which are utilized across the API.
* :doc:`Cooler API <cooler>` -- Documents the :py:class:`hictkpy.cooler.SingleCellFile` and :py:class:`hictkpy.cooler.FileWriter` classes, which are used to access ``.scool`` files and create ``.cool`` files, respectively.
* :doc:`.hic API <hic>` -- Documents the :py:class:`hictkpy.hic.FileWriter` class, which is used to create ``.hic`` files.
* :doc:`Logging API <logging>` -- Documents how to tweak hictkpy's :py:class:`logging.Logger`.

.. toctree::
   :hidden:
   :maxdepth: 1

   generic
   logging
   cooler
   hic
