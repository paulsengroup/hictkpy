..
   Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Generic API
===========

.. py:module:: hictkpy
.. py:currentmodule:: hictkpy


.. autofunction:: is_cooler
.. autofunction:: is_mcool_file
.. autofunction:: is_scool_file

.. autofunction:: is_hic

.. autoclass:: MultiResFile

   .. automethod:: __init__
   .. automethod:: chromosomes
   .. automethod:: path
   .. automethod:: resolutions

.. autoclass:: File

   .. automethod:: __init__
   .. automethod:: attributes
   .. automethod:: avail_normalizations
   .. automethod:: bins
   .. automethod:: chromosomes
   .. automethod:: fetch
   .. automethod:: has_normalization
   .. automethod:: is_cooler
   .. automethod:: is_hic
   .. automethod:: nbins
   .. automethod:: nchroms
   .. automethod:: path
   .. automethod:: resolution
   .. automethod:: uri
   .. automethod:: weights

.. autoclass:: PixelSelector

   .. automethod:: coord1
   .. automethod:: coord2
   .. automethod:: nnz
   .. automethod:: sum
   .. automethod:: to_coo
   .. automethod:: to_df
   .. automethod:: to_numpy

.. autoclass:: Bin

.. autoclass:: BinTable

   .. automethod:: __init__
   .. automethod:: chromosomes
   .. automethod:: get
   .. automethod:: get_id
   .. automethod:: get_ids
   .. automethod:: merge
   .. automethod:: resolution
   .. automethod:: to_df
   .. automethod:: type
