..
   Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Cooler API
##########

.. py:module:: hictkpy.cooler
.. py:currentmodule:: hictkpy.cooler

.. autoclass:: SingleCellFile

   .. automethod:: __init__
   .. automethod:: __getitem__
   .. automethod:: __enter__
   .. automethod:: __exit__
   .. automethod:: attributes
   .. automethod:: bins
   .. automethod:: cells
   .. automethod:: chromosomes
   .. automethod:: close
   .. automethod:: path
   .. automethod:: resolution

.. autoclass:: FileWriter

   .. automethod:: __init__
   .. automethod:: __enter__
   .. automethod:: __exit__
   .. automethod:: add_pixels
   .. automethod:: add_pixels_from_dict
   .. automethod:: bins
   .. automethod:: chromosomes
   .. automethod:: finalize
   .. automethod:: path
   .. automethod:: resolution
