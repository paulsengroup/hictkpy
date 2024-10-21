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
   .. automethod:: __getitem__
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
   .. automethod:: to_arrow
   .. automethod:: to_coo
   .. automethod:: to_csr
   .. automethod:: to_df
   .. automethod:: to_numpy
   .. automethod:: to_pandas

   .. automethod:: __iter__

    .. code-block:: ipythonconsole

      In [1]: import hictkpy as htk

      In [2]: f = htk.File("file.cool")

      In [3]: sel = f.fetch("chr2L:10,000,000-20,000,000")

      In [4]: for i, pixel in enumerate(sel):
         ...:     print(pixel.bin1_id, pixel.bin2_id, pixel.count)
         ...:     if i > 10:
         ...:         break
         ...:
      1000 1000 6759
      1000 1001 3241
      1000 1002 760
      1000 1003 454
      1000 1004 289
      1000 1005 674
      1000 1006 354
      1000 1007 124
      1000 1008 130
      1000 1009 105
      1000 1010 99
      1000 1011 120

    It is also possible to iterate over pixels together with their genomic coordinates by specifying ``join=True`` when calling :py:meth:`hictkpy.File.fetch()`:

      .. code-block:: ipythonconsole

        In [5]: sel = f.fetch("chr2L:10,000,000-20,000,000", join=True)

        In [6]: for i, pixel in enumerate(sel):
           ...:     print(
           ...:         pixel.chrom1, pixel.start1, pixel.end1,
           ...:         pixel.chrom2, pixel.start2, pixel.end2,
           ...:         pixel.count
           ...:     )
           ...:     if i > 10:
           ...:         break
           ...:
        chr2L 10000000 10010000 chr2L 10000000 10010000 6759
        chr2L 10000000 10010000 chr2L 10010000 10020000 3241
        chr2L 10000000 10010000 chr2L 10020000 10030000 760
        chr2L 10000000 10010000 chr2L 10030000 10040000 454
        chr2L 10000000 10010000 chr2L 10040000 10050000 289
        chr2L 10000000 10010000 chr2L 10050000 10060000 674
        chr2L 10000000 10010000 chr2L 10060000 10070000 354
        chr2L 10000000 10010000 chr2L 10070000 10080000 124
        chr2L 10000000 10010000 chr2L 10080000 10090000 130
        chr2L 10000000 10010000 chr2L 10090000 10100000 105
        chr2L 10000000 10010000 chr2L 10100000 10110000 99
        chr2L 10000000 10010000 chr2L 10110000 10120000 120

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

   .. automethod:: __iter__
