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
   .. automethod:: __enter__
   .. automethod:: __exit__
   .. automethod:: attributes
   .. automethod:: chromosomes
   .. automethod:: close
   .. automethod:: is_hic
   .. automethod:: is_mcool
   .. automethod:: path
   .. automethod:: resolutions

.. autoclass:: File

   .. automethod:: __init__
   .. automethod:: __enter__
   .. automethod:: __exit__
   .. automethod:: attributes
   .. automethod:: avail_normalizations
   .. automethod:: bins
   .. automethod:: chromosomes
   .. automethod:: close
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
   .. automethod:: dtype
   .. automethod:: to_arrow
   .. automethod:: to_coo
   .. automethod:: to_csr
   .. automethod:: to_df
   .. automethod:: to_numpy
   .. automethod:: to_pandas
   .. automethod:: size

   **Statistics**

   :py:class:`hictkpy.PixelSelector` exposes several methods to compute or estimate several statistics efficiently.

   The main features of these methods are:

   * All statistics are computed by traversing the data only once and without caching interactions.
   * All methods can be tweaked to include or exclude non-finite values.
   * All functions implemented using short-circuiting to detect scenarios where the required statistics can be computed without traversing all pixels.

   The following statistics are guaranteed to be exact:

   * nnz
   * sum
   * min
   * max
   * mean

   The rest of the supported statistics (currently variance, skewness, and kurtosis) are estimated and are thus not guaranteed to be exact.
   However, in practice, the estimation is usually very accurate (relative error < 1.0e-6).

   You can instruct hictkpy to compute the exact statistics by passing ``exact=True`` to :py:meth:`hictkpy.PixelSelector.describe()` and related methods.
   It should be noted that for large queries this will result in slower computations and higher memory usage.

   .. automethod:: describe
   .. automethod:: kurtosis
   .. automethod:: max
   .. automethod:: mean
   .. automethod:: min
   .. automethod:: nnz
   .. automethod:: skewness
   .. automethod:: sum
   .. automethod:: variance

   **Iteration**

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

  .. autoproperty:: id
  .. autoproperty:: rel_id
  .. autoproperty:: chrom
  .. autoproperty:: start
  .. autoproperty:: end

.. autoclass:: BinTable

   .. automethod:: __init__
   .. automethod:: chromosomes
   .. automethod:: get
   .. automethod:: get_id
   .. automethod:: get_ids
   .. automethod:: merge
   .. automethod:: resolution
   .. automethod:: to_arrow
   .. automethod:: to_df
   .. automethod:: to_pandas
   .. automethod:: type

   .. automethod:: __iter__

.. autoclass:: Pixel

  .. autoproperty:: bin1_id
  .. autoproperty:: bin2_id
  .. autoproperty:: count

  The following properties are only available when pixels are in BG2 format.

  .. autoproperty:: bin1
  .. autoproperty:: bin2
  .. autoproperty:: chrom1
  .. autoproperty:: start1
  .. autoproperty:: end1
  .. autoproperty:: chrom2
  .. autoproperty:: start2
  .. autoproperty:: end2
