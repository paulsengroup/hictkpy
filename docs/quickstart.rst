..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Quickstart
##########

hictkpy provides Python bindings for hictk through pybind11.

``hictk.File()`` can open .cool and .hic files and allows retrieval of interactions as well as file metadata.

The example use file `4DNFIOTPSS3L.hic <https://data.4dnucleome.org/files-processed/4DNFIOTPSS3L>`_, which can be downloaded from `here <https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/7386f953-8da9-47b0-acb2-931cba810544/4DNFIOTPSS3L.hic>`_.

Read file metadata
------------------

.. code-block:: ipythonconsole

  In [1]: import hictkpy as htk

  In [2]: f = htk.File("4DNFIOTPSS3L.hic", 10_000)

  In [3]: f.path()
  Out[3]: '4DNFIOTPSS3L.hic'

  In [4]: f.bin_size()
  Out[4]: 10000

  In [5]: f.chromosomes()
  Out[5]:
  {'2L': 23513712,
   '2R': 25286936,
   '3L': 28110227,
   '3R': 32079331,
   '4': 1348131,
   'X': 23542271,
   'Y': 3667352}

  In [6]: f.attributes()
  Out[6]:
  {'bin_size': 10000,
   'format': 'HIC',
   'format_version': 8,
   'assembly': '/var/lib/cwl/stgb25a903a-ebb6-4a56-bf3f-90bd84a40bf4/4DNFIBEEN92C.chrom.sizes',
   'format-url': 'https://github.com/aidenlab/hic-format',
   'nbins': 13758,
   'nchroms': 8}


Fetch interactions
------------------

.. code-block:: ipythonconsole

  In [1]: import hictkpy as htk

  In [2]: f = htk.File("4DNFIOTPSS3L.hic", 10_000)

  In [3]: sel = f.fetch("2L:10,000,000-20,000,000", join=True)

  In [4]: sel.nnz()
  Out[4]: 339041

  In [5]: sel.sum()
  Out[5]: 7163361

  In [6]: sel.to_df()
  Out[6]:
         chrom1    start1      end1 chrom2    start2      end2  count
  0          2L  10000000  10010000     2L  10000000  10010000   6759
  1          2L  10000000  10010000     2L  10010000  10020000   3241
  2          2L  10000000  10010000     2L  10020000  10030000    760
  3          2L  10000000  10010000     2L  10030000  10040000    454
  4          2L  10000000  10010000     2L  10040000  10050000    289
  ...       ...       ...       ...    ...       ...       ...    ...
  339036     2L  19970000  19980000     2L  19980000  19990000    407
  339037     2L  19970000  19980000     2L  19990000  20000000    221
  339038     2L  19980000  19990000     2L  19980000  19990000    391
  339039     2L  19980000  19990000     2L  19990000  20000000    252
  339040     2L  19990000  20000000     2L  19990000  20000000    266

  [339041 rows x 7 columns]

  In [7]: sel.to_coo()
  Out[7]:
  <1000x1000 sparse matrix of type '<class 'numpy.int32'>'
          with 339041 stored elements in COOrdinate format>

  In [8]: m = sel.to_numpy()

  In [9]: import matplotlib.pyplot as plt

  In [10]: from matplotlib.colors import LogNorm

  In [11]: plt.imshow(m, norm=LogNorm())

  In [12]: plt.show()

.. image:: assets/heatmap_001.avif

.. code-block:: ipythonconsole

  In [13]: plt.clf()

  In [37]: sel = f.fetch("2L:10,000,000-20,000,000", "X")

  In [38]: m = sel.to_numpy()

  In [39]: plt.imshow(m, norm=LogNorm())
  Out[39]: <matplotlib.image.AxesImage at 0x7faadbcb1150>

  In [40]: plt.savefig("/tmp/test.png", dpi=600)

.. image:: assets/heatmap_002.avif
