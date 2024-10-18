..
   Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Creating .cool and .hic files
#############################

hictkpy supports creating .cool and .hic files from pre-binned interactions in COO or BedGraph2 format.

The example in this section use file `4DNFIOTPSS3L.hic <https://data.4dnucleome.org/files-processed/4DNFIOTPSS3L>`_, which can be downloaded from `here <https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/7386f953-8da9-47b0-acb2-931cba810544/4DNFIOTPSS3L.hic>`_.

Preparation
-----------

The first step consists of converting interactions from ``4DNFIOTPSS3L.hic`` to bedGraph2 format.
This can be achieved using ``hictk dump`` (or alternatively with :py:meth:hictkpy.File.fetch()`.

.. code-block:: console

  user@dev:/tmp$ hictk dump --join 4DNFIOTPSS3L.hic --resolution 50000 > pixels.bg2

  user@dev:/tmp$ head pixels.bg2

    2L	0	50000	2L	0	50000	30211
    2L	0	50000	2L	50000	100000	13454
    2L	0	50000	2L	100000	150000	2560
    2L	0	50000	2L	150000	200000	911
    2L	0	50000	2L	200000	250000	753
    2L	0	50000	2L	250000	300000	846
    2L	0	50000	2L	300000	350000	530
    2L	0	50000	2L	350000	400000	378
    2L	0	50000	2L	400000	450000	630
    2L	0	50000	2L	450000	500000	756


Next, we also generate the list of chromosomes to use as reference.

.. code-block:: console

  user@dev:/tmp$ hictk dump -t chroms 4DNFIOTPSS3L.hic > chrom.sizes

  user@dev:/tmp$ head chrom.sizes.bg2

    2L	23513712
    2R	25286936
    3L	28110227
    3R	32079331
    4	1348131
    X	23542271
    Y	3667352


Ingesting interactions in a .cool file
--------------------------------------

.. code-block:: ipythonconsole

  In [1]: import hictkpy as htk

  In [2]: import pandas as pd

  # Create a dictionary mapping chromosome names to chromosome sizes
  In [3]: chroms = pd.read_table("chrom.sizes", names=["name", "length"])
      ...            .set_index("name")["length"]
      ...            .to_dict()

  In [4]: chroms
  Out[4]:
  {'2L': 23513712,
   '2R': 25286936,
   '3L': 28110227,
   '3R': 32079331,
   '4': 1348131,
   'X': 23542271,
   'Y': 3667352}

  # Initialize an empty .cool file
  In [5]: f = htk.cooler.FileWriter("out.cool", chroms, resolution=50_000)

  In [6]: cols = ["chrom1", "start1", "end1",
      ...         "chrom2", "start2", "end2",
      ...         "count"]

  # Loop over chunks of interactions and progressively add them to "out.cool"
  In [7]: for df in pd.read_table("pixels.bg2", names=cols, chunksize=1_000_000):
    ...:      f.add_pixels(df)
    ...:

  # Important! If you forget to call f.finalize() the resulting .cool file will be empty
  In [8]: f.finalize()

  # Check that the resulting file has some interactions
  In [9]: htk.File("out.cool").attributes()["nnz"]
  Out[9]: 3118456


Ingesting interactions in a .hic file
-------------------------------------

Follow the same step as in the previous section and replace ``htk.cooler.File`` with ``htk.hic.File``.

Tips and tricks
---------------

When loading interactions into a .cool or .hic file, interactions are initially stored in a temporary file. When loading a large number of interactions, this temporary file can grow to be quite large. When this is the case, it is wise to pass a custom temporary folder where temporary files will be created:


.. code-block:: ipythonconsole

  In [1]: f = htk.cooler.FileWriter("out.cool", chroms, resolution=50_000, tmpdir="/var/tmp/hictk")

When ingesting interactions in a .hic file, performance can be improved by using multiple threads:

.. code-block:: ipythonconsole

  In [1]: f = htk.hic.FileWriter("out.hic", chroms, resolution=50_000, n_threads=8)

When memory allows it, it is possible to bypass temporary files by specifying a very large chunk size and ingesting all interactions at once. This can significantly speed up file creation.

.. code-block:: ipythonconsole

  # Initialize an empty .cool file

  In [1]: cols = ["chrom1", "start1", "end1",
      ...         "chrom2", "start2", "end2",
      ...         "count"]

  In [2]: df = pd.read_table("pixels.bg2", names=cols)

  In [3]: f = htk.cooler.FileWriter("out.cool", chroms, resolution=50_000, chunk_size=len(df) + 1)

  In [4]: f.add_pixels(df)

  In [5]: f.finalize()
