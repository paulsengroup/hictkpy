..
   Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Creating .cool and .hic files
#############################

hictkpy supports creating .cool and .hic files from pre-binned interactions in COO or BedGraph2 format.

The examples in this section use file `4DNFIOTPSS3L.hic <https://data.4dnucleome.org/files-processed/4DNFIOTPSS3L>`_,
which can be downloaded from the 4D Nucleome Data Portal
`here <https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/7386f953-8da9-47b0-acb2-931cba810544/4DNFIOTPSS3L.hic>`_.

Preparation
-----------

The first step involves converting interactions from ``4DNFIOTPSS3L.hic`` to bedGraph2 format.
This can be achieved using ``hictk dump`` (or alternatively with :py:meth:`hictkpy.File.fetch()`).

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

  user@dev:/tmp$ head chrom.sizes

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
     ...:            .set_index("name")["length"]
     ...:            .to_dict()

  In [4]: chroms
  Out[4]:
  {'2L': 23513712,
   '2R': 25286936,
   '3L': 28110227,
   '3R': 32079331,
   '4': 1348131,
   'X': 23542271,
   'Y': 3667352}

  # Define the name of the columns for later use
  In [5]: cols = ["chrom1", "start1", "end1",
     ...:         "chrom2", "start2", "end2",
     ...:         "count"]

  # Initialize an empty .cool file
  In [6]: with htk.cooler.FileWriter("out.cool", chroms, resolution=50_000) as writer:
     ...:     # Lazily load pixels in chunks to reduce memory usage
     ...:     pixels = pd.read_table("pixels.bg2", names=cols, chunksize=1_000_000)
     ...:     # Add chunks of pixels one by one
     ...:     for i, df in enumerate(pixels):
     ...:         print(f"adding chunk #{i}...")
     ...:         writer.add_pixels(df)
     ...:
  adding chunk #0...
  adding chunk #1...
  adding chunk #2...
  adding chunk #3...

  # Check that the resulting file has some interactions
  In [7]: htk.File("out.cool").attributes()["nnz"]
  Out[7]: 3118456


Ingesting interactions in a .hic file
-------------------------------------

Follow the same steps as above for ``.cool`` files, but replace ``htk.cooler.FileWriter`` with ``htk.hic.FileWriter``.

Tips and tricks
---------------

When loading interactions into a .cool or .hic file, interactions are initially stored in a temporary file.
For a large number of interactions, this temporary file can become quite large.
In such cases, it may be appropriate to pass a custom temporary folder where these files will be created:

.. code-block:: ipythonconsole

  In [1]: f = htk.cooler.FileWriter("out.cool", chroms, resolution=50_000, tmpdir="/var/tmp/hictk")

When ingesting interactions in a .hic file, performance can be improved by using multiple threads:

.. code-block:: ipythonconsole

  In [1]: f = htk.hic.FileWriter("out.hic", chroms, resolution=50_000, n_threads=8)

When memory allows, it is possible to bypass temporary file creation by specifying a very large chunk size and ingesting all interactions at once.
This can significantly speed up file creation:

.. code-block:: ipythonconsole

  # Initialize an empty .cool file

  In [1]: cols = ["chrom1", "start1", "end1",
     ...:         "chrom2", "start2", "end2",
     ...:         "count"]

  In [2]: df = pd.read_table("pixels.bg2", names=cols)

  In [3]: with htk.cooler.FileWriter("out.cool", chroms, resolution=50_000, chunk_size=len(df) + 1) as writer:
     ...:     writer.add_pixels(df)
     ...:

In case it is not possible to install a compatible version of pandas or pyarrow, the ``FileWriter``
classes support ingesting interactions from dictionaries of iterables
(e.g., a dictionary mapping keys ``bin1_id``, ``bin2_id``, and ``count`` to iterables yielding numbers of the appropriate type).

For more details, refer to the documentation for
:py:meth:`hictkpy.cooler.FileWriter.add_pixels_from_dict()`
and
:py:meth:`hictkpy.hic.FileWriter.add_pixels_from_dict()`.
