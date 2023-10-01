..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Installation
############

hictkpy can be installed in various ways.

PIP
---

.. code-block:: bash

  pip install hictkpy



Conda (bioconda)
----------------

.. code-block:: bash

  conda install -c conda-forge -c bioconda hictkpy

From source
-----------

.. code-block:: bash

  pip install 'git+https://github.com/paulsengroup/hictkpy.git@main'


On Windows you will have to manually install some of hictk dependencies, namely hdf5 (with zlib support) and libdeflate.
