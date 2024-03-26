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

Note that this will install hictk's build dependencies under ``~/.conan2``, if you don't need Conan for other purposes feel free to delete this ``~/.conan2`` after installing hictkpy from git.
