..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Installation
############

hictkpy can be installed in various ways.

PIP
---

.. code-block:: bash

  pip install 'hictkpy[all]'

This will install hictkpy together with all its third-party dependencies.

It is also possible to install hictkpy with a minimal set of dependencies with one of the following commands:

.. code-block:: bash

  pip install hictkpy  # this target has no runtime dependencies!
  pip install 'hictkpy[numpy]'
  pip install 'hictkpy[pandas]'
  pip install 'hictkpy[pyarrow]'
  pip install 'hictkpy[scipy]'


Pandas is required when calling :py:meth:`hictkpy.PixelSelector.to_df()`, :py:meth:`hictkpy.PixelSelector.to_pandas()`, and :py:meth:`hictkpy.File.bins()`.

SciPy is required when fetching interactions as sparse matrix with e.g. :py:meth:`hictkpy.PixelSelector.to_coo()` ot :py:meth:`hictkpy.PixelSelector.to_csr()`.


Conda (bioconda)
----------------

.. code-block:: bash

  conda install -c conda-forge -c bioconda hictkpy

From source
-----------

.. code-block:: bash

  pip install 'git+https://github.com/paulsengroup/hictkpy.git@main'

Note that this will install hictk's build dependencies under ``~/.conan2``, if you don't need Conan for other purposes feel free to delete th ``~/.conan2`` folder after installing hictkpy from git.
