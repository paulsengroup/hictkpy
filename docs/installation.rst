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


In general, pandas and pyarrow are required when hictkpy is returning data using pandas.DataFrame or arrow.Table, such as when calling :py:meth:`hictkpy.PixelSelector.to_pandas()` or :py:meth:`hictkpy.BinTable.to_df()`.

NumPy is required when calling methods returning data as np.array, such as :py:meth:`hictkpy.PixelSelector.to_numpy()` or certain overloads of :py:meth:`hictkpy.BinTable.get_ids()`.

SciPy is required when fetching interactions as sparse matrix with, such as :py:meth:`hictkpy.PixelSelector.to_coo()` and :py:meth:`hictkpy.PixelSelector.to_csr()`

In case applications call methods depending on missing third-party dependencies, hictkpy will raise an exception like the following:

.. code-block:: ipythonconsole

  In [3] f.fetch().to_numpy()

  ModuleNotFoundError: No module named 'numpy'

  The above exception was the direct cause of the following exception:

  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
  ModuleNotFoundError: To enable numpy support, please install numpy with: pip install 'hictkpy[numpy]'
  Alternatively, you can install hictkpy with all its dependencies by running: pip install 'hictkpy[all]'

Conda (bioconda)
----------------

.. code-block:: bash

  conda install -c conda-forge -c bioconda hictkpy

From source
-----------

Building hictkpy from source should not be necessary for regular users, as we publish pre-built wheels for Linux, MacOS, and Windows for all Python versions we support (currently these include all CPython versions from 3.9 up until 3.13). For a complete and up-to-date list of available wheels refer to the `Download files <https://pypi.org/project/hictkpy/#files>`_ page on hictkpy's `homepage <https://pypi.org/project/hictkpy/>`_ on PyPI.

Building hictkpy's wheels from source requires a compiler toolchain supporting C++17, such as:

* GCC 8+
* Clang 8+
* Apple-Clang 10.0+
* MSVC 19.12+

Furthermore, the following tools are required:

* CMake 3.25+
* git 2.7+
* make or ninja


Installing the latest version from the main branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  pip install 'hictkpy[all] @ git+https://github.com/paulsengroup/hictkpy.git@main'

Installing version corresponding to a git tag
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  pip install 'hictkpy[all] @ git+https://github.com/paulsengroup/hictkpy.git@v1.2.0'

Installing from a release archive
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  pip install 'hictkpy[all] @ https://pypi.python.org/packages/source/h/hictkpy/hictkpy-1.2.0.tar.gz'

Running the automated tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^

When building hictkpy from source we highly recommend running the automated test suite before using hictkpy in production.

This can be achieved in several ways. Here is an example:

.. code-block:: bash

  git clone https://github.com/paulsengroup/hictkpy.git

  # make sure to run tests for the same version/tag/commit used to build hictkpy
  git checkout v1.2.0

  python -m venv venv

  # On Windows use venv\Scripts\pip.exe instead
  venv/bin/pip install pytest

  venv/bin/pytest test/

**All tests are expected to pass. Do not ignore test failures!**

However, it is expected that some test cases are skipped (especially if not all optional dependencies where installed).

Notes
^^^^^

Building hictkpy requires several dependencies that are not needed after the build process.
Some of these dependencies are installed using Conan, which creates several files under ``~/.conan2``. if you don't need Conan for other purposes feel free to delete the ``~/.conan2`` once the build process completes successfully.
