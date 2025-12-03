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


In general, Pandas and PyArrow are required when hictkpy is returning data using :py:class:`pandas.DataFrame` or :py:class:`pyarrow.Table`, such as when calling :py:meth:`hictkpy.PixelSelector.to_pandas()` or :py:meth:`hictkpy.BinTable.to_df()`.

NumPy is required when calling methods returning data as :py:class:`numpy.array`, such as :py:meth:`hictkpy.PixelSelector.to_numpy()` or certain overloads of :py:meth:`hictkpy.BinTable.get_ids()`.

SciPy is required when fetching interactions as sparse matrix, such as :py:meth:`hictkpy.PixelSelector.to_coo()` and :py:meth:`hictkpy.PixelSelector.to_csr()`

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

Building hictkpy from source should not be necessary for regular users, as we publish pre-built wheels for Linux, macOS, and Windows for all Python versions we support (currently these include all CPython versions from 3.10 up until 3.14). For a complete and up-to-date list of available wheels refer to the `Download files <https://pypi.org/project/hictkpy/#files>`_ page on hictkpy's `homepage <https://pypi.org/project/hictkpy/>`_ on PyPI.

Building hictkpy's wheels from source requires a compiler toolchain supporting C++17, such as:

* GCC 8+
* Clang 8+
* Apple-Clang 10.0+
* MSVC 19.12+

Based on our testing, hictkpy's wheels compiled on Linux using Clang are noticeably faster than those compiled with GCC.
For this reason we recommend building hictkpy using a modern version of Clang whenever possible.
This can be achieved by redefining the ``CC`` and ``CXX`` environment variables before running pip (e.g. ``CC=clang CXX=clang++ pip install ...``).

Furthermore, the following tools are required:

* git 2.7+
* make or ninja


Installing the latest version from the main branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  pip install 'hictkpy[all] @ git+https://github.com/paulsengroup/hictkpy.git@main'

Installing version corresponding to a git tag
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  pip install 'hictkpy[all] @ git+https://github.com/paulsengroup/hictkpy.git@v1.4.0'

Installing from a release archive
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  pip install 'hictkpy[all] @ https://pypi.python.org/packages/source/h/hictkpy/hictkpy-1.4.0.tar.gz'

Running the automated tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^

When building hictkpy from source we highly recommend running the automated test suite before using hictkpy in production.

This can be achieved in several ways. Here is an example:

.. code-block:: bash

  git clone https://github.com/paulsengroup/hictkpy.git

  cd hictkpy

  # make sure to run tests for the same version/tag/commit used to build hictkpy
  git checkout v1.4.0

  # if you installed hictkpy in a venv make sure to install pytest in the venv
  pip install pytest

  pytest test/

**All tests are expected to pass. Do not ignore test failures!**

However, it is expected that some test cases will be skipped (especially if not all optional dependencies were installed)

Notes
^^^^^

Building hictkpy requires several dependencies that are not needed after the build process.
Some of these dependencies are installed using Conan, which creates several files under ``~/.conan2``.
If you don't need Conan for other purposes feel free to delete the ``~/.conan2`` once the build process completes successfully.

If you do not want to use Conan for dependency management you can set the ``HICTKPY_PROJECT_TOP_LEVEL_INCLUDES`` environment variable to an empty string.
See section ``[tool.scikit-build.cmake.define]`` in the `pyproject.toml <https://github.com/paulsengroup/hictkpy/blob/main/pyproject.toml>`__ file for the list of CMake variables that can be overridden by defining the appropriate environment variables.
