..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Introduction
============

hictkpy provides Python bindings to `hictk <https://github.com/paulsengroup/hictk>`_, a blazing fast toolkit to work with .hic and .cool files.

.. only:: not latex

   Documentation formats
   ---------------------

   You are reading the HTML version of the documentation. An alternative `PDF
   version <https://hictkpy.readthedocs.io/_/downloads/en/stable/pdf/>`_ is
   also available.

   Installation
   ------------

.. only:: latex

   Documentation formats
   ---------------------

   You are reading the PDF version of the documentation.

   The live HTML version of the documentation is available at `<https://hictkpy.readthedocs.io/en/stable/>`_.

   .. rubric:: Installation

hictkpy can be installed using pip or conda with e.g., ``pip install 'hictkpy[all]'``.
Refer to :doc:`Installation <./installation>` for more details.

.. only:: not latex

   How to cite this project?
   -------------------------

.. only:: latex

   .. rubric:: How to cite this project?

Please use the following BibTeX template to cite hictkpy in scientific
discourse:

.. code-block:: bibtex

  @article{hictk,
      author = {Rossini, Roberto and Paulsen, Jonas},
      title = "{hictk: blazing fast toolkit to work with .hic and .cool files}",
      journal = {Bioinformatics},
      volume = {40},
      number = {7},
      pages = {btae408},
      year = {2024},
      month = {06},
      issn = {1367-4811},
      doi = {10.1093/bioinformatics/btae408},
      url = {https://doi.org/10.1093/bioinformatics/btae408},
      eprint = {https://academic.oup.com/bioinformatics/article-pdf/40/7/btae408/58385157/btae408.pdf},
  }

.. only:: not latex

   Table of contents
   -----------------

.. toctree::
   :maxdepth: 1

   installation
   quickstart
   creating_cool_hic_files
   multithreading_and_multiprocessing

.. toctree::
   :caption: API Reference
   :maxdepth: 1

   api/index
