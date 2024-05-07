..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Introduction
============

hictkpy provides Python bindings to `hictk <https://github/paulsengroup/hictk>`_, a blazing fast toolkit to work with .hic and .cool files.

.. only:: not latex

   Documentation formats
   ---------------------

   You are reading the HTML version of the documentation. An alternative `PDF
   version <https://hictkpy.readthedocs.io/_/downloads/en/latest/pdf/>`__ is
   also available.

   Installation
   ------------

.. only:: latex

   .. rubric:: Installation

Python bindings for hictk can be installed using pip or conda. See :doc:`here <./installation>` for more details.

.. only:: not latex

   How to cite this project?
   -------------------------

.. only:: latex

   .. rubric:: How to cite this project?

Please use the following BibTeX template to cite hictkpy in scientific
discourse:

.. code-block:: bibtex

  @article {hictk,
	  author = {Roberto Rossini and Jonas Paulsen},
	  title = {hictk: blazing fast toolkit to work with .hic and .cool files},
	  elocation-id = {2023.11.26.568707},
	  year = {2023},
	  doi = {10.1101/2023.11.26.568707},
	  publisher = {Cold Spring Harbor Laboratory},
	  URL = {https://www.biorxiv.org/content/early/2023/11/27/2023.11.26.568707},
	  eprint = {https://www.biorxiv.org/content/early/2023/11/27/2023.11.26.568707.full.pdf},
	  journal = {bioRxiv}
  }

.. only:: not latex

   Table of contents
   -----------------

.. toctree::
   :maxdepth: 1

   installation
   quickstart
   creating_cool_hic_files

.. toctree::
   :caption: API Reference
   :maxdepth: 1

   api/index
