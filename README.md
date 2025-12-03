<!--
Copyright (C) 2023 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# hictkpy

---

<!-- markdownlint-disable MD033 -->

<table>
    <tr>
      <td>Downloads</td>
      <td>
        <a href="https://pypi.org/project/hictkpy/">
          <img src="https://img.shields.io/pypi/v/hictkpy" alt="PyPI">
        </a>
        &nbsp
        <a href="https://anaconda.org/bioconda/hictkpy">
          <img src="https://img.shields.io/conda/vn/bioconda/hictkpy?label=bioconda&logo=Anaconda" alt="Bioconda">
        </a>
        &nbsp
        <a href="https://doi.org/10.5281/zenodo.8220299">
          <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.8220299.svg" alt="Zenodo">
        </a>
      </td>
    </tr>
    <tr>
      <td>Documentation</td>
      <td>
        <a href="https://hictkpy.readthedocs.io/">
          <img src="https://readthedocs.org/projects/hictkpy/badge/?version=latest" alt="Documentation">
        </a>
      </td>
    </tr>
    <tr>
      <td>License</td>
      <td>
        <a href="https://github.com/paulsengroup/hictkpy/blob/main/LICENSE">
          <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
        </a>
      </td>
    </tr>
    <tr>
      <td>CI</td>
      <td>
        <a href="https://github.com/paulsengroup/hictkpy/actions/workflows/ci.yml">
          <img src="https://github.com/paulsengroup/hictkpy/actions/workflows/ci.yml/badge.svg" alt="CI Status">
        </a>
        &nbsp
        <a href="https://github.com/paulsengroup/hictkpy/actions/workflows/build-wheels.yml">
          <img src="https://github.com/paulsengroup/hictkpy/actions/workflows/build-wheels.yml/badge.svg" alt="Build wheels Status">
        </a>
      </td>
    </tr>
    <tr>
      <td>CodeQL</td>
      <td>
        <a href="https://github.com/paulsengroup/hictkpy/actions/workflows/codeql-cpp.yml">
          <img src="https://github.com/paulsengroup/hictkpy/actions/workflows/codeql-cpp.yml/badge.svg" alt="CodeQL (C++) Status">
        </a>
        &nbsp
        <a href="https://github.com/paulsengroup/hictkpy/actions/workflows/codeql-python.yml">
          <img src="https://github.com/paulsengroup/hictkpy/actions/workflows/codeql-python.yml/badge.svg" alt="CodeQL (Python) Status">
        </a>
        &nbsp
        <a href="https://github.com/paulsengroup/hictkpy/actions/workflows/codeql-actions.yml">
          <img src="https://github.com/paulsengroup/hictkpy/actions/workflows/codeql-actions.yml/badge.svg" alt="CodeQL (GH Actions) Status">
        </a>
      </td>
    </tr>
    <tr>
      <td>Fuzzy Testing</td>
      <td>
        <a href="https://github.com/paulsengroup/hictkpy/actions/workflows/fuzzy-testing.yml">
          <img src="https://github.com/paulsengroup/hictkpy/actions/workflows/fuzzy-testing.yml/badge.svg" alt="Fuzzy Testing Status">
        </a>
      </td>
    </tr>
</table>

<!-- markdownlint-enable MD033 -->

---

Python bindings for [hictk](https://github.com/paulsengroup/hictk), a blazing fast toolkit to work with .hic and .cool files.

If you are looking for the R API, checkout the [hictkR](https://github.com/paulsengroup/hictkR) repository.

## Installing hictkpy

hictkpy can be installed in various ways.
The simplest method is using pip: `pip install 'hictkpy[all]'`.

Refer to [Installation](https://hictkpy.readthedocs.io/en/stable/installation.html) for alternative methods.

## Using hictkpy

```python3
import hictkpy

path_to_clr = "file.mcool"  # "file.hic"

clr = hictkpy.File(path_to_clr, 100_000)
sel = clr.fetch("chr1")

df = sel.to_df()     # Get interactions as a pandas.DataFrame
m1 = sel.to_numpy()  # Get interactions as a numpy matrix
m2 = sel.to_csr()    # Get interactions as a scipy.sparse.csr_matrix
```

For more detailed examples refer to the [Quickstart](https://hictkpy.readthedocs.io/en/stable/quickstart.html) section in the documentation.

<!-- markdownlint-disable MD059 -->

The complete documentation for the hictkpy API is available [here](https://hictkpy.readthedocs.io/en/stable/api/index.html).

<!-- markdownlint-enable MD059 -->

## Citing

If you use hictkpy in your research, please cite the following publication:

Roberto Rossini, Jonas Paulsen, hictk: blazing fast toolkit to work with .hic and .cool files\
_Bioinformatics_, Volume 40, Issue 7, July 2024, btae408, [https://doi.org/10.1093/bioinformatics/btae408](https://doi.org/10.1093/bioinformatics/btae408)

<details>
<summary>BibTex</summary>

```bibtex
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
```

</details>
