<!--
Copyright (C) 2023 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# hictkpy

[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)
[![CI](https://github.com/paulsengroup/hictkpy/actions/workflows/pip.yml/badge.svg)](https://github.com/paulsengroup/hictkpy/actions/workflows/pip.yml)
[![Download from Bioconda](https://img.shields.io/conda/vn/bioconda/hictkpy?label=bioconda&logo=Anaconda)](https://anaconda.org/bioconda/hictkpy)
[![docs](https://readthedocs.org/projects/hictkpy/badge/?version=latest)](https://hictkpy.readthedocs.io/en/latest/?badge=latest)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8220299.svg)](https://doi.org/10.5281/zenodo.8220299)

---

Python bindings for hictk, a blazing fast toolkit to work with .hic and .cool files.

## Installing hictkpy

hictkpy can be installed in various ways. The simplest method is using pip: `pip install hictkpy[all]`.

Refer to [Installation](https://hictkpy.readthedocs.io/en/stable/installation.html) for alternative methods.

## Using hictkpy

```python3
import hictkpy

path_to_clr = "file.mcool"  # "file.hic"

clr = hictkpy.File(path_to_clr, 100_000)
sel = clr.fetch("chr1")

df = sel.to_df()     # Get interactions as a pd.DataFrame
m1 = sel.to_numpy()  # Get interactions as a numpy matrix
m2 = sel.to_csr()    # Get interactions as a scipy.sparse.csr_matrix
```

For more detailed examples refer to [Quickstart](https://hictkpy.readthedocs.io/en/stable/quickstart.html).

The complete documentation for hictkpy API is available [here](https://hictkpy.readthedocs.io/en/stable/api/index.html).

## Citing

If you use hictkpy in you research, please cite the following publication:

Roberto Rossini, Jonas Paulsen, hictk: blazing fast toolkit to work with .hic and .cool files
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
