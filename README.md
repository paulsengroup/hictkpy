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

hictkpy can be installed in various ways. The simples method is using pip: `pip install hictkpy`.

Refer to [Installation](https://hictkpy.readthedocs.io/en/latest/installation.html) for alternative methods.

## Using hictkpy

```python3
import hictkpy

path_to_clr = "file.mcool"  # "file.hic"

clr = hictkpy.File(path_to_clr, 100_000)
sel = clr.fetch("chr1")

df = sel.to_df()     # Get interactions as a pd.DataFrame
m1 = sel.to_numpy()  # Get interactions as a numpy matrix
m2 = sel.to_coo()    # Get interactions as a scipy.sparse.coo_matrix
```

For more detailed examples refer to [Quickstart](https://hictkpy.readthedocs.io/en/latest/quickstart.html).

The complete documentation for hictkpy API is available [here](https://hictkpy.readthedocs.io/en/latest/hictkpy.html).

## Citing

If you use hictkpy in you reaserch, please cite the following publication:

Preprint available soon.

<details>
<summary>BibTex</summary>

```bibtex
@misc{hictk,
author = {Roberto Rossini},
year = {2023},
note = {https://github.com/paulsengroup/hictk},
title = {hictk: blazing fast toolkit to work with .hic and .cool files}
}
```

</details>
