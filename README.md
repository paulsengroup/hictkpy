<!--
Copyright (C) 2023 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# hictkpy

[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)
[![Download from Bioconda](https://img.shields.io/conda/vn/bioconda/hictkpy?label=bioconda&logo=Anaconda)](https://anaconda.org/bioconda/hictkpy)

<!-- [![Zenodo DOI]()]() -->
---

Python bindings for hictk, a blazing fast toolkit to work with .hic and .cool files.

## Installing hictkpy

hictkpy can be installed in various ways.

### PIP

`pip install hictkpy`

### Conda (bioconda)

`conda install -c conda-forge -c bioconda hictkpy`

### From source

`pip install 'git+https://github.com/paulsengroup/hictkpy.git@main'`

On Windows you will have to manually install some of hictk dependencies, namely hdf5 (with zlib support) and libdeflate.

## Using hictkpy

```python3
import hictkpy

path_to_clr = "file.mcool"  # "file.hic"

clr = hictkpy.File(path_to_clr, 100_000)
sel = clr.fetch("chr1")

df = sel.to_df()     # Get interactions as a pd.DataFrame
m1 = sel.to_numpy()  # Get interactions as a numpy matrix
m2 = sel.to_coo()    # Get interactions as a scipy.sparse.coo_matrix

# Loop over interactions
for bin1_id, bin2_id, count in clr.fetch("chr1"):
  print(bin1_id, ...)

# Loop over interactions
for chrom1, start1, end1, chrom2, start2, end2, count in clr.fetch("chr1", join=True):
  print(chrom1, ...)

# Fetch interactions using UCSC queries
clr.fetch("chr1:0-10,000,000").to_df()
clr.fetch("chr1:0-10,000,000", "chr2:100,000,000-105,000,000").to_df()

# Fetch interactions using BED queries
clr.fetch("chr1\t0\t10000000", query_type="BED").to_df()

# Fetch balanced interactions
clr.fetch("chr1", normalization="weight").to_df()
clr.fetch("chr1", normalization="VC").to_df()

# Sum interactions overlapping query
clr.fetch("chr1").sum()

# Count non-zero entries
clr.fetch("chr1").nnz()
```
