<!--
Copyright (C) 2025 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# Documentation README

## How to build hictkpy's documentation

The instructions in this README assume all commands are being run from the root of hictkpy's repository.

```bash
venv/bin/pip install --upgrade pip  # --group option requires a modern version of pip
venv/bin/pip install . --group docs -v

# Activate venv
. venv/bin/activate

# Clean old build files (optional)
make -C docs clean

make -C docs linkcheck html latexpdf
```

Open the HTML documentation:

```bash
# Linux
xdg-open docs/_build/html/index.html

# macOS
open docs/_build/html/index.html
```

Open the PDF documentation:

```bash
# Linux
xdg-open docs/_build/latex/hictkpy.pdf

# macOS
open docs/_build/latex/hictkpy.pdf
```
