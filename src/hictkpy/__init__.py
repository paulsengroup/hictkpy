# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

def _load_arrow_python_lib():
    import pyarrow


_load_arrow_python_lib()


from ._hictkpy import (
    __doc__,
    File,
    MultiResFile,
    PixelSelector,
    is_cooler,
    is_mcool_file,
    is_scool_file,
    is_hic,
    cooler,
    hic,
    __hictk_version__,
)
from importlib.metadata import version

__version__ = version("hictkpy")
__all__ = [
    "__doc__",
    "File",
    "MultiResFile",
    "PixelSelector",
    "is_cooler",
    "is_mcool_file",
    "is_scool_file",
    "is_hic",
    "cooler",
    "hic",
    "__hictk_version__",
]
