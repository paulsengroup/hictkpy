# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


from importlib.metadata import version

from ._hictkpy import (
    File,
    MultiResFile,
    PixelSelector,
    __doc__,
    __hictk_version__,
    cooler,
    hic,
    is_cooler,
    is_hic,
    is_mcool_file,
    is_scool_file,
)

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
