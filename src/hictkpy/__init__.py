# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


def _get_hictkpy_version() -> str:
    from importlib.metadata import version

    return version("hictkpy")


from ._hictkpy import (
    Bin,
    BinTable,
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

__version__ = _get_hictkpy_version()
__all__ = [
    "__doc__",
    "Bin",
    "BinTable",
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
