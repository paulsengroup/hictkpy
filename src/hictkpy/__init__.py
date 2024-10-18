# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


def _load_pyarrow_and_check_abi_compat():
    import pyarrow as pa

    from ._hictkpy import __hictkpy_arrow_version__

    major, minor, patch = __hictkpy_arrow_version__

    if not pa.__version__.startswith(f"{major}.{minor}"):
        raise ImportError(
            "Detected Arrow ABI version mismatch!\n"
            f"hictkpy was compiled with Arrow v{major}.{minor}.{patch}, which is not ABI compatible with the currently "
            f"installed version of pyarrow (v{pa.__version__}).\n"
            'Please install a compatible version of pyarrow with e.g. "pip install '
            f'pyarrow=={major}.{minor}".'
        )


_load_pyarrow_and_check_abi_compat()


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
