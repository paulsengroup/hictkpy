# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


from ._hictkpy import __doc__, File, PixelSelector, is_cooler, is_hic, cooler, __hictk_version__
from importlib.metadata import version

__version__ = version("hictkpy")
__all__ = ["__doc__", "File", "PixelSelector", "is_cooler", "is_hic", "cooler", "__hictk_version__"]
