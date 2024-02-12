# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


from .hictkpy import File, PixelSelector, is_cooler, is_hic, cooler, __hictk_version__

try:
    from importlib.metadata import version
except ModuleNotFoundError:
    from importlib_metadata import version

__version__ = version("hictkpy")
