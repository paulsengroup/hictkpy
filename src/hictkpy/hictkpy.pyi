# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from . import (
    _hictkpy.cooler as _hictkpy.cooler,
    _hictkpy.hic as _hictkpy.hic
)
from ._hictkpy import (
    File as File,
    MultiResFile as MultiResFile,
    PixelSelector as PixelSelector,
    is_cooler as is_cooler,
    is_hic as is_hic,
    is_mcool_file as is_mcool_file,
    is_scool_file as is_scool_file
)


__all__: list = ['__doc__', 'File', 'MultiResFile', 'PixelSelector', 'is_cooler', 'is_mcool_file', 'is_scool_file', 'is_hic', 'cooler', 'hic', '__hictk_version__']

__hictk_version__: str
