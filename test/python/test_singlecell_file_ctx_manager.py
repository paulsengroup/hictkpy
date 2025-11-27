# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pytest

import hictkpy

from .helpers import get_test_dir


class TestSingleCellFileCTXManager:
    pattern = "caught an attempt to access file .*, which has already been closed"

    def test_single_cell(self):
        path = get_test_dir() / "data" / "cooler_test_file.scool"
        with hictkpy.cooler.SingleCellFile(path) as f:
            assert f.path() == path

        with pytest.raises(RuntimeError, match=self.pattern):
            f.path()

        with hictkpy.cooler.SingleCellFile(path) as f:
            f.close()
            f.close()  # no-op

            with pytest.raises(RuntimeError, match=self.pattern):
                f.path()
