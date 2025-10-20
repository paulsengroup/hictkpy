# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

import pytest

import hictkpy

testdir = pathlib.Path(__file__).resolve().parent


class TestClass:
    pattern = "caught an attempt to access file .*, which has already been closed"

    def test_single_cell(self):
        path = testdir / "data" / "cooler_test_file.scool"
        with hictkpy.cooler.SingleCellFile(path) as f:
            assert f.path() == path

        with pytest.raises(RuntimeError, match=self.pattern):
            f.path()

        with hictkpy.cooler.SingleCellFile(path) as f:
            f.close()
            f.close()  # no-op

            with pytest.raises(RuntimeError, match=self.pattern):
                f.path()
