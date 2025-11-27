# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

import pytest

import hictkpy

testdir = pathlib.Path(__file__).resolve().parent


class TestMultiResFileCTXManager:
    pattern = "caught an attempt to access file .*, which has already been closed"

    def test_multi_resolution(self):
        path = testdir / "data" / "cooler_test_file.mcool"
        with hictkpy.MultiResFile(path) as f:
            assert f.path() == path

        with pytest.raises(RuntimeError, match=self.pattern):
            f.path()

        with hictkpy.MultiResFile(path) as f:
            f.close()
            f.close()  # no-op

            with pytest.raises(RuntimeError, match=self.pattern):
                f.path()
