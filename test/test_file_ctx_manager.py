# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

import pytest

import hictkpy

testdir = pathlib.Path(__file__).resolve().parent


class TestFileCTXManager:
    pattern = "caught an attempt to access file .*, which has already been closed"

    def test_single_resolution(self):
        path = testdir / "data" / "cooler_test_file.mcool"
        resolution = 100_000
        with hictkpy.File(path, resolution) as f:
            assert f.path() == path

        with pytest.raises(RuntimeError, match=self.pattern):
            f.uri()

        with hictkpy.File(path, resolution) as f:
            f.close()
            f.close()  # no-op

            with pytest.raises(RuntimeError, match=self.pattern):
                f.path()
