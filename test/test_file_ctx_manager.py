# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

import pytest

import hictkpy

testdir = pathlib.Path(__file__).resolve().parent


class TestClass:
    pattern = "caught an attempt to access file .*, which has already been closed"

    def test_single_resolution(self):
        uri = testdir / "data" / "cooler_test_file.mcool::/resolutions/100000"
        with hictkpy.File(uri) as f:
            assert f.uri() == str(uri)

        with pytest.raises(RuntimeError, match=self.pattern):
            f.uri()

        with hictkpy.File(uri) as f:
            f.close()
            f.close()  # no-op

            with pytest.raises(RuntimeError, match=self.pattern):
                f.uri()
