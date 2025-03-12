# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

import pytest

import hictkpy

testdir = pathlib.Path(__file__).resolve().parent

pytestmark = pytest.mark.parametrize(
    "file",
    [
        testdir / "data" / "cooler_test_file.mcool",
        testdir / "data" / "hic_test_file.hic",
    ],
)


class TestClass:
    def test_negative_resolution(self, file):
        with pytest.raises(ValueError, match="resolution cannot be negative"):
            hictkpy.File(file, -1)
