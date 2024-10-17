# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

import pytest

import hictkpy

testdir = pathlib.Path(__file__).resolve().parent


pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (testdir / "data" / "cooler_test_file.mcool", 100_000),
        (testdir / "data" / "hic_test_file.hic", 100_000),
    ],
)


class TestClass:
    def test_fetch_sum(self, file, resolution):
        f = hictkpy.File(file, resolution)
        assert f.fetch().sum() == 119_208_613
        assert f.fetch("chr2L").sum() == 19_968_156
        assert f.fetch("chr2L", "chr2R").sum() == 1_483_112
