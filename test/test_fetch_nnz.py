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
    def test_fetch_nnz(self, file, resolution):
        f = hictkpy.File(file, resolution)
        assert f.fetch().nnz() == 890_384
        assert f.fetch("chr2R").nnz() == 31_900
        assert f.fetch("chr2L", "chr2R").nnz() == 59_124
