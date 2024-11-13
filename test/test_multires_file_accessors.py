# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
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
    ],
)


class TestClass:
    def test_attributes(self, file):
        f = hictkpy.MultiResFile(file)

        assert f.path() == file
        assert (f.resolutions() == [100_000, 1_000_000]).all()
        assert len(f.chromosomes()) == 8

        assert f[100_000].resolution() == 100_000
