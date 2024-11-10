# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
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
    def test_resolution_mismatch(self, file, resolution):
        f = hictkpy.MultiResFile(file)

        if not f[resolution].is_cooler():
            pytest.skip(f'File "{file}" is not in .mcool format')

        assert 1_000_000 in f.resolutions()
        with pytest.raises(RuntimeError, match="resolution mismatch"):
            hictkpy.File(f"{file}::/resolutions/{resolution}", 1_000_000)

    def test_missing_resolution_param(self, file, resolution):
        f = hictkpy.MultiResFile(file)

        if f[resolution].is_cooler():
            ext = ".mcool"
        else:
            ext = ".hic"

        with pytest.raises(RuntimeError, match=f"resolution is required and cannot be None when opening {ext} files"):
            hictkpy.File(file)
