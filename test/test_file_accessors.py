# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os

import pytest

import hictkpy

testdir = os.path.dirname(os.path.abspath(__file__))

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (os.path.join(testdir, "data", "cooler_test_file.cool"), 100_000),
        (os.path.join(testdir, "data", "hic_test_file.hic"), 100_000),
    ],
)


class TestClass:
    def test_attributes_cooler(self, file, resolution):
        f = hictkpy.File(file, resolution)
        assert f.bin_size() == 100_000
        assert f.nbins() == 1380

        assert "chr2L" in f.chromosomes()
        assert len(f.bins()) == 1380
        assert len(f.chromosomes()) == 8

        if f.is_cooler():
            assert f.attributes()["format"] == "HDF5::Cooler"
        else:
            assert f.attributes()["format"] == "HIC"
