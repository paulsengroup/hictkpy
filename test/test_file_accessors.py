# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

import pytest

import hictkpy

from .helpers import pandas_avail

testdir = pathlib.Path(__file__).resolve().parent

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        ((testdir / "data" / "cooler_test_file.mcool::/resolutions/100000").as_posix(), None),
        (testdir / "data" / "cooler_test_file.mcool", 100_000),
        (testdir / "data" / "hic_test_file.hic", 100_000),
    ],
)


class TestClass:
    @pytest.mark.skipif(not pandas_avail(), reason="pandas is not available")
    def test_attributes(self, file, resolution):
        f = hictkpy.File(file, resolution)
        assert f.resolution() == 100_000
        assert f.nchroms() == 8
        assert f.nbins() == 1380

        assert "chr2L" in f.chromosomes()

        assert len(f.bins()) == 1380
        assert len(f.chromosomes()) == 8

        if f.is_cooler():
            assert f.attributes()["format"] == "HDF5::Cooler"
        else:
            assert f.attributes()["format"] == "HIC"

    def test_normalizations(self, file, resolution):
        f = hictkpy.File(file, resolution)

        cooler_weights = ["KR", "SCALE", "VC", "VC_SQRT", "weight"]
        hic_weights = ["ICE"]

        if f.is_cooler():
            assert f.avail_normalizations() == cooler_weights
        else:
            assert f.avail_normalizations() == hic_weights

        assert not f.has_normalization("foo")

        name = "weight" if f.is_cooler() else "ICE"
        weights = f.weights(name)
        assert len(weights) == f.nbins()

        if f.is_cooler():
            df = f.weights(cooler_weights)
            assert len(df.columns) == len(cooler_weights)
        else:
            df = f.weights(hic_weights)
            assert len(df.columns) == len(hic_weights)
        assert len(df) == f.nbins()
