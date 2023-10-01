# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os

import numpy as np
import pytest

import hictkpy

testdir = os.path.dirname(os.path.abspath(__file__))

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (os.path.join(testdir, "data", "cooler_test_file.mcool"), 100_000),
        (os.path.join(testdir, "data", "hic_test_file.hic"), 100_000),
    ],
)


class TestClass:
    def test_genome_wide(self, file, resolution):
        f = hictkpy.File(file, resolution)

        df = f.fetch().to_df()
        assert df["count"].sum() == 119_208_613
        assert len(df) == 890_384

    def test_cis(self, file, resolution):
        f = hictkpy.File(file, resolution)

        df = f.fetch("chr2R:10,000,000-15,000,000").to_df()
        assert df["count"].sum() == 4_519_080
        assert len(df.columns) == 3

        df = f.fetch("chr2R:10,000,000-15,000,000", join=True).to_df()
        assert df["count"].sum() == 4_519_080
        assert len(df.columns) == 7

        df = f.fetch("chr2R:10,000,000-15,000,000", count_type="int").to_df()
        assert df["count"].dtype == np.int32

        df = f.fetch("chr2R:10,000,000-15,000,000", count_type="float").to_df()
        assert df["count"].dtype == np.float64

        df = f.fetch("chr2R\t10000000\t15000000", query_type="BED").to_df()
        assert len(df) == 1275

    def test_trans(self, file, resolution):
        f = hictkpy.File(file, resolution)

        df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000").to_df()
        assert df["count"].sum() == 83_604
        assert len(df.columns) == 3

        df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", join=True).to_df()
        assert df["count"].sum() == 83_604
        assert len(df.columns) == 7

        df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="int").to_df()
        assert df["count"].dtype == np.int32

        df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="float").to_df()
        assert df["count"].dtype == np.float64

        df = f.fetch("chr2R\t10000000\t15000000", "chrX\t0\t10000000", query_type="BED").to_df()
        assert len(df) == 4995

    def test_balanced(self, file, resolution):
        f = hictkpy.File(file, resolution)

        if f.is_cooler():
            df = f.fetch("chr2R:10,000,000-15,000,000", normalization="weight").to_df()
        else:
            df = f.fetch("chr2R:10,000,000-15,000,000", normalization="ICE").to_df()

        assert np.isclose(59.349524704033215, df["count"].sum())
