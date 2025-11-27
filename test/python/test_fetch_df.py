# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import math

import pytest

import hictkpy

from .helpers import get_test_dir, numpy_avail, pandas_avail, pyarrow_avail

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (get_test_dir() / "data" / "cooler_test_file.mcool", 100_000),
        (get_test_dir() / "data" / "hic_test_file.hic", 100_000),
    ],
)


@pytest.mark.skipif(not pandas_avail() or not pyarrow_avail(), reason="either pandas or pyarrow are not available")
class TestFetchDF:
    def test_genome_wide(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            df = f.fetch().to_df()

        assert df["count"].sum() == 119_208_613
        assert len(df) == 890_384

    def test_cis(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            df = f.fetch("chr2L").to_df()
            assert df["count"].sum() == 19_968_156

            df = f.fetch("chr2R:10,000,000-15,000,000").to_df()
            assert df["count"].sum() == 4_519_080
            assert len(df.columns) == 3

            df = f.fetch("chr2R:10,000,000-15,000,000", join=True).to_df()
            assert df["count"].sum() == 4_519_080
            assert len(df.columns) == 7
            assert df["chrom1"].dtype.name == "category"
            assert df["chrom2"].dtype.name == "category"

            if numpy_avail():
                import numpy as np

                df = f.fetch("chr2R:10,000,000-15,000,000", count_type="int32").to_df()
                assert df["count"].dtype == np.int32

                df = f.fetch("chr2R:10,000,000-15,000,000", count_type=float).to_df()
                assert df["count"].dtype == np.float64

            df = f.fetch("chr2R\t10000000\t15000000", query_type="BED").to_df()
            assert len(df) == 1275

            df = f.fetch("chr2L:0-10,000,000", "chr2L:10,000,000-20,000,000").to_df()
            assert df["count"].sum() == 761_223
            assert len(df) == 9_999

            df = f.fetch("chr2L:0-10,000,000", "chr2L:0-15,000,000").to_df()
            assert df["count"].sum() == 9_270_385
            assert len(df) == 10_050

    def test_trans(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            df = f.fetch("chr2L", "chr2R").to_df()
            assert df["count"].sum() == 1_483_112

            df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000").to_df()
            assert df["count"].sum() == 83_604
            assert len(df.columns) == 3

            df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", join=True).to_df()
            assert df["count"].sum() == 83_604
            assert len(df.columns) == 7
            assert df["chrom1"].dtype.name == "category"
            assert df["chrom2"].dtype.name == "category"

            if numpy_avail():
                import numpy as np

                df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="int32").to_df()
                assert df["count"].dtype == np.int32

                df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type=float).to_df()
                assert df["count"].dtype == np.float64

            df = f.fetch("chr2R\t10000000\t15000000", "chrX\t0\t10000000", query_type="BED").to_df()
            assert len(df) == 4995

            df1 = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000").to_df()
            df2 = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000").to_df("lower_triangle")

            df2["bin1_id"], df2["bin2_id"] = df2["bin2_id"], df2["bin1_id"]
            df2.sort_values(by=["bin1_id", "bin2_id"], inplace=True)
            df2.reset_index(inplace=True, drop=True)
            assert df1.equals(df2)

    def test_balanced(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            if f.is_cooler():
                df = f.fetch("chr2R:10,000,000-15,000,000", normalization="weight").to_df()
            else:
                df = f.fetch("chr2R:10,000,000-15,000,000", normalization="ICE").to_df()

        assert math.isclose(59.349524704033215, df["count"].sum(), rel_tol=1.0e-5, abs_tol=1.0e-8)

    def test_diagonal_band(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            df = f.fetch(diagonal_band_width=100).to_df()

        assert len(df) == 129018
        assert df["count"].sum() == 104681024
        for i, dff in df.groupby("bin1_id"):
            assert len(dff) <= 100
