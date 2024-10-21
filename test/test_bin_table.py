# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


import pytest

import hictkpy

from .helpers import numpy_avail, pandas_avail, pyarrow_avail


class TestClass:
    def test_accessors(self):
        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)

        assert len(bins.chromosomes()) == 2
        assert bins.resolution() == 100
        assert bins.type() == "fixed"
        assert len(bins) == 15

    def test_getters(self):
        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)

        assert bins.get(1).chrom == "chr1"
        assert bins.get(1).start == 100
        assert bins.get(1).end == 200
        with pytest.raises(Exception):
            bins.get(9999)

        assert bins.get_id("chr1", 153) == 1
        with pytest.raises(Exception):
            bins.get_id("abc", 100)

    @pytest.mark.skipif(not numpy_avail(), reason="numpy is not available")
    def test_vectorized_getters(self):
        import numpy as np

        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)

        assert len(bins.get(np.array([1, 1]))) == 2
        assert len(bins.get_ids(np.array(["chr1", "chr1"]), np.array([1, 1]))) == 2

    @pytest.mark.skipif(not pandas_avail() or not pyarrow_avail(), reason="pandas is not available")
    def test_merge(self):
        import pandas as pd

        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)

        expected = pd.DataFrame(
            {
                "chrom1": pd.Categorical(["chr1", "chr2"], categories=["chr1", "chr2"]),
                "start1": [0, 0],
                "end1": [100, 100],
                "chrom2": pd.Categorical(["chr1", "chr2"], categories=["chr1", "chr2"]),
                "start2": [0, 0],
                "end2": [100, 100],
            }
        )

        bin_ids = pd.DataFrame({"bin1_id": [0, 10], "bin2_id": [0, 10]})

        assert (expected == bins.merge(bin_ids).drop(columns=["bin1_id", "bin2_id"])).all().all()

    def test_to_df(self):
        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)

        assert len(bins.to_df()) == len(bins)
        assert len(bins.to_df("chr1")) == 10
        assert len(bins.to_df("chr2:0-200")) == 2
        assert len(bins.to_df("chr2\t0\t200", "BED")) == 2
        with pytest.raises(RuntimeError):
            bins.to_df("chr0")
