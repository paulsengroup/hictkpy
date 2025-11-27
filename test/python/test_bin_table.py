# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


import itertools

import pytest

import hictkpy

from .helpers import numpy_avail, pandas_avail, pyarrow_avail


class TestBinTable:
    def test_ctor_fixed_bins(self):
        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)
        assert bins.type() == "fixed"
        assert len(bins) == 15

    @pytest.mark.skipif(not pandas_avail(), reason="pandas is not available")
    def test_ctor_variable_bins(self):
        import pandas as pd

        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2", "chr2", "chr3"],
                "start": [0, 100, 0, 15, 0],
                "end": [100, 127, 15, 75, 12],
            }
        )

        bins = hictkpy.BinTable(data)
        assert bins.type() == "variable"
        assert len(bins) == len(data)

    def test_accessors(self):
        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)

        assert len(bins.chromosomes()) == 2
        assert bins.resolution() == 100
        assert bins.type() == "fixed"
        assert len(bins) == 15

        assert str(bins).startswith("BinTable(")

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

    @pytest.mark.skipif(
        not pandas_avail() or not pyarrow_avail(),
        reason="pandas or pyarrow are not available",
    )
    def test_vectorized_getters(self):

        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)

        assert len(bins.get([1, 1])) == 2
        assert len(bins.get_ids(["chr1", "chr1"], [1, 1])) == 2

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

    @pytest.mark.skipif(not pandas_avail() or not pyarrow_avail(), reason="pandas is not available")
    def test_to_df(self):
        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)

        assert len(bins.to_df()) == len(bins)
        assert len(bins.to_df("chr1")) == 10
        assert len(bins.to_df("chr2:0-200")) == 2
        assert len(bins.to_df("chr2\t0\t200", "BED")) == 2
        with pytest.raises(RuntimeError):
            bins.to_df("chr0")

    def test_iters(self):
        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)

        expected_chroms = []
        expected_starts = []
        expected_ends = []

        for chrom, size in chroms.items():
            num_bins = (size + bins.resolution() - 1) // bins.resolution()
            expected_chroms.extend([chrom] * num_bins)
            starts = list(range(0, size, bins.resolution()))
            ends = [min(pos + bins.resolution(), size) for pos in starts]
            expected_starts.extend(starts)
            expected_ends.extend(ends)

        for chrom, start, end, bin in itertools.zip_longest(expected_chroms, expected_starts, expected_ends, bins):
            assert bin.chrom == chrom
            assert bin.start == start
            assert bin.end == end
