# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


import pandas as pd

import hictkpy


class TestClass:
    def test_accessors(self):
        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)

        assert len(bins.chromosomes()) == 2
        assert bins.bin_size() == 100
        assert len(bins) == 15

    def test_getters(self):
        chroms = {"chr1": 1000, "chr2": 500}
        bins = hictkpy.BinTable(chroms, 100)

        assert bins[1].chrom == "chr1"
        assert bins[1].start == 100
        assert bins[1].end == 200

        assert bins.get(1).chrom == "chr1"
        assert bins.get(1).start == 100
        assert bins.get(1).end == 200
        assert bins.get(9999) is None

        assert bins.get("chr1", 153).id == 1
        assert bins.get("abc", 100) is None

    def test_merge(self):
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
