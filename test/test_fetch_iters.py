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


def compute_sum(sel):
    return sum(x.count for x in sel)


def compute_nnz(sel):
    return sum(1 for _ in sel)


class TestClass:
    def test_genome_wide(self, file, resolution):
        f = hictkpy.File(file, resolution)

        sel = f.fetch()
        assert compute_sum(sel) == 119_208_613
        assert compute_nnz(sel) == 890_384

    def test_cis(self, file, resolution):
        f = hictkpy.File(file, resolution)

        sel = f.fetch("chr2R:10,000,000-15,000,000")
        assert compute_sum(sel) == 4_519_080

        sel = f.fetch("chr2R:10,000,000-15,000,000", join=True)
        assert compute_sum(sel) == 4_519_080
        sel = f.fetch("chr2R:10,000,000-15,000,000", join=True)
        assert all(x.chrom1 == "chr2R" for x in sel)

        sel = f.fetch("chr2R:10,000,000-15,000,000", count_type="int")
        assert isinstance(compute_sum(sel), int)

        sel = f.fetch("chr2R:10,000,000-15,000,000", count_type="float")
        assert isinstance(compute_sum(sel), float)

        sel = f.fetch("chr2R\t10000000\t15000000", query_type="BED")
        assert compute_nnz(sel) == 1275

    def test_trans(self, file, resolution):
        f = hictkpy.File(file, resolution)

        sel = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000")
        assert compute_sum(sel) == 83_604

        sel = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", join=True)
        assert compute_sum(sel) == 83_604
        assert all(x.chrom1 == "chr2R" and x.chrom2 == "chrX" for x in sel)

        sel = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="int")
        assert isinstance(compute_sum(sel), int)

        sel = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="float")
        assert isinstance(compute_sum(sel), float)

        sel = f.fetch("chr2R\t10000000\t15000000", "chrX\t0\t10000000", query_type="BED")
        assert compute_nnz(sel) == 4995

        with pytest.raises(Exception):
            for p in f.fetch("chrX:0-10,000,000", "chr2R:10,000,000-15,000,000"):
                pass

    def test_balanced(self, file, resolution):
        f = hictkpy.File(file, resolution)

        if f.is_cooler():
            sel = f.fetch("chr2R:10,000,000-15,000,000", normalization="weight")
        else:
            sel = f.fetch("chr2R:10,000,000-15,000,000", normalization="ICE")

        assert np.isclose(59.349524704033215, compute_sum(sel))
