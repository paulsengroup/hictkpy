# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import math

import pytest

import hictkpy

from .helpers import get_test_dir

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (get_test_dir() / "data" / "cooler_test_file.mcool", 100_000),
        (get_test_dir() / "data" / "hic_test_file.hic", 100_000),
    ],
)


def compute_sum(sel):
    return sum(x.count for x in sel)


def compute_nnz(sel):
    return sum(1 for _ in sel)


class TestFetchIterator:
    def test_genome_wide(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            sel = f.fetch()
            assert compute_sum(sel) == 119_208_613
            assert compute_nnz(sel) == 890_384

    def test_cis(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            sel = f.fetch("chr2R:10,000,000-15,000,000")
            assert compute_sum(sel) == 4_519_080

            sel = f.fetch("chr2R:10,000,000-15,000,000", join=True)
            assert compute_sum(sel) == 4_519_080
            sel = f.fetch("chr2R:10,000,000-15,000,000", join=True)
            assert all(x.chrom1 == "chr2R" for x in sel)

            sel = f.fetch("chr2R:10,000,000-15,000,000", count_type="int32")
            assert isinstance(compute_sum(sel), int)

            sel = f.fetch("chr2R:10,000,000-15,000,000", count_type=float)
            assert isinstance(compute_sum(sel), float)

            sel = f.fetch("chr2R\t10000000\t15000000", query_type="BED")
            assert compute_nnz(sel) == 1275

    def test_trans(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            sel = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000")
            assert compute_sum(sel) == 83_604

            sel = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", join=True)
            assert compute_sum(sel) == 83_604
            assert all(x.chrom1 == "chr2R" and x.chrom2 == "chrX" for x in sel)

            sel = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type=int)
            assert isinstance(compute_sum(sel), int)

            sel = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type=float)
            assert isinstance(compute_sum(sel), float)

            sel = f.fetch("chr2R\t10000000\t15000000", "chrX\t0\t10000000", query_type="BED")
            assert compute_nnz(sel) == 4995

            with pytest.raises(Exception):
                for _ in f.fetch("chrX:0-10,000,000", "chr2R:10,000,000-15,000,000"):
                    pass

    def test_balanced(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            if f.is_cooler():
                sel = f.fetch("chr2R:10,000,000-15,000,000", normalization="weight")
            else:
                sel = f.fetch("chr2R:10,000,000-15,000,000", normalization="ICE")

            assert math.isclose(59.349524704033215, compute_sum(sel), rel_tol=1.0e-5, abs_tol=1.0e-8)
