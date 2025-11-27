# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import math

import pytest

import hictkpy

from .helpers import get_test_dir

file = get_test_dir() / "data" / "cooler_test_file.mcool"
resolution = 100_000


class TestPixel:
    def test_coo_int_pixel(self):
        p = hictkpy.Pixel(0, 1, 123)

        assert p.bin1_id == 0
        assert p.bin2_id == 1
        assert isinstance(p.count, int)
        assert p.count == 123

        invalid_attrs = ("bin1", "bin2", "chrom1", "start1", "end1", "chrom2", "start2", "end2")

        for attr in invalid_attrs:
            with pytest.raises(AttributeError, match="does not have.*genomic coordinates"):
                _ = getattr(p, attr)

    def test_coo_float_pixel(self):
        p = hictkpy.Pixel(0, 1, 1.23)

        assert p.bin1_id == 0
        assert p.bin2_id == 1
        assert isinstance(p.count, float)
        assert math.isclose(p.count, 1.23)

        invalid_attrs = ("bin1", "bin2", "chrom1", "start1", "end1", "chrom2", "start2", "end2")

        for attr in invalid_attrs:
            with pytest.raises(AttributeError, match="does not have.*genomic coordinates"):
                _ = getattr(p, attr)

    def test_bg2_int_pixel(self):
        with hictkpy.File(file, resolution) as f:
            bins = f.bins()

        bin1 = bins.get(0)
        bin2 = bins.get(1)

        p = hictkpy.Pixel(bin1, bin2, 123)

        assert p.bin1_id == 0
        assert p.bin1.id == 0
        assert p.bin1.chrom == "chr2L"
        assert p.bin1.start == 0
        assert p.bin1.end == 100_000

        assert p.bin2_id == 1
        assert p.bin2.id == 1
        assert p.bin2.chrom == "chr2L"
        assert p.bin2.start == 100_000
        assert p.bin2.end == 200_000

        assert isinstance(p.count, int)
        assert p.count == 123

    def test_bg2_float_pixel(self):
        with hictkpy.File(file, resolution) as f:
            bins = f.bins()

        bin1 = bins.get(0)
        bin2 = bins.get(1)

        p = hictkpy.Pixel(bin1, bin2, 1.23)

        assert p.bin1_id == 0
        assert p.bin1.id == 0
        assert p.bin1.chrom == "chr2L"
        assert p.bin1.start == 0
        assert p.bin1.end == 100_000

        assert p.bin2_id == 1
        assert p.bin2.id == 1
        assert p.bin2.chrom == "chr2L"
        assert p.bin2.start == 100_000
        assert p.bin2.end == 200_000

        assert isinstance(p.count, float)
        assert p.count == 1.23
