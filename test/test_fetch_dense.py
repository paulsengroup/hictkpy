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
        m = f.fetch().to_numpy()
        assert m.shape == (1380, 1380)
        assert m.sum() == 178_263_235

    def test_cis(self, file, resolution):
        f = hictkpy.File(file, resolution)
        m = f.fetch("chr2R:10,000,000-15,000,000").to_numpy()
        assert m.shape == (50, 50)
        assert m.sum() == 6_029_333

        m = f.fetch("chr2R:10,000,000-15,000,000", count_type="int").to_numpy()
        assert m.dtype == np.int32

        m = f.fetch("chr2R:10,000,000-15,000,000", count_type="float").to_numpy()
        assert m.dtype == np.float64

        m = f.fetch("chr2R\t10000000\t15000000", query_type="BED").to_numpy()
        assert m.shape == (50, 50)

    def test_trans(self, file, resolution):
        f = hictkpy.File(file, resolution)
        m = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000").to_numpy()
        assert m.shape == (50, 100)
        assert m.sum() == 83_604

        m = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="int").to_numpy()
        assert m.dtype == np.int32

        m = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="float").to_numpy()
        assert m.dtype == np.float64

        m = f.fetch("chr2R\t10000000\t15000000", "chrX\t0\t10000000", query_type="BED").to_numpy()
        assert m.shape == (50, 100)
