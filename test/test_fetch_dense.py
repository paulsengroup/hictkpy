# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

import pytest

import hictkpy

from .helpers import numpy_avail

testdir = pathlib.Path(__file__).resolve().parent

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (testdir / "data" / "cooler_test_file.mcool", 100_000),
        (testdir / "data" / "hic_test_file.hic", 100_000),
    ],
)


def issymmetric(m) -> bool:
    import numpy as np

    assert m.ndim == 2
    if m.size == 0:
        return True
    if m.shape[0] != m.shape[1]:
        return False

    return np.allclose(m, m.T, atol=0.0, rtol=0.0)


@pytest.mark.skipif(not numpy_avail(), reason="numpy is not available")
class TestClass:
    def test_genome_wide(self, file, resolution):
        f = hictkpy.File(file, resolution)
        m = f.fetch().to_numpy()
        assert m.shape == (1380, 1380)
        assert m.sum() == 178_263_235
        assert issymmetric(m)

    def test_cis(self, file, resolution):
        import numpy as np

        f = hictkpy.File(file, resolution)
        m = f.fetch("chr2R:10,000,000-15,000,000").to_numpy()
        assert m.shape == (50, 50)
        assert m.sum() == 6_029_333
        assert issymmetric(m)

        m = f.fetch("chr2R:10,000,000-15,000,000", count_type="int").to_numpy()
        assert m.dtype == np.int32

        m = f.fetch("chr2R:10,000,000-15,000,000", count_type="float").to_numpy()
        assert m.dtype == np.float64

        m = f.fetch("chr2R\t10000000\t15000000", query_type="BED").to_numpy()
        assert m.shape == (50, 50)

        m = f.fetch("chr2L:0-10,000,000", "chr2L:5,000,000-20,000,000").to_numpy()
        assert m.shape == (100, 150)
        assert m.sum() == 6_287_451

        m = f.fetch("chr2L:0-10,000,000", "chr2L:10,000,000-20,000,000").to_numpy()
        assert m.shape == (100, 100)
        assert m.sum() == 761_223

        m = f.fetch("chr2L:0-10,000,000", "chr2L:0-15,000,000").to_numpy()
        assert m.shape == (100, 150)
        assert m.sum() == 12_607_205

    def test_trans(self, file, resolution):
        import numpy as np

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
