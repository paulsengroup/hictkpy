# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pytest

import hictkpy

from .helpers import get_test_dir, numpy_avail

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (get_test_dir() / "data" / "cooler_test_file.mcool", 100_000),
        (get_test_dir() / "data" / "hic_test_file.hic", 100_000),
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
class TestFetchDense:
    def test_genome_wide(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            m = f.fetch().to_numpy()
        assert m.shape == (1380, 1380)
        assert m.sum() == 178_263_235
        assert issymmetric(m)

    def test_cis(self, file, resolution):
        import numpy as np

        with hictkpy.File(file, resolution) as f:
            m = f.fetch("chr2R:10,000,000-15,000,000").to_numpy()
            assert m.shape == (50, 50)
            assert m.sum() == 6_029_333
            assert issymmetric(m)

            m = f.fetch("chr2R:10,000,000-15,000,000", count_type="int32").to_numpy()
            assert m.dtype == np.int32

            m = f.fetch("chr2R:10,000,000-15,000,000", count_type=float).to_numpy()
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

        with hictkpy.File(file, resolution) as f:
            m = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000").to_numpy()
            assert m.shape == (50, 100)
            assert m.sum() == 83_604

            m = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="int32").to_numpy()
            assert m.dtype == np.int32

            m = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type=float).to_numpy()
            assert m.dtype == np.float64

            m = f.fetch("chr2R\t10000000\t15000000", "chrX\t0\t10000000", query_type="BED").to_numpy()
            assert m.shape == (50, 100)

    def test_diagonal_band(self, file, resolution):
        import numpy as np

        with hictkpy.File(file, resolution) as f:
            m = f.fetch(diagonal_band_width=100).to_numpy()

        assert m.sum() == 149208057
        assert m.shape[0] == m.shape[1]
        assert m.shape[0] == 1380

        m_triu = np.triu(m)
        m_tril = np.tril(m)

        for i in range(m.shape[0]):
            row1 = m_triu[i, :]
            row2 = m_tril[i, :]
            assert (row1 != 0).sum() <= 100
            assert (row2 != 0).sum() <= 100
