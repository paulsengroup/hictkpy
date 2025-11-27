# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import math

import pytest

import hictkpy

from .helpers import get_test_dir, numpy_avail, scipy_avail

pytestmark = pytest.mark.parametrize(
    "file,resolution,low_memory",
    [
        (get_test_dir() / "data" / "cooler_test_file.mcool", 100_000, False),
        (get_test_dir() / "data" / "hic_test_file.hic", 100_000, False),
        (get_test_dir() / "data" / "cooler_test_file.mcool", 100_000, True),
        (get_test_dir() / "data" / "hic_test_file.hic", 100_000, True),
    ],
)


@pytest.mark.skipif(not scipy_avail(), reason="scipy is not available")
class TestFetchSparse:
    def test_genome_wide(self, file, resolution, low_memory):
        with hictkpy.File(file, resolution) as f:
            m = f.fetch().to_coo(low_memory=low_memory)

        assert m.shape == (1380, 1380)
        assert m.sum() == 119_208_613

    def test_cis(self, file, resolution, low_memory):
        with hictkpy.File(file, resolution) as f:
            m = f.fetch("chr2R:10,000,000-15,000,000").to_coo(low_memory=low_memory)

            assert m.shape == (50, 50)
            assert m.sum() == 4_519_080

            if numpy_avail():
                import numpy as np

                m = f.fetch("chr2R:10,000,000-15,000,000", count_type="int32").to_coo(low_memory=low_memory)
                assert m.dtype == np.int32

                m = f.fetch("chr2R:10,000,000-15,000,000", count_type=float).to_coo(low_memory=low_memory)
                assert m.dtype == np.float64

            m = f.fetch("chr2R\t10000000\t15000000", query_type="BED").to_coo(low_memory=low_memory)
            assert m.shape == (50, 50)

            m = f.fetch("chr2L:0-10,000,000", "chr2L:10,000,000-20,000,000").to_coo(low_memory=low_memory)
            assert m.shape == (100, 100)
            assert m.sum() == 761_223

            m = f.fetch("chr2L:0-10,000,000", "chr2L:0-15,000,000").to_coo(low_memory=low_memory)
            assert m.shape == (100, 150)
            assert m.sum() == 9_270_385

    def test_trans(self, file, resolution, low_memory):
        with hictkpy.File(file, resolution) as f:
            m = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000").to_coo(low_memory=low_memory)
            assert m.shape == (50, 100)
            assert m.sum() == 83_604

            if numpy_avail():
                import numpy as np

                m = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="int32").to_coo(
                    low_memory=low_memory
                )
                assert m.dtype == np.int32

                m = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type=float).to_coo(
                    low_memory=low_memory
                )
                assert m.dtype == np.float64

            m = f.fetch("chr2R\t10000000\t15000000", "chrX\t0\t10000000", query_type="BED").to_coo(
                low_memory=low_memory
            )
            assert m.shape == (50, 100)

    def test_balanced(self, file, resolution, low_memory):
        with hictkpy.File(file, resolution) as f:
            if f.is_cooler():
                m = f.fetch("chr2R:10,000,000-15,000,000", normalization="weight").to_coo(low_memory=low_memory)
            else:
                m = f.fetch("chr2R:10,000,000-15,000,000", normalization="ICE").to_coo(low_memory=low_memory)

        assert math.isclose(59.349524704033215, m.sum(), rel_tol=1.0e-5, abs_tol=1.0e-8)

    def test_diagonal_band(self, file, resolution, low_memory):
        import scipy.sparse as ss

        with hictkpy.File(file, resolution) as f:
            m = f.fetch(diagonal_band_width=100).to_csr(low_memory=low_memory)

        assert m.sum() == 104681024
        assert m.shape[0] == m.shape[1]
        assert m.shape[0] == 1380

        m_triu = ss.triu(m, format="csr")
        m_tril = ss.tril(m, format="csr")

        for i in range(m.shape[0]):
            row1 = m_triu[i, :].toarray()
            row2 = m_tril[i, :].toarray()
            assert (row1 != 0).sum() <= 100
            assert (row2 != 0).sum() <= 100
