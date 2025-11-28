# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

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


class TestFetchStatsNNZ:
    def test_fetch_nnz(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            # venv/bin/python test/scripts/compute_stats_for_testing.py test/data/cooler_test_file.mcool 100000 --metrics nnz --range "" chr2R chr2L --range2 "" chr2R chr2R --normalization NONE VC weight
            assert f.fetch().nnz() == 890_384
            assert f.fetch().nnz(keep_infs=True) == 890_384
            assert f.fetch().nnz(keep_nans=True) == 890_384
            assert f.fetch().nnz(keep_nans=True, keep_infs=True) == 890_384

            assert f.fetch("chr2R").nnz() == 31_900
            assert f.fetch("chr2R").nnz(keep_infs=True) == 31_900
            assert f.fetch("chr2R").nnz(keep_nans=True) == 31_900
            assert f.fetch("chr2R").nnz(keep_nans=True, keep_infs=True) == 31_900

            assert f.fetch("chr2L", "chr2R").nnz() == 59_124
            assert f.fetch("chr2L", "chr2R").nnz(keep_infs=True) == 59_124
            assert f.fetch("chr2L", "chr2R").nnz(keep_nans=True) == 59_124
            assert f.fetch("chr2L", "chr2R").nnz(keep_nans=True, keep_infs=True) == 59_124

            if "VC" in f.avail_normalizations():
                assert f.fetch(normalization="VC").nnz() == 887_644
                assert f.fetch(normalization="VC").nnz(keep_infs=True) == 890_384
                assert f.fetch(normalization="VC").nnz(keep_nans=True) == 887_644
                assert f.fetch(normalization="VC").nnz(keep_nans=True, keep_infs=True) == 890_384

                assert f.fetch("chr2R", normalization="VC").nnz() == 31_900
                assert f.fetch("chr2R", normalization="VC").nnz(keep_infs=True) == 31_900
                assert f.fetch("chr2R", normalization="VC").nnz(keep_nans=True) == 31_900
                assert f.fetch("chr2R", normalization="VC").nnz(keep_nans=True, keep_infs=True) == 31_900

                assert f.fetch("chr2L", "chr2R", normalization="VC").nnz() == 59_124
                assert f.fetch("chr2L", "chr2R", normalization="VC").nnz(keep_infs=True) == 59_124
                assert f.fetch("chr2L", "chr2R", normalization="VC").nnz(keep_nans=True) == 59_124
                assert f.fetch("chr2L", "chr2R", normalization="VC").nnz(keep_nans=True, keep_infs=True) == 59_124

            norm = "weight" if f.is_cooler() else "ICE"
            assert norm in f.avail_normalizations()

            assert f.fetch(normalization=norm).nnz() == 834_993
            assert f.fetch(normalization=norm).nnz(keep_infs=True) == 834_993
            assert f.fetch(normalization=norm).nnz(keep_nans=True) == 890_384
            assert f.fetch(normalization=norm).nnz(keep_nans=True, keep_infs=True) == 890_384

            assert f.fetch("chr2R", normalization=norm).nnz() == 28_919
            assert f.fetch("chr2R", normalization=norm).nnz(keep_infs=True) == 28_919
            assert f.fetch("chr2R", normalization=norm).nnz(keep_nans=True) == 31_900
            assert f.fetch("chr2R", normalization=norm).nnz(keep_nans=True, keep_infs=True) == 31_900

            assert f.fetch("chr2L", "chr2R", normalization=norm).nnz() == 55_432
            assert f.fetch("chr2L", "chr2R", normalization=norm).nnz(keep_infs=True) == 55_432
            assert f.fetch("chr2L", "chr2R", normalization=norm).nnz(keep_nans=True) == 59_124
            assert f.fetch("chr2L", "chr2R", normalization=norm).nnz(keep_nans=True, keep_infs=True) == 59_124
