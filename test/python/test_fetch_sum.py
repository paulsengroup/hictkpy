# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from math import isclose, isinf, isnan

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


class TestFetchStatsSum:
    def test_fetch_sum(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            # venv/bin/python test/scripts/compute_stats_for_testing.py test/data/cooler_test_file.mcool 100000 --metrics sum --range "" chr2R chr2L --range2 "" chr2R chr2R --normalization NONE VC weight
            assert f.fetch().sum() == 119_208_613
            assert f.fetch().sum(keep_infs=True) == 119_208_613
            assert f.fetch().sum(keep_nans=True) == 119_208_613
            assert f.fetch().sum(keep_nans=True, keep_infs=True) == 119_208_613

            assert f.fetch("chr2R").sum() == 21_513_984
            assert f.fetch("chr2R").sum(keep_infs=True) == 21_513_984
            assert f.fetch("chr2R").sum(keep_nans=True) == 21_513_984
            assert f.fetch("chr2R").sum(keep_nans=True, keep_infs=True) == 21_513_984

            assert f.fetch("chr2L", "chr2R").sum() == 1_483_112
            assert f.fetch("chr2L", "chr2R").sum(keep_infs=True) == 1_483_112
            assert f.fetch("chr2L", "chr2R").sum(keep_nans=True) == 1_483_112
            assert f.fetch("chr2L", "chr2R").sum(keep_nans=True, keep_infs=True) == 1_483_112

            if "VC" in f.avail_normalizations():
                assert isclose(f.fetch(normalization="VC").sum(), 118404896.7857023)
                assert isclose(f.fetch(normalization="VC").sum(keep_nans=True), 118404896.7857023)
                assert isinf(f.fetch(normalization="VC").sum(keep_infs=True))
                assert isinf(f.fetch(normalization="VC").sum(keep_nans=True, keep_infs=True))

                assert isclose(f.fetch("chr2R", normalization="VC").sum(), 20970618.85739681)
                assert isclose(f.fetch("chr2R", normalization="VC").sum(keep_infs=True), 20970618.85739681)
                assert isclose(f.fetch("chr2R", normalization="VC").sum(keep_nans=True), 20970618.85739681)
                assert isclose(
                    f.fetch("chr2R", normalization="VC").sum(keep_nans=True, keep_infs=True), 20970618.85739681
                )

                assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").sum(), 1726683.719399437)
                assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").sum(keep_infs=True), 1726683.719399437)
                assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").sum(keep_nans=True), 1726683.719399437)
                assert isclose(
                    f.fetch("chr2L", "chr2R", normalization="VC").sum(keep_nans=True, keep_infs=True), 1726683.719399437
                )

            norm = "weight" if f.is_cooler() else "ICE"
            assert norm in f.avail_normalizations()

            assert isclose(f.fetch(normalization=norm).sum(), 1830.577410801785)
            assert isclose(f.fetch(normalization=norm).sum(keep_infs=True), 1830.577410801785)
            assert isnan(f.fetch(normalization=norm).sum(keep_nans=True))
            assert isnan(f.fetch(normalization=norm).sum(keep_nans=True, keep_infs=True))

            assert isclose(f.fetch("chr2R", normalization=norm).sum(), 304.2987255987707)
            assert isclose(f.fetch("chr2R", normalization=norm).sum(keep_infs=True), 304.2987255987707)
            assert isnan(f.fetch("chr2R", normalization=norm).sum(keep_nans=True))
            assert isnan(f.fetch("chr2R", normalization=norm).sum(keep_nans=True, keep_infs=True))

            assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).sum(), 21.06706555733309)
            assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).sum(keep_infs=True), 21.06706555733309)
            assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).sum(keep_nans=True))
            assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).sum(keep_nans=True, keep_infs=True))
