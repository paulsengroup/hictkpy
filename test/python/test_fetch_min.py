# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from math import isclose, isnan

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


class TestFetchStatsMin:
    def test_fetch_min(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            # venv/bin/python test/scripts/compute_stats_for_testing.py test/data/cooler_test_file.mcool 100000 --metrics min --range "" chr2R chr2L --range2 "" chr2R chr2R --normalization NONE VC weight
            assert f.fetch().min() == 1
            assert f.fetch().min(keep_infs=True) == 1
            assert f.fetch().min(keep_nans=True) == 1
            assert f.fetch().min(keep_nans=True, keep_infs=True) == 1

            assert f.fetch("chr2R").min() == 1
            assert f.fetch("chr2R").min(keep_infs=True) == 1
            assert f.fetch("chr2R").min(keep_nans=True) == 1
            assert f.fetch("chr2R").min(keep_nans=True, keep_infs=True) == 1

            assert f.fetch("chr2L", "chr2R").min() == 1
            assert f.fetch("chr2L", "chr2R").min(keep_infs=True) == 1
            assert f.fetch("chr2L", "chr2R").min(keep_nans=True) == 1
            assert f.fetch("chr2L", "chr2R").min(keep_nans=True, keep_infs=True) == 1

            if "VC" in f.avail_normalizations():
                assert isclose(f.fetch(normalization="VC").min(), 0.006248577402113111)
                assert isclose(f.fetch(normalization="VC").min(keep_infs=True), 0.006248577402113111)
                assert isclose(f.fetch(normalization="VC").min(keep_nans=True), 0.006248577402113111)
                assert isclose(f.fetch(normalization="VC").min(keep_nans=True, keep_infs=True), 0.006248577402113111)

                assert isclose(f.fetch("chr2R", normalization="VC").min(), 3.024889833701339)
                assert isclose(f.fetch("chr2R", normalization="VC").min(keep_infs=True), 3.024889833701339)
                assert isclose(f.fetch("chr2R", normalization="VC").min(keep_nans=True), 3.024889833701339)
                assert isclose(
                    f.fetch("chr2R", normalization="VC").min(keep_nans=True, keep_infs=True), 3.024889833701339
                )

                assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").min(), 1.179576020123557)
                assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").min(keep_infs=True), 1.179576020123557)
                assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").min(keep_nans=True), 1.179576020123557)
                assert isclose(
                    f.fetch("chr2L", "chr2R", normalization="VC").min(keep_nans=True, keep_infs=True), 1.179576020123557
                )

            # Even though the test .mcool and .hic files are supposed to be identical, normalized interactions are expected to be very close, but not identical.
            # This is due to differences in how balancing weights are stored Cooler and .hic files:
            # - Cooler uses float64 multiplicative weights
            # - .hic v8 and older use float64 divisive weights
            # - .hic v9 uses float32 divisive weights
            # Usually the differences are small enough to not being picked up by isclose.
            # However, this is not the case for very small or very large values (like those returned by min() and max())
            def isclose_(n1, n2) -> bool:
                import math

                return math.isclose(n1, n2, rel_tol=1.0e-6)

            norm = "weight" if f.is_cooler() else "ICE"
            assert norm in f.avail_normalizations()

            assert isclose_(f.fetch(normalization=norm).min(), 1.514271829598328e-05)
            assert isclose_(f.fetch(normalization=norm).min(keep_infs=True), 1.514271829598328e-05)
            assert isnan(f.fetch(normalization=norm).min(keep_nans=True))
            assert isnan(f.fetch(normalization=norm).min(keep_nans=True, keep_infs=True))

            assert isclose_(f.fetch("chr2R", normalization=norm).min(), 3.324204473054833e-05)
            assert isclose_(f.fetch("chr2R", normalization=norm).min(keep_infs=True), 3.324204473054833e-05)
            assert isnan(f.fetch("chr2R", normalization=norm).min(keep_nans=True))
            assert isnan(f.fetch("chr2R", normalization=norm).min(keep_nans=True, keep_infs=True))

            assert isclose_(f.fetch("chr2L", "chr2R", normalization=norm).min(), 2.117120750101533e-05)
            assert isclose_(f.fetch("chr2L", "chr2R", normalization=norm).min(keep_infs=True), 2.117120750101533e-05)
            assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).min(keep_nans=True))
            assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).min(keep_nans=True, keep_infs=True))
