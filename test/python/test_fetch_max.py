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


class TestFetchStatsMax:
    def test_fetch_max(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            assert f.fetch().max() == 660_210
            assert f.fetch().max(keep_infs=True) == 660_210
            assert f.fetch().max(keep_nans=True) == 660_210
            assert f.fetch().max(keep_nans=True, keep_infs=True) == 660_210

            assert f.fetch("chr2R").max() == 229_239
            assert f.fetch("chr2R").max(keep_infs=True) == 229_239
            assert f.fetch("chr2R").max(keep_nans=True) == 229_239
            assert f.fetch("chr2R").max(keep_nans=True, keep_infs=True) == 229_239

            assert f.fetch("chr2L", "chr2R").max() == 3_051
            assert f.fetch("chr2L", "chr2R").max(keep_infs=True) == 3_051
            assert f.fetch("chr2L", "chr2R").max(keep_nans=True) == 3_051
            assert f.fetch("chr2L", "chr2R").max(keep_nans=True, keep_infs=True) == 3_051

            if "VC" in f.avail_normalizations():
                assert isclose(f.fetch(normalization="VC").max(), 2929440.776073241)
                assert isclose(f.fetch(normalization="VC").max(keep_nans=True), 2929440.776073241)
                assert isinf(f.fetch(normalization="VC").max(keep_infs=True))
                assert isinf(f.fetch(normalization="VC").max(keep_nans=True, keep_infs=True))

                assert isclose(f.fetch("chr2R", normalization="VC").max(), 322150.9086956783)
                assert isclose(f.fetch("chr2R", normalization="VC").max(keep_infs=True), 322150.9086956783)
                assert isclose(f.fetch("chr2R", normalization="VC").max(keep_nans=True), 322150.9086956783)
                assert isclose(
                    f.fetch("chr2R", normalization="VC").max(keep_nans=True, keep_infs=True), 322150.9086956783
                )

                assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").max(), 20736.04261680762)
                assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").max(keep_infs=True), 20736.04261680762)
                assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").max(keep_nans=True), 20736.04261680762)
                assert isclose(
                    f.fetch("chr2L", "chr2R", normalization="VC").max(keep_nans=True, keep_infs=True), 20736.04261680762
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

            assert isclose_(f.fetch(normalization=norm).max(), 9.307500711565524)
            assert isclose_(f.fetch(normalization=norm).max(keep_infs=True), 9.307500711565524)
            assert isnan(f.fetch(normalization=norm).max(keep_nans=True))
            assert isnan(f.fetch(normalization=norm).max(keep_nans=True, keep_infs=True))

            assert isclose_(f.fetch("chr2R", normalization=norm).max(), 1.386381144237653)
            assert isclose_(f.fetch("chr2R", normalization=norm).max(keep_infs=True), 1.386381144237653)
            assert isnan(f.fetch("chr2R", normalization=norm).max(keep_nans=True))
            assert isnan(f.fetch("chr2R", normalization=norm).max(keep_nans=True, keep_infs=True))

            assert isclose_(f.fetch("chr2L", "chr2R", normalization=norm).max(), 0.02017488661930269)
            assert isclose_(f.fetch("chr2L", "chr2R", normalization=norm).max(keep_infs=True), 0.02017488661930269)
            assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).max(keep_nans=True))
            assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).max(keep_nans=True, keep_infs=True))
