# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from math import isnan

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


class TestFetchStatsSkewness:
    @staticmethod
    def _test_fetch_skewness(file, resolution, exact):
        def isclose(n1, n2) -> bool:
            import math

            return math.isclose(n1, n2, abs_tol=1.0e-5)

        with hictkpy.File(file, resolution) as f:
            # venv/bin/python test/scripts/compute_stats_for_testing.py test/data/cooler_test_file.mcool 100000 --metrics skewness --range "" chr2R chr2L --range2 "" chr2R chr2R --normalization NONE VC weight
            assert isclose(f.fetch().skewness(exact=exact), 59.84043704958881)
            assert isclose(f.fetch().skewness(exact=exact, keep_infs=True), 59.84043704958881)
            assert isclose(f.fetch().skewness(exact=exact, keep_nans=True), 59.84043704958881)
            assert isclose(f.fetch().skewness(exact=exact, keep_nans=True, keep_infs=True), 59.84043704958881)

            assert isclose(f.fetch("chr2R").skewness(exact=exact), 16.74745540670616)
            assert isclose(f.fetch("chr2R").skewness(exact=exact, keep_infs=True), 16.74745540670616)
            assert isclose(f.fetch("chr2R").skewness(exact=exact, keep_nans=True), 16.74745540670616)
            assert isclose(f.fetch("chr2R").skewness(exact=exact, keep_nans=True, keep_infs=True), 16.74745540670616)

            assert isclose(f.fetch("chr2L", "chr2R").skewness(exact=exact), 39.00616966835411)
            assert isclose(f.fetch("chr2L", "chr2R").skewness(exact=exact, keep_infs=True), 39.00616966835411)
            assert isclose(f.fetch("chr2L", "chr2R").skewness(exact=exact, keep_nans=True), 39.00616966835411)
            assert isclose(
                f.fetch("chr2L", "chr2R").skewness(exact=exact, keep_nans=True, keep_infs=True), 39.00616966835411
            )

            if "VC" in f.avail_normalizations():
                assert isclose(f.fetch(normalization="VC").skewness(exact=exact), 383.0994807473643)
                assert isclose(f.fetch(normalization="VC").skewness(exact=exact, keep_nans=True), 383.0994807473643)
                assert isnan(f.fetch(normalization="VC").skewness(exact=exact, keep_infs=True))
                assert isnan(f.fetch(normalization="VC").skewness(exact=exact, keep_nans=True, keep_infs=True))

                assert isclose(f.fetch("chr2R", normalization="VC").skewness(exact=exact), 26.3819375118102)
                assert isclose(
                    f.fetch("chr2R", normalization="VC").skewness(exact=exact, keep_infs=True), 26.3819375118102
                )
                assert isclose(
                    f.fetch("chr2R", normalization="VC").skewness(exact=exact, keep_nans=True), 26.3819375118102
                )
                assert isclose(
                    f.fetch("chr2R", normalization="VC").skewness(exact=exact, keep_nans=True, keep_infs=True),
                    26.3819375118102,
                )

                assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").skewness(exact=exact), 91.43201554494111)
                assert isclose(
                    f.fetch("chr2L", "chr2R", normalization="VC").skewness(exact=exact, keep_infs=True),
                    91.43201554494111,
                )
                assert isclose(
                    f.fetch("chr2L", "chr2R", normalization="VC").skewness(exact=exact, keep_nans=True),
                    91.43201554494111,
                )
                assert isclose(
                    f.fetch("chr2L", "chr2R", normalization="VC").skewness(exact=exact, keep_nans=True, keep_infs=True),
                    91.43201554494111,
                )

            norm = "weight" if f.is_cooler() else "ICE"
            assert norm in f.avail_normalizations()

            assert isclose(f.fetch(normalization=norm).skewness(exact=exact), 52.08028528233407)
            assert isclose(f.fetch(normalization=norm).skewness(exact=exact, keep_infs=True), 52.08028528233407)
            assert isnan(f.fetch(normalization=norm).skewness(exact=exact, keep_nans=True))
            assert isnan(f.fetch(normalization=norm).skewness(exact=exact, keep_nans=True, keep_infs=True))

            assert isclose(f.fetch("chr2R", normalization=norm).skewness(exact=exact), 10.84820776645202)
            assert isclose(
                f.fetch("chr2R", normalization=norm).skewness(exact=exact, keep_infs=True), 10.84820776645202
            )
            assert isnan(f.fetch("chr2R", normalization=norm).skewness(exact=exact, keep_nans=True))
            assert isnan(f.fetch("chr2R", normalization=norm).skewness(exact=exact, keep_nans=True, keep_infs=True))

            assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).skewness(exact=exact), 16.53465790282013)
            assert isclose(
                f.fetch("chr2L", "chr2R", normalization=norm).skewness(exact=exact, keep_infs=True), 16.53465790282013
            )
            assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).skewness(exact=exact, keep_nans=True))
            assert isnan(
                f.fetch("chr2L", "chr2R", normalization=norm).skewness(exact=exact, keep_nans=True, keep_infs=True)
            )

    def test_skewness_exact(self, file, resolution):
        TestFetchStatsSkewness._test_fetch_skewness(file, resolution, exact=True)

    def test_skewness_single_pass(self, file, resolution):
        TestFetchStatsSkewness._test_fetch_skewness(file, resolution, exact=False)
