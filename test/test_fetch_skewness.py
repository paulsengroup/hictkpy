# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from math import isnan

import pytest

import hictkpy

testdir = pathlib.Path(__file__).resolve().parent

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (testdir / "data" / "cooler_test_file.mcool", 100_000),
        (testdir / "data" / "hic_test_file.hic", 100_000),
    ],
)


class TestClass:
    def test_fetch_skewness(self, file, resolution):
        f = hictkpy.File(file, resolution)

        def isclose(n1, n2) -> bool:
            import math

            return math.isclose(n1, n2, rel_tol=1.0e-4)

        # venv/bin/python test/scripts/compute_stats_for_testing.py test/data/cooler_test_file.mcool 100000 --metrics skewness --range "" chr2R chr2L --range2 "" chr2R chr2R --normalization NONE VC weight
        assert isclose(f.fetch().skewness(), 59.84043704958881)
        assert isclose(f.fetch().skewness(keep_infs=True), 59.84043704958881)
        assert isclose(f.fetch().skewness(keep_nans=True), 59.84043704958881)
        assert isclose(f.fetch().skewness(keep_nans=True, keep_infs=True), 59.84043704958881)

        assert isclose(f.fetch("chr2R").skewness(), 16.74745540670616)
        assert isclose(f.fetch("chr2R").skewness(keep_infs=True), 16.74745540670616)
        assert isclose(f.fetch("chr2R").skewness(keep_nans=True), 16.74745540670616)
        assert isclose(f.fetch("chr2R").skewness(keep_nans=True, keep_infs=True), 16.74745540670616)

        assert isclose(f.fetch("chr2L", "chr2R").skewness(), 39.00616966835411)
        assert isclose(f.fetch("chr2L", "chr2R").skewness(keep_infs=True), 39.00616966835411)
        assert isclose(f.fetch("chr2L", "chr2R").skewness(keep_nans=True), 39.00616966835411)
        assert isclose(f.fetch("chr2L", "chr2R").skewness(keep_nans=True, keep_infs=True), 39.00616966835411)

        if "VC" in f.avail_normalizations():
            assert isclose(f.fetch(normalization="VC").skewness(), 383.0994807473643)
            assert isclose(f.fetch(normalization="VC").skewness(keep_nans=True), 383.0994807473643)
            assert isnan(f.fetch(normalization="VC").skewness(keep_infs=True))
            assert isnan(f.fetch(normalization="VC").skewness(keep_nans=True, keep_infs=True))

            assert isclose(f.fetch("chr2R", normalization="VC").skewness(), 26.3819375118102)
            assert isclose(f.fetch("chr2R", normalization="VC").skewness(keep_infs=True), 26.3819375118102)
            assert isclose(f.fetch("chr2R", normalization="VC").skewness(keep_nans=True), 26.3819375118102)
            assert isclose(
                f.fetch("chr2R", normalization="VC").skewness(keep_nans=True, keep_infs=True), 26.3819375118102
            )

            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").skewness(), 91.43201554494111)
            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").skewness(keep_infs=True), 91.43201554494111)
            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").skewness(keep_nans=True), 91.43201554494111)
            assert isclose(
                f.fetch("chr2L", "chr2R", normalization="VC").skewness(keep_nans=True, keep_infs=True),
                91.43201554494111,
            )

        norm = "weight" if f.is_cooler() else "ICE"
        assert norm in f.avail_normalizations()

        assert isclose(f.fetch(normalization=norm).skewness(), 52.08028528233407)
        assert isclose(f.fetch(normalization=norm).skewness(keep_infs=True), 52.08028528233407)
        assert isnan(f.fetch(normalization=norm).skewness(keep_nans=True))
        assert isnan(f.fetch(normalization=norm).skewness(keep_nans=True, keep_infs=True))

        assert isclose(f.fetch("chr2R", normalization=norm).skewness(), 10.84820776645202)
        assert isclose(f.fetch("chr2R", normalization=norm).skewness(keep_infs=True), 10.84820776645202)
        assert isnan(f.fetch("chr2R", normalization=norm).skewness(keep_nans=True))
        assert isnan(f.fetch("chr2R", normalization=norm).skewness(keep_nans=True, keep_infs=True))

        assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).skewness(), 16.53465790282013)
        assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).skewness(keep_infs=True), 16.53465790282013)
        assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).skewness(keep_nans=True))
        assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).skewness(keep_nans=True, keep_infs=True))
