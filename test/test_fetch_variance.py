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
    def test_fetch_variance(self, file, resolution):
        f = hictkpy.File(file, resolution)

        def isclose(n1, n2) -> bool:
            import math

            return math.isclose(n1, n2, rel_tol=1.0e-4)

        # venv/bin/python test/scripts/compute_stats_for_testing.py test/data/cooler_test_file.mcool 100000 --metrics variance --range "" chr2R chr2L --range2 "" chr2R chr2R --normalization NONE VC weight
        assert isclose(f.fetch().variance(), 5373287.890326086)
        assert isclose(f.fetch().variance(keep_infs=True), 5373287.890326086)
        assert isclose(f.fetch().variance(keep_nans=True), 5373287.890326086)
        assert isclose(f.fetch().variance(keep_nans=True, keep_infs=True), 5373287.890326086)

        assert isclose(f.fetch("chr2R").variance(), 30250660.97409304)
        assert isclose(f.fetch("chr2R").variance(keep_infs=True), 30250660.97409304)
        assert isclose(f.fetch("chr2R").variance(keep_nans=True), 30250660.97409304)
        assert isclose(f.fetch("chr2R").variance(keep_nans=True, keep_infs=True), 30250660.97409304)

        assert isclose(f.fetch("chr2L", "chr2R").variance(), 775.8613589939481)
        assert isclose(f.fetch("chr2L", "chr2R").variance(keep_infs=True), 775.8613589939481)
        assert isclose(f.fetch("chr2L", "chr2R").variance(keep_nans=True), 775.8613589939481)
        assert isclose(f.fetch("chr2L", "chr2R").variance(keep_nans=True, keep_infs=True), 775.8613589939481)

        if "VC" in f.avail_normalizations():
            assert isclose(f.fetch(normalization="VC").variance(), 19247119.46018428)
            assert isclose(f.fetch(normalization="VC").variance(keep_nans=True), 19247119.46018428)
            assert isnan(f.fetch(normalization="VC").variance(keep_infs=True))
            assert isnan(f.fetch(normalization="VC").variance(keep_nans=True, keep_infs=True))

            assert isclose(f.fetch("chr2R", normalization="VC").variance(), 37984857.23777649)
            assert isclose(f.fetch("chr2R", normalization="VC").variance(keep_infs=True), 37984857.23777649)
            assert isclose(f.fetch("chr2R", normalization="VC").variance(keep_nans=True), 37984857.23777649)
            assert isclose(
                f.fetch("chr2R", normalization="VC").variance(keep_nans=True, keep_infs=True), 37984857.23777649
            )

            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").variance(), 16780.3290659542)
            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").variance(keep_infs=True), 16780.3290659542)
            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").variance(keep_nans=True), 16780.3290659542)
            assert isclose(
                f.fetch("chr2L", "chr2R", normalization="VC").variance(keep_nans=True, keep_infs=True), 16780.3290659542
            )

        norm = "weight" if f.is_cooler() else "ICE"
        assert norm in f.avail_normalizations()

        assert isclose(f.fetch(normalization=norm).variance(), 0.001012475106404277)
        assert isclose(f.fetch(normalization=norm).variance(keep_infs=True), 0.001012475106404277)
        assert isnan(f.fetch(normalization=norm).variance(keep_nans=True))
        assert isnan(f.fetch(normalization=norm).variance(keep_nans=True, keep_infs=True))

        assert isclose(f.fetch("chr2R", normalization=norm).variance(), 0.005035086513187063)
        assert isclose(f.fetch("chr2R", normalization=norm).variance(keep_infs=True), 0.005035086513187063)
        assert isnan(f.fetch("chr2R", normalization=norm).variance(keep_nans=True))
        assert isnan(f.fetch("chr2R", normalization=norm).variance(keep_nans=True, keep_infs=True))

        assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).variance(), 1.026472030208059e-07)
        assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).variance(keep_infs=True), 1.026472030208059e-07)
        assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).variance(keep_nans=True))
        assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).variance(keep_nans=True, keep_infs=True))
