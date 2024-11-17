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
    def test_fetch_kurtosis(self, file, resolution):
        f = hictkpy.File(file, resolution)

        def isclose(n1, n2) -> bool:
            import math

            return math.isclose(n1, n2, rel_tol=1.0e-4)

        # venv/bin/python test/scripts/compute_stats_for_testing.py test/data/cooler_test_file.mcool 100000 --metrics kurtosis --range "" chr2R chr2L --range2 "" chr2R chr2R --normalization NONE VC weight
        assert isclose(f.fetch().kurtosis(), 9102.714436703003)
        assert isclose(f.fetch().kurtosis(keep_infs=True), 9102.714436703003)
        assert isclose(f.fetch().kurtosis(keep_nans=True), 9102.714436703003)
        assert isclose(f.fetch().kurtosis(keep_nans=True, keep_infs=True), 9102.714436703003)

        assert isclose(f.fetch("chr2R").kurtosis(), 373.3395734325773)
        assert isclose(f.fetch("chr2R").kurtosis(keep_infs=True), 373.3395734325773)
        assert isclose(f.fetch("chr2R").kurtosis(keep_nans=True), 373.3395734325773)
        assert isclose(f.fetch("chr2R").kurtosis(keep_nans=True, keep_infs=True), 373.3395734325773)

        assert isclose(f.fetch("chr2L", "chr2R").kurtosis(), 3688.676195306675)
        assert isclose(f.fetch("chr2L", "chr2R").kurtosis(keep_infs=True), 3688.676195306675)
        assert isclose(f.fetch("chr2L", "chr2R").kurtosis(keep_nans=True), 3688.676195306675)
        assert isclose(f.fetch("chr2L", "chr2R").kurtosis(keep_nans=True, keep_infs=True), 3688.676195306675)

        if "VC" in f.avail_normalizations():
            assert isclose(f.fetch(normalization="VC").kurtosis(), 231607.900103242)
            assert isclose(f.fetch(normalization="VC").kurtosis(keep_nans=True), 231607.900103242)
            assert isnan(f.fetch(normalization="VC").kurtosis(keep_infs=True))
            assert isnan(f.fetch(normalization="VC").kurtosis(keep_nans=True, keep_infs=True))

            assert isclose(f.fetch("chr2R", normalization="VC").kurtosis(), 969.1654059364474)
            assert isclose(f.fetch("chr2R", normalization="VC").kurtosis(keep_infs=True), 969.1654059364474)
            assert isclose(f.fetch("chr2R", normalization="VC").kurtosis(keep_nans=True), 969.1654059364474)
            assert isclose(
                f.fetch("chr2R", normalization="VC").kurtosis(keep_nans=True, keep_infs=True), 969.1654059364474
            )

            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").kurtosis(), 12761.80816176265)
            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").kurtosis(keep_infs=True), 12761.80816176265)
            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").kurtosis(keep_nans=True), 12761.80816176265)
            assert isclose(
                f.fetch("chr2L", "chr2R", normalization="VC").kurtosis(keep_nans=True, keep_infs=True),
                12761.80816176265,
            )

        norm = "weight" if f.is_cooler() else "ICE"
        assert norm in f.avail_normalizations()

        assert isclose(f.fetch(normalization=norm).kurtosis(), 9402.779596517719)
        assert isclose(f.fetch(normalization=norm).kurtosis(keep_infs=True), 9402.779596517719)
        assert isnan(f.fetch(normalization=norm).kurtosis(keep_nans=True))
        assert isnan(f.fetch(normalization=norm).kurtosis(keep_nans=True, keep_infs=True))

        assert isclose(f.fetch("chr2R", normalization=norm).kurtosis(), 129.2767572383961)
        assert isclose(f.fetch("chr2R", normalization=norm).kurtosis(keep_infs=True), 129.2767572383961)
        assert isnan(f.fetch("chr2R", normalization=norm).kurtosis(keep_nans=True))
        assert isnan(f.fetch("chr2R", normalization=norm).kurtosis(keep_nans=True, keep_infs=True))

        assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(), 577.4015771468706)
        assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(keep_infs=True), 577.4015771468706)
        assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(keep_nans=True))
        assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(keep_nans=True, keep_infs=True))
