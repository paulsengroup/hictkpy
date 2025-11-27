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


class TestFetchStatsKurtosis:
    @staticmethod
    def _test_fetch_kurtosis(file, resolution, exact):
        def isclose(n1, n2) -> bool:
            import math

            if exact:
                return math.isclose(n1, n2)

            return math.isclose(n1, n2, abs_tol=1.0e-8)

        with hictkpy.File(file, resolution) as f:
            # venv/bin/python test/scripts/compute_stats_for_testing.py test/data/cooler_test_file.mcool 100000 --metrics kurtosis --range "" chr2R chr2L --range2 "" chr2R chr2R --normalization NONE VC weight
            assert isclose(f.fetch().kurtosis(exact=exact), 9102.714436703003)
            assert isclose(f.fetch().kurtosis(exact=exact, keep_infs=True), 9102.714436703003)
            assert isclose(f.fetch().kurtosis(exact=exact, keep_nans=True), 9102.714436703003)
            assert isclose(f.fetch().kurtosis(exact=exact, keep_nans=True, keep_infs=True), 9102.714436703003)

            assert isclose(f.fetch("chr2R").kurtosis(exact=exact), 373.3395734325773)
            assert isclose(f.fetch("chr2R").kurtosis(exact=exact, keep_infs=True), 373.3395734325773)
            assert isclose(f.fetch("chr2R").kurtosis(exact=exact, keep_nans=True), 373.3395734325773)
            assert isclose(f.fetch("chr2R").kurtosis(exact=exact, keep_nans=True, keep_infs=True), 373.3395734325773)

            assert isclose(f.fetch("chr2L", "chr2R").kurtosis(exact=exact), 3688.676195306675)
            assert isclose(f.fetch("chr2L", "chr2R").kurtosis(exact=exact, keep_infs=True), 3688.676195306675)
            assert isclose(f.fetch("chr2L", "chr2R").kurtosis(exact=exact, keep_nans=True), 3688.676195306675)
            assert isclose(
                f.fetch("chr2L", "chr2R").kurtosis(exact=exact, keep_nans=True, keep_infs=True), 3688.676195306675
            )

            if "VC" in f.avail_normalizations():
                assert isclose(f.fetch(normalization="VC").kurtosis(exact=exact), 231607.900103242)
                assert isclose(f.fetch(normalization="VC").kurtosis(exact=exact, keep_nans=True), 231607.900103242)
                assert isnan(f.fetch(normalization="VC").kurtosis(exact=exact, keep_infs=True))
                assert isnan(f.fetch(normalization="VC").kurtosis(exact=exact, keep_nans=True, keep_infs=True))

                assert isclose(f.fetch("chr2R", normalization="VC").kurtosis(exact=exact), 969.1654059364474)
                assert isclose(
                    f.fetch("chr2R", normalization="VC").kurtosis(exact=exact, keep_infs=True), 969.1654059364474
                )
                assert isclose(
                    f.fetch("chr2R", normalization="VC").kurtosis(exact=exact, keep_nans=True), 969.1654059364474
                )
                assert isclose(
                    f.fetch("chr2R", normalization="VC").kurtosis(exact=exact, keep_nans=True, keep_infs=True),
                    969.1654059364474,
                )

                assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").kurtosis(exact=exact), 12761.80816176265)
                assert isclose(
                    f.fetch("chr2L", "chr2R", normalization="VC").kurtosis(exact=exact, keep_infs=True),
                    12761.80816176265,
                )
                assert isclose(
                    f.fetch("chr2L", "chr2R", normalization="VC").kurtosis(exact=exact, keep_nans=True),
                    12761.80816176265,
                )
                assert isclose(
                    f.fetch("chr2L", "chr2R", normalization="VC").kurtosis(exact=exact, keep_nans=True, keep_infs=True),
                    12761.80816176265,
                )

            if f.is_cooler():
                norm = "weight"
                assert isclose(f.fetch(normalization=norm).kurtosis(exact=exact), 9402.779596517719)
                assert isclose(f.fetch(normalization=norm).kurtosis(exact=exact, keep_infs=True), 9402.779596517719)
                assert isnan(f.fetch(normalization=norm).kurtosis(exact=exact, keep_nans=True))
                assert isnan(f.fetch(normalization=norm).kurtosis(exact=exact, keep_nans=True, keep_infs=True))

                assert isclose(f.fetch("chr2R", normalization=norm).kurtosis(exact=exact), 129.2767572383961)
                assert isclose(
                    f.fetch("chr2R", normalization=norm).kurtosis(exact=exact, keep_infs=True), 129.2767572383961
                )
                assert isnan(f.fetch("chr2R", normalization=norm).kurtosis(exact=exact, keep_nans=True))
                assert isnan(f.fetch("chr2R", normalization=norm).kurtosis(exact=exact, keep_nans=True, keep_infs=True))

                assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(exact=exact), 577.4015771468706)
                assert isclose(
                    f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(exact=exact, keep_infs=True),
                    577.4015771468706,
                )
                assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(exact=exact, keep_nans=True))
                assert isnan(
                    f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(exact=exact, keep_nans=True, keep_infs=True)
                )
            else:
                norm = "ICE"
                assert isclose(f.fetch(normalization=norm).kurtosis(exact=exact), 9402.783270262977)
                assert isclose(f.fetch(normalization=norm).kurtosis(exact=exact, keep_infs=True), 9402.783270262977)
                assert isnan(f.fetch(normalization=norm).kurtosis(exact=exact, keep_nans=True))
                assert isnan(f.fetch(normalization=norm).kurtosis(exact=exact, keep_nans=True, keep_infs=True))

                assert isclose(f.fetch("chr2R", normalization=norm).kurtosis(exact=exact), 129.27675506429753)
                assert isclose(
                    f.fetch("chr2R", normalization=norm).kurtosis(exact=exact, keep_infs=True), 129.27675506429753
                )
                assert isnan(f.fetch("chr2R", normalization=norm).kurtosis(exact=exact, keep_nans=True))
                assert isnan(f.fetch("chr2R", normalization=norm).kurtosis(exact=exact, keep_nans=True, keep_infs=True))

                assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(exact=exact), 577.4016006162126)
                assert isclose(
                    f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(exact=exact, keep_infs=True),
                    577.4016006162126,
                )
                assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(exact=exact, keep_nans=True))
                assert isnan(
                    f.fetch("chr2L", "chr2R", normalization=norm).kurtosis(exact=exact, keep_nans=True, keep_infs=True)
                )

    def test_fetch_kurtosis_exact(self, file, resolution):
        TestFetchStatsKurtosis._test_fetch_kurtosis(file, resolution, exact=True)

    def test_fetch_kurtosis_single_pass(self, file, resolution):
        TestFetchStatsKurtosis._test_fetch_kurtosis(file, resolution, exact=False)
