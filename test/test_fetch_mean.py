# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from math import isclose, isinf, isnan

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
    def test_fetch_mean(self, file, resolution):
        f = hictkpy.File(file, resolution)

        # venv/bin/python test/scripts/compute_stats_for_testing.py test/data/cooler_test_file.mcool 100000 --metrics mean --range "" chr2R chr2L --range2 "" chr2R chr2R --normalization NONE VC weight
        assert isclose(f.fetch().mean(), 133.8844959028913)
        assert isclose(f.fetch().mean(keep_infs=True), 133.8844959028913)
        assert isclose(f.fetch().mean(keep_nans=True), 133.8844959028913)
        assert isclose(f.fetch().mean(keep_nans=True, keep_infs=True), 133.8844959028913)

        assert isclose(f.fetch("chr2R").mean(), 674.4195611285267)
        assert isclose(f.fetch("chr2R").mean(keep_infs=True), 674.4195611285267)
        assert isclose(f.fetch("chr2R").mean(keep_nans=True), 674.4195611285267)
        assert isclose(f.fetch("chr2R").mean(keep_nans=True, keep_infs=True), 674.4195611285267)

        assert isclose(f.fetch("chr2L", "chr2R").mean(), 25.08477098978418)
        assert isclose(f.fetch("chr2L", "chr2R").mean(keep_infs=True), 25.08477098978418)
        assert isclose(f.fetch("chr2L", "chr2R").mean(keep_nans=True), 25.08477098978418)
        assert isclose(f.fetch("chr2L", "chr2R").mean(keep_nans=True, keep_infs=True), 25.08477098978418)

        if "VC" in f.avail_normalizations():
            assert isclose(f.fetch(normalization="VC").mean(), 133.3923248348463)
            assert isclose(f.fetch(normalization="VC").mean(keep_nans=True), 133.3923248348463)
            assert isinf(f.fetch(normalization="VC").mean(keep_infs=True))
            assert isinf(f.fetch(normalization="VC").mean(keep_nans=True, keep_infs=True))

            assert isclose(f.fetch("chr2R", normalization="VC").mean(), 657.38617107827)
            assert isclose(f.fetch("chr2R", normalization="VC").mean(keep_infs=True), 657.38617107827)
            assert isclose(f.fetch("chr2R", normalization="VC").mean(keep_nans=True), 657.38617107827)
            assert isclose(f.fetch("chr2R", normalization="VC").mean(keep_nans=True, keep_infs=True), 657.38617107827)

            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").mean(), 29.20444691494886)
            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").mean(keep_infs=True), 29.20444691494886)
            assert isclose(f.fetch("chr2L", "chr2R", normalization="VC").mean(keep_nans=True), 29.20444691494886)
            assert isclose(
                f.fetch("chr2L", "chr2R", normalization="VC").mean(keep_nans=True, keep_infs=True), 29.20444691494886
            )

        norm = "weight" if f.is_cooler() else "ICE"
        assert norm in f.avail_normalizations()

        assert isclose(f.fetch(normalization=norm).mean(), 0.002192326655195654)
        assert isclose(f.fetch(normalization=norm).mean(keep_infs=True), 0.002192326655195654)
        assert isnan(f.fetch(normalization=norm).mean(keep_nans=True))
        assert isnan(f.fetch(normalization=norm).mean(keep_nans=True, keep_infs=True))

        assert isclose(f.fetch("chr2R", normalization=norm).mean(), 0.0105224497942104)
        assert isclose(f.fetch("chr2R", normalization=norm).mean(keep_infs=True), 0.0105224497942104)
        assert isnan(f.fetch("chr2R", normalization=norm).mean(keep_nans=True))
        assert isnan(f.fetch("chr2R", normalization=norm).mean(keep_nans=True, keep_infs=True))

        assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).mean(), 0.0003800524166065285)
        assert isclose(f.fetch("chr2L", "chr2R", normalization=norm).mean(keep_infs=True), 0.0003800524166065285)
        assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).mean(keep_nans=True))
        assert isnan(f.fetch("chr2L", "chr2R", normalization=norm).mean(keep_nans=True, keep_infs=True))
