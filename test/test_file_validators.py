# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

import pytest

import hictkpy

testdir = pathlib.Path(__file__).resolve().parent

cool_file = (testdir / "data" / "cooler_test_file.mcool::/resolutions/100000").as_posix()
mcool_file = testdir / "data" / "cooler_test_file.mcool"
scool_file = testdir / "data" / "cooler_test_file.scool"
hic_file = testdir / "data" / "hic_test_file.hic"


class TestClass:
    def test_validators(self):
        assert hictkpy.is_cooler(cool_file)
        assert not hictkpy.is_cooler(hic_file)

        assert hictkpy.is_mcool_file(mcool_file)
        assert not hictkpy.is_mcool_file(cool_file)

        assert hictkpy.is_scool_file(scool_file)
        assert not hictkpy.is_scool_file(cool_file)
