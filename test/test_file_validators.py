# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os

import pytest

import hictkpy

testdir = os.path.dirname(os.path.abspath(__file__))

cool_file = os.path.join(testdir, "data", "cooler_test_file.mcool::/resolutions/100000")
mcool_file = os.path.join(testdir, "data", "cooler_test_file.mcool")
scool_file = os.path.join(testdir, "data", "cooler_test_file.scool")
hic_file = os.path.join(testdir, "data", "hic_test_file.hic")


class TestClass:
    def test_validators(self):
        assert hictkpy.is_cooler(cool_file)
        assert not hictkpy.is_cooler(hic_file)

        assert hictkpy.is_mcool_file(mcool_file)
        assert not hictkpy.is_mcool_file(cool_file)

        assert hictkpy.is_scool_file(scool_file)
        assert not hictkpy.is_scool_file(cool_file)
