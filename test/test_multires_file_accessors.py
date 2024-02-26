# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os

import pytest

import hictkpy

testdir = os.path.dirname(os.path.abspath(__file__))

pytestmark = pytest.mark.parametrize(
    "file",
    [
        (os.path.join(testdir, "data", "cooler_test_file.mcool")),
    ],
)


class TestClass:
    def test_attributes(self, file):
        f = hictkpy.MultiResFile(file)

        assert f.path() == file
        assert f.resolutions() == [100_000, 1_000_000]
        assert len(f.chromosomes()) == 8

        assert f[100_000].resolution() == 100_000
