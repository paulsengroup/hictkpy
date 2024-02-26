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
        (os.path.join(testdir, "data", "cooler_test_file.scool")),
    ],
)


class TestClass:
    def test_attributes(self, file):
        f = hictkpy.cooler.SingleCellFile(file)

        assert f.path() == file
        assert f.resolution() == 100_000
        assert len(f.chromosomes()) == 20
        assert len(f.cells()) == 5

        assert f.attributes()["format"] == "HDF5::SCOOL"
        assert f["GSM2687248_41669_ACAGTG-R1-DpnII.100000.cool"].resolution() == 100_000
