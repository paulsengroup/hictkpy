# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pytest

import hictkpy

from .helpers import get_test_dir

pytestmark = pytest.mark.parametrize(
    "file",
    [
        get_test_dir() / "data" / "cooler_test_file.scool",
    ],
)


class TestSingleCellFileAccessors:
    def test_accessors(self, file):
        with hictkpy.cooler.SingleCellFile(file) as f:
            assert str(f).startswith("SingleCellFile(")

            assert f.path() == file
            assert f.resolution() == 100_000
            assert len(f.chromosomes()) == 20
            assert len(f.bins()) == 26398
            assert len(f.cells()) == 5

            assert f.attributes()["format"] == "HDF5::SCOOL"
            assert f["GSM2687248_41669_ACAGTG-R1-DpnII.100000.cool"].resolution() == 100_000

            with pytest.raises(Exception):
                f["ABC"]  # noqa
