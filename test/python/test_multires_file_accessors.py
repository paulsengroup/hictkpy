# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


import pytest

import hictkpy

from .helpers import get_test_dir

pytestmark = pytest.mark.parametrize(
    "file,format",
    [
        (get_test_dir() / "data" / "cooler_test_file.mcool", "mcool"),
        (get_test_dir() / "data" / "hic_test_file.hic", "hic"),
    ],
)


class TestMultiResFileAccessors:
    def test_accessors(self, file, format):
        with hictkpy.MultiResFile(file) as f:
            assert str(f).startswith("MultiResFile(")

            assert f.path() == file
            assert f.is_mcool() == (format == "mcool")
            assert f.is_hic() == (format == "hic")
            assert len(f.chromosomes()) == 8

            if f.is_hic():
                resolutions = [100_000]
                assert (f.resolutions() == resolutions).all()
                assert f.attributes()["format"] == "HIC"
                assert f.attributes()["format-version"] == 9
                assert (f.attributes()["resolutions"] == resolutions).all()
            else:
                resolutions = [100_000, 1_000_000]
                assert (f.resolutions() == resolutions).all()
                assert f.attributes()["format"] == "HDF5::MCOOL"
                assert f.attributes()["format-version"] == 2
                assert (f.attributes()["resolutions"] == resolutions).all()

            assert f[100_000].resolution() == 100_000

            with pytest.raises(Exception):
                f[1234]  # noqa
