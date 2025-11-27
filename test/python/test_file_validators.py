# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

import hictkpy

from .helpers import get_test_dir

cool_file = (get_test_dir() / "data" / "cooler_test_file.mcool::/resolutions/100000").as_posix()
mcool_file = get_test_dir() / "data" / "cooler_test_file.mcool"
scool_file = get_test_dir() / "data" / "cooler_test_file.scool"
hic_file = get_test_dir() / "data" / "hic_test_file.hic"


class TestFileValidators:
    def test_valid_formats(self):
        assert hictkpy.is_cooler(cool_file)
        assert not hictkpy.is_cooler(hic_file)

        assert hictkpy.is_mcool_file(mcool_file)
        assert not hictkpy.is_mcool_file(cool_file)

        assert hictkpy.is_scool_file(scool_file)
        assert not hictkpy.is_scool_file(cool_file)

        assert hictkpy.is_hic(hic_file)
        assert not hictkpy.is_hic(cool_file)

    def test_invalid_formats(self):
        path = pathlib.Path(__file__).resolve()

        assert not hictkpy.is_cooler(path)
        assert not hictkpy.is_mcool_file(path)
        assert not hictkpy.is_scool_file(path)
        assert not hictkpy.is_hic(path)

    def test_invalid_files(self):
        non_existing_file = get_test_dir() / "foobar.123"
        assert not non_existing_file.exists()

        assert not hictkpy.is_cooler(non_existing_file)
        assert not hictkpy.is_mcool_file(non_existing_file)
        assert not hictkpy.is_scool_file(non_existing_file)
        assert not hictkpy.is_hic(non_existing_file)

        folder = get_test_dir()
        assert folder.is_dir()

        assert not hictkpy.is_cooler(folder)
        assert not hictkpy.is_mcool_file(folder)
        assert not hictkpy.is_scool_file(folder)
        assert not hictkpy.is_hic(folder)
