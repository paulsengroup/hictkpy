# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import gc
import logging
import pathlib

import pytest

import hictkpy

from .helpers import pandas_avail, pyarrow_avail

testdir = pathlib.Path(__file__).resolve().parent


pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (testdir / "data" / "cooler_test_file.mcool", 100_000),
        (testdir / "data" / "cooler_variable_bins_test_file.cool", None),
    ],
)


@pytest.mark.skipif(not pandas_avail() or not pyarrow_avail(), reason="either pandas or pyarrow are not available")
class TestClass:
    @staticmethod
    def setup_method():
        logging.basicConfig(level="INFO", force=True)
        logging.getLogger().setLevel("INFO")

    def test_accessors(self, file, resolution, tmpdir):
        bins = hictkpy.File(file, resolution).bins()
        if bins.type() != "fixed":
            pytest.skip(f'BinTable of file "{file}" does not have fixed bins.')

        path = tmpdir / "test.hic"
        w = hictkpy.hic.FileWriter(path, bins)

        assert str(w).startswith("HiCFileWriter(")
        assert w.path() == path
        assert w.resolutions() == [resolution]
        assert w.chromosomes() == bins.chromosomes()
        assert len(w.bins(resolution).to_df().compare(bins.to_df())) == 0

    def test_file_creation_empty(self, file, resolution, tmpdir):
        if resolution is None:
            pytest.skip()

        path = tmpdir / "test.hic"
        w = hictkpy.hic.FileWriter(path, {"chr1": 100, "chr2": 50}, 10)
        f = w.finalize("info")

        del w
        gc.collect()

        assert f.fetch().nnz() == 0

    def test_file_creation_thin_pixel(self, file, resolution, tmpdir):
        f = hictkpy.File(file, resolution)
        if f.bins().type() != "fixed":
            pytest.skip(f'BinTable of file "{file}" does not have fixed bins.')

        df = f.fetch(join=False).to_df()
        expected_sum = df["count"].sum()

        path = tmpdir / "test.hic"
        w = hictkpy.hic.FileWriter(path, f.chromosomes(), f.resolution())

        chunk_size = 1000
        for start in range(0, len(df), chunk_size):
            end = start + chunk_size
            w.add_pixels(df[start:end])

        f = w.finalize()
        with pytest.raises(Exception):
            w.add_pixels(df)
        with pytest.raises(Exception):
            w.finalize()

        del w
        gc.collect()

        assert f.fetch().sum() == expected_sum

    def test_file_creation(self, file, resolution, tmpdir):
        f = hictkpy.File(file, resolution)
        if f.bins().type() != "fixed":
            pytest.skip(f'BinTable of file "{file}" does not have fixed bins.')

        df = f.fetch(join=True).to_df()
        expected_sum = df["count"].sum()

        path = tmpdir / "test.hic"
        w = hictkpy.hic.FileWriter(path, f.chromosomes(), f.resolution())

        chunk_size = 1000
        for start in range(0, len(df), chunk_size):
            end = start + chunk_size
            w.add_pixels(df[start:end])

        f = w.finalize()
        with pytest.raises(Exception):
            w.add_pixels(df)
        with pytest.raises(Exception):
            w.finalize()

        del w
        gc.collect()

        assert f.fetch().sum() == expected_sum

    def test_file_creation_bin_table(self, file, resolution, tmpdir):
        f = hictkpy.File(file, resolution)

        df = f.fetch(join=True).to_df()
        expected_sum = df["count"].sum()

        path = tmpdir / "test.hic"
        if f.bins().type() != "fixed":
            with pytest.raises(Exception):
                hictkpy.hic.FileWriter(path, f.bins())
            return

        w = hictkpy.hic.FileWriter(path, f.bins())

        chunk_size = 1000
        for start in range(0, len(df), chunk_size):
            end = start + chunk_size
            w.add_pixels(df[start:end])

        f = w.finalize()
        with pytest.raises(Exception):
            w.add_pixels(df)
        with pytest.raises(Exception):
            w.finalize()

        del w
        gc.collect()

        assert f.fetch().sum() == expected_sum
