# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import gc
import logging
import pathlib

import pytest

import hictkpy

testdir = pathlib.Path(__file__).resolve().parent


pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (testdir / "data" / "hic_test_file.hic", 100_000),
    ],
)


def pandas_avail() -> bool:
    try:
        import pandas
    except ModuleNotFoundError:
        return False

    return True


def pyarrow_avail() -> bool:
    try:
        import pyarrow
    except ModuleNotFoundError:
        return False

    return True


@pytest.mark.skipif(not pandas_avail() or not pyarrow_avail(), reason="either pandas or pyarrow are not available")
class TestClass:
    @staticmethod
    def setup_method():
        logging.basicConfig(level="INFO", force=True)
        logging.getLogger().setLevel("INFO")

    def test_file_creation_thin_pixel(self, file, resolution, tmpdir):
        f = hictkpy.File(file, resolution)

        df = f.fetch(join=False).to_df()
        expected_sum = df["count"].sum()

        path = tmpdir / "test1.hic"
        w = hictkpy.hic.FileWriter(path, f.chromosomes(), f.resolution())

        chunk_size = 1000
        for start in range(0, len(df), chunk_size):
            end = start + chunk_size
            w.add_pixels(df[start:end])

        w.finalize()
        with pytest.raises(Exception):
            w.add_pixels(df)
        with pytest.raises(Exception):
            w.finalize()

        del w
        gc.collect()

        f = hictkpy.File(path, resolution)
        assert f.fetch().sum() == expected_sum

    def test_file_creation(self, file, resolution, tmpdir):
        f = hictkpy.File(file, resolution)

        df = f.fetch(join=True).to_df()
        expected_sum = df["count"].sum()

        path = tmpdir / "test2.hic"
        w = hictkpy.hic.FileWriter(path, f.chromosomes(), f.resolution())

        chunk_size = 1000
        for start in range(0, len(df), chunk_size):
            end = start + chunk_size
            w.add_pixels(df[start:end])

        w.finalize()
        with pytest.raises(Exception):
            w.add_pixels(df)
        with pytest.raises(Exception):
            w.finalize()

        del w
        gc.collect()

        f = hictkpy.File(path, resolution)
        assert f.fetch().sum() == expected_sum
