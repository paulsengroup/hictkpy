# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import gc
import logging
import os
import tempfile

import pytest

import hictkpy

testdir = os.path.dirname(os.path.abspath(__file__))


pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (os.path.join(testdir, "data", "cooler_test_file.mcool"), 100_000),
    ],
)


class TestClass:
    @staticmethod
    def setup_method():
        logging.basicConfig(level="INFO", force=True)
        logging.getLogger().setLevel("INFO")

    def test_file_creation_thin_pixel(self, file, resolution, tmpdir):

        f = hictkpy.File(file, resolution)

        df = f.fetch(join=False).to_df()
        expected_sum = df["count"].sum()

        path = os.path.join(tmpdir, "test1.cool")
        w = hictkpy.cooler.FileWriter(path, f.chromosomes(), f.resolution())

        chunk_size = 1000
        for start in range(0, len(df), chunk_size):
            end = start + chunk_size
            w.add_pixels(df[start:end])

        w.finalize("info", 100_000, 100_000)
        del w
        gc.collect()

        f = hictkpy.File(path, resolution)
        assert f.fetch().sum() == expected_sum

    def test_file_creation(self, file, resolution, tmpdir):
        f = hictkpy.File(file, resolution)

        df = f.fetch(join=True).to_df()
        expected_sum = df["count"].sum()

        path = os.path.join(tmpdir, "test2.cool")
        w = hictkpy.cooler.FileWriter(path, f.chromosomes(), f.resolution())

        chunk_size = 1000
        for start in range(0, len(df), chunk_size):
            end = start + chunk_size
            w.add_pixels(df[start:end])

        w.finalize("info", 100_000, 100_000)
        del w
        gc.collect()

        f = hictkpy.File(path, resolution)
        assert f.fetch().sum() == expected_sum

    def test_file_creation_float_counts(self, file, resolution, tmpdir):
        f = hictkpy.File(file, resolution)

        df = f.fetch(join=True, count_type="float").to_df()
        df["count"] += 0.12345
        expected_sum = df["count"].sum()

        path = os.path.join(tmpdir, "test2.cool")
        w = hictkpy.cooler.FileWriter(path, f.chromosomes(), f.resolution())

        chunk_size = 1000
        for start in range(0, len(df), chunk_size):
            end = start + chunk_size
            w.add_pixels(df[start:end])

        w.finalize("info", 100_000, 100_000)
        del w
        gc.collect()

        f = hictkpy.File(path, resolution)
        assert pytest.approx(f.fetch(count_type="float").sum()) == expected_sum
