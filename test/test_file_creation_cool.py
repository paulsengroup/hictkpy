# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

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
        with hictkpy.File(file, resolution) as f:
            bins = f.bins()

        path = tmpdir / "test.cool"

        with hictkpy.cooler.FileWriter(path, bins) as w:
            assert str(w).startswith("CoolFileWriter(")
            assert w.path() == path
            if resolution is None:
                assert w.resolution() == 0
            else:
                assert w.resolution() == resolution
            assert w.chromosomes() == bins.chromosomes()
            assert len(w.bins().to_df().compare(bins.to_df())) == 0

    def test_ctx_manager(self, file, resolution, tmpdir):
        import pandas as pd

        if resolution is None:
            pytest.skip()

        path = pathlib.Path(tmpdir) / "test.cool"
        with hictkpy.cooler.FileWriter(path, {"chr1": 100}, 10) as w:
            pass

        assert w.path() == path  # Shouldn't throw
        assert path.is_file()

        with pytest.raises(RuntimeError, match=r"finalize\(\) was already called.*"):
            w.finalize()

        with pytest.raises(
            RuntimeError,
            match=r"caught attempt to add_pixels\(\) to a \.cool file that has already been finalized",
        ):
            w.add_pixels(pd.DataFrame())

        with pytest.raises(
            RuntimeError,
            match="caught an attempt to access file .*, which has already been closed",
        ):
            w.chromosomes()

        path.unlink(missing_ok=True)
        with hictkpy.cooler.FileWriter(path, {"chr1": 100}, 10) as w:
            w.finalize()
            # leaving the context manager should not raise

        assert path.is_file()
        path.unlink(missing_ok=True)

        tmpdir_ = pathlib.Path(tmpdir) / "tmp"
        tmpdir_.mkdir()
        with pytest.raises(RuntimeError, match="foo"):
            with hictkpy.cooler.FileWriter(path, {"chr1": 100}, 10, tmpdir=tmpdir_):
                raise RuntimeError("foo")
        assert len(list(tmpdir_.iterdir())) == 0
        assert not hictkpy.is_cooler(path)

    def test_file_creation_empty(self, file, resolution, tmpdir):
        if resolution is None:
            pytest.skip()

        path = tmpdir / "test.cool"
        with hictkpy.cooler.FileWriter(path, {"chr1": 100, "chr2": 50}, 10):
            pass

        with hictkpy.File(path) as f:
            assert f.fetch().nnz() == 0

    def test_file_creation_thin_pixel(self, file, resolution, tmpdir):
        with hictkpy.File(file, resolution) as f:
            if f.bins().type() != "fixed":
                pytest.skip(f'BinTable of file "{file}" does not have fixed bins.')

            chroms = f.chromosomes()
            df = f.fetch(join=False).to_df()
            expected_sum = df["count"].sum()

        path = tmpdir / "test.cool"
        chunk_size = 1000
        with hictkpy.cooler.FileWriter(path, chroms, resolution) as w:
            for start in range(0, len(df), chunk_size):
                end = start + chunk_size
                w.add_pixels(df[start:end])

            w.finalize("info", 100_000, 100_000)
            with pytest.raises(
                RuntimeError,
                match=r"caught attempt to add_pixels\(\) to a \.cool file that has already been finalized",
            ):
                w.add_pixels(df)
            with pytest.raises(RuntimeError, match=r"finalize\(\) was already called.*"):
                w.finalize()

        with hictkpy.File(path) as f:
            assert f.fetch().sum() == expected_sum

    def test_file_creation(self, file, resolution, tmpdir):
        with hictkpy.File(file, resolution) as f:
            if f.bins().type() != "fixed":
                pytest.skip(f'BinTable of file "{file}" does not have fixed bins.')

            chroms = f.chromosomes()
            df = f.fetch(join=True).to_df()
            expected_sum = df["count"].sum()

        path = tmpdir / "test.cool"
        chunk_size = 1000
        with hictkpy.cooler.FileWriter(path, chroms, resolution) as w:
            for start in range(0, len(df), chunk_size):
                end = start + chunk_size
                w.add_pixels(df[start:end])

            w.finalize("info", 100_000, 100_000)
            with pytest.raises(
                RuntimeError,
                match=r"caught attempt to add_pixels\(\) to a \.cool file that has already been finalized",
            ):
                w.add_pixels(df)
            with pytest.raises(RuntimeError, match=r"finalize\(\) was already called.*"):
                w.finalize()

        with hictkpy.File(path) as f:
            assert f.fetch().sum() == expected_sum

    def test_file_creation_bin_table(self, file, resolution, tmpdir):
        with hictkpy.File(file, resolution) as f:
            bins = f.bins()
            df = f.fetch(join=True).to_df()
            expected_sum = df["count"].sum()

        path = tmpdir / "test.cool"
        chunk_size = 1000
        with hictkpy.cooler.FileWriter(path, bins) as w:
            for start in range(0, len(df), chunk_size):
                end = start + chunk_size
                w.add_pixels(df[start:end])

            w.finalize("info", 100_000, 100_000)
            with pytest.raises(
                RuntimeError,
                match=r"caught attempt to add_pixels\(\) to a \.cool file that has already been finalized",
            ):
                w.add_pixels(df)
            with pytest.raises(RuntimeError, match=r"finalize\(\) was already called.*"):
                w.finalize()

        with hictkpy.File(path) as f:
            assert f.fetch().sum() == expected_sum

    def test_file_creation_float_counts(self, file, resolution, tmpdir):
        with hictkpy.File(file, resolution) as f:
            if f.bins().type() != "fixed":
                pytest.skip(f'BinTable of file "{file}" does not have fixed bins.')

            chroms = f.chromosomes()
            df = f.fetch(join=True, count_type=float).to_df()
            df["count"] += 0.12345
            expected_sum = df["count"].sum()

        path = tmpdir / "test.cool"
        chunk_size = 1000

        with hictkpy.cooler.FileWriter(path, chroms, resolution) as w:
            for start in range(0, len(df), chunk_size):
                end = start + chunk_size
                w.add_pixels(df[start:end])

            w.finalize("info", 100_000, 100_000)
            with pytest.raises(
                RuntimeError,
                match=r"caught attempt to add_pixels\(\) to a \.cool file that has already been finalized",
            ):
                w.add_pixels(df)
            with pytest.raises(RuntimeError, match=r"finalize\(\) was already called.*"):
                w.finalize()

        with hictkpy.File(path) as f:
            assert pytest.approx(f.fetch(count_type="float").sum()) == expected_sum
