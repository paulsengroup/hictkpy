# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import pathlib

import pytest

import hictkpy

from .helpers import get_test_dir, pandas_avail, pyarrow_avail

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (get_test_dir() / "data" / "cooler_test_file.mcool", 100_000),
        (get_test_dir() / "data" / "cooler_variable_bins_test_file.cool", None),
    ],
)


@pytest.mark.skipif(not pandas_avail() or not pyarrow_avail(), reason="either pandas or pyarrow are not available")
class TestFileCreationHiC:
    @staticmethod
    def setup_method():
        logging.basicConfig(level="INFO", force=True)
        logging.getLogger().setLevel("INFO")
        hictkpy.logging.setLevel("INFO")

    def test_accessors(self, file, resolution, tmpdir):
        with hictkpy.File(file, resolution) as f:
            bins = f.bins()

        if bins.type() != "fixed":
            pytest.skip(f'BinTable of file "{file}" does not have fixed bins.')

        path = tmpdir / "test.hic"
        with hictkpy.hic.FileWriter(path, bins) as w:
            assert str(w).startswith("HiCFileWriter(")
            assert w.path() == path
            assert w.resolutions() == [resolution]
            assert w.chromosomes() == bins.chromosomes()
            assert len(w.bins(resolution).to_df().compare(bins.to_df())) == 0

    def test_ctx_manager(self, file, resolution, tmpdir):
        import pandas as pd

        if resolution is None:
            pytest.skip()

        path = pathlib.Path(tmpdir) / "test.hic"
        with hictkpy.hic.FileWriter(path, {"chr1": 100}, 10) as w:
            pass

        assert w.path() == path  # Shouldn't throw
        assert path.is_file()

        with pytest.raises(RuntimeError, match=r"finalize\(\) was already called.*"):
            w.finalize()

        with pytest.raises(
            RuntimeError,
            match=r"caught attempt to add_pixels\(\) to a \.hic file that has already been finalized",
        ):
            w.add_pixels(pd.DataFrame())

        with pytest.raises(
            RuntimeError,
            match="caught an attempt to access file .*, which has already been closed",
        ):
            w.chromosomes()

        path.unlink(missing_ok=True)
        with hictkpy.hic.FileWriter(path, {"chr1": 100}, 10) as w:
            w.finalize()
            # leaving the context manager should not raise

        assert path.is_file()
        path.unlink(missing_ok=True)

        tmpdir_ = pathlib.Path(tmpdir) / "tmp"
        tmpdir_.mkdir()
        with pytest.raises(RuntimeError, match="foo"):
            with hictkpy.hic.FileWriter(path, {"chr1": 100}, 10, tmpdir=tmpdir_):
                raise RuntimeError("foo")
        assert len(list(tmpdir_.iterdir())) == 0
        assert not hictkpy.is_hic(path)

    def test_finalizer_warnings(self, file, resolution, tmpdir):
        path = tmpdir / "test.hic"
        with hictkpy.cooler.FileWriter(path, {"chr1": 100, "chr2": 50}, 10) as w:
            with pytest.deprecated_call():
                w.finalize("INFO")

    def test_file_creation_empty(self, file, resolution, tmpdir):
        if resolution is None:
            pytest.skip()

        path = tmpdir / "test.hic"
        with hictkpy.hic.FileWriter(path, {"chr1": 100, "chr2": 50}, 10) as w:
            pass

        with hictkpy.File(path) as f:
            assert f.fetch().nnz() == 0

    def test_file_creation_coo(self, file, resolution, tmpdir):
        with hictkpy.File(file, resolution) as f:
            if f.bins().type() != "fixed":
                pytest.skip(f'BinTable of file "{file}" does not have fixed bins.')

            chroms = f.chromosomes()
            df = f.fetch(join=False).to_df()
            expected_sum = df["count"].sum()

        path = tmpdir / "test.hic"
        chunk_size = 1000
        with hictkpy.hic.FileWriter(path, chroms, resolution) as w:
            for start in range(0, len(df), chunk_size):
                end = start + chunk_size
                w.add_pixels(df[start:end])

            w.finalize()
            with pytest.raises(
                RuntimeError,
                match=r"caught attempt to add_pixels\(\) to a \.hic file that has already been finalized",
            ):
                w.add_pixels(df)
            with pytest.raises(RuntimeError, match=r"finalize\(\) was already called.*"):
                w.finalize()

        with hictkpy.File(path) as f:
            assert f.fetch().sum() == expected_sum

    def test_file_creation_bg2(self, file, resolution, tmpdir):
        with hictkpy.File(file, resolution) as f:
            if f.bins().type() != "fixed":
                pytest.skip(f'BinTable of file "{file}" does not have fixed bins.')

            chroms = f.chromosomes()
            df = f.fetch(join=True).to_df()
            expected_sum = df["count"].sum()

        path = tmpdir / "test.hic"
        chunk_size = 1000
        with hictkpy.hic.FileWriter(path, chroms, resolution) as w:
            for start in range(0, len(df), chunk_size):
                end = start + chunk_size
                w.add_pixels(df[start:end])

            w.finalize()
            with pytest.raises(
                RuntimeError,
                match=r"caught attempt to add_pixels\(\) to a \.hic file that has already been finalized",
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

        path = tmpdir / "test.hic"
        if bins.type() != "fixed":
            with pytest.raises(RuntimeError):
                hictkpy.hic.FileWriter(path, bins)
            return

        chunk_size = 1000
        with hictkpy.hic.FileWriter(path, bins) as w:
            for start in range(0, len(df), chunk_size):
                end = start + chunk_size
                w.add_pixels(df[start:end])

            w.finalize()
            with pytest.raises(
                RuntimeError,
                match=r"caught attempt to add_pixels\(\) to a \.hic file that has already been finalized",
            ):
                w.add_pixels(df)
            with pytest.raises(RuntimeError, match=r"finalize\(\) was already called.*"):
                w.finalize()

        with hictkpy.File(path) as f:
            assert f.fetch().sum() == expected_sum
