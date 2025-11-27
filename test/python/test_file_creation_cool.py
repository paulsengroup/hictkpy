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
class TestFileCreationCooler:
    @staticmethod
    def setup_method():
        logging.basicConfig(level="INFO", force=True)
        logging.getLogger().setLevel("INFO")
        hictkpy.logging.setLevel("INFO")

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

    def test_finalizer_warnings(self, file, resolution, tmpdir):
        path = tmpdir / "test.cool"
        with hictkpy.cooler.FileWriter(path, {"chr1": 100, "chr2": 50}, 10) as w:
            with pytest.deprecated_call():
                w.finalize("INFO")

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

            w.finalize()
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

            w.finalize()
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

            w.finalize()
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

            w.finalize()
            with pytest.raises(
                RuntimeError,
                match=r"caught attempt to add_pixels\(\) to a \.cool file that has already been finalized",
            ):
                w.add_pixels(df)
            with pytest.raises(RuntimeError, match=r"finalize\(\) was already called.*"):
                w.finalize()

        with hictkpy.File(path) as f:
            assert pytest.approx(f.fetch(count_type="float").sum()) == expected_sum

    def test_non_uniform_dtypes_coo(self, file, resolution, tmpdir):
        import pyarrow as pa

        chroms = {"chr1": 10, "chr2": 5}
        resolution = 2
        path = tmpdir / "test.cool"

        df = pa.Table.from_pydict(
            {
                "bin1_id": pa.array([0], type=pa.uint64()),
                "bin2_id": pa.array([0], type=pa.int64()),
                "count": pa.array([1], type=pa.float64()),
            }
        )

        with hictkpy.cooler.FileWriter(path, chroms, resolution) as w:
            w.add_pixels(df)

        with hictkpy.File(path) as f:
            assert f.fetch().nnz() == 1

    def test_non_uniform_dtypes_bg2(self, file, resolution, tmpdir):
        import pyarrow as pa

        chroms = {"chr1": 10, "chr2": 5}
        resolution = 2
        path = tmpdir / "test.cool"

        df = pa.Table.from_pydict(
            {
                "chrom1": pa.array(["chr1"], type=pa.string()),
                "start1": pa.array([0], type=pa.uint8()),
                "end1": pa.array([resolution], type=pa.int16()),
                "chrom2": pa.array(["chr1"], type=pa.large_string()),
                "start2": pa.array([0], type=pa.int32()),
                "end2": pa.array([resolution], type=pa.uint64()),
                "count": pa.array([1], type=pa.float64()),
            }
        )

        with hictkpy.cooler.FileWriter(path, chroms, resolution) as w:
            w.add_pixels(df)

        with hictkpy.File(path) as f:
            assert f.fetch().nnz() == 1
        pass

    def test_invalid_pixels_coo(self, file, resolution, tmpdir):
        import pandas as pd
        import pyarrow as pa

        chroms = {"chr1": 10, "chr2": 5}
        resolution = 2
        path = tmpdir / "test.cool"

        with hictkpy.cooler.FileWriter(path, chroms, resolution) as w:
            with pytest.raises(
                ValueError, match=r"expected table to be of type pandas\.DataFrame or pyarrow\.Table, found .*"
            ):
                w.add_pixels({})

            df = pd.DataFrame({"a": [0]})
            with pytest.raises(ValueError, match="DataFrame is not in COO or BG2 format.*"):
                w.add_pixels(df)

            df = pa.Table.from_pandas(df)
            with pytest.raises(ValueError, match="DataFrame is not in COO or BG2 format.*"):
                w.add_pixels(df)

            df = pd.DataFrame({"bin1_id": [0], "bin2_id": [0], "count": [1], "count_": [2]}).rename(
                columns={"count_": "count"}
            )
            with pytest.raises(ValueError, match="Duplicate column names found.*"):
                w.add_pixels(df)

            df = pd.DataFrame({"bin1_id": [-1], "bin2_id": [0], "count": [1]})
            with pytest.raises(ValueError, match=".*found negative value in bin1_id column"):
                w.add_pixels(df)

    def test_invalid_pixels_bg2(self, file, resolution, tmpdir):
        import pandas as pd

        chroms = {"chr1": 10, "chr2": 5}
        resolution = 2
        path = tmpdir / "test.cool"

        with hictkpy.cooler.FileWriter(path, chroms, resolution) as w:
            df = pd.DataFrame(
                {
                    "chrom1": ["chr1"],
                    "start1": [0],
                    "end1": [2],
                    "chrom2": ["chr1"],
                    "start2": [0],
                    "end2": [2],
                    "count": [1],
                    "chrom1_": ["0"],
                }
            ).rename(columns={"chrom1_": "chrom1"})

            with pytest.raises(ValueError, match="Duplicate column names found.*"):
                w.add_pixels(df)

            df = pd.DataFrame(
                {
                    "chrom1": ["chr1"],
                    "start1": [0],
                    "end1": [2],
                    "chrom2": ["chr1"],
                    "start2": [0],
                    "end2": [2],
                    "count": [1],
                }
            )

            for col in ("start1", "end1", "start2", "end2"):
                dff = df.copy()
                dff[col] = [-1]
                with pytest.raises(ValueError, match="genomic coordinates cannot be negative"):
                    w.add_pixels(dff)

            for col in ("start1", "start2"):
                dff = df.copy()
                dff[col] = [3]
                with pytest.raises(ValueError, match="end position of a bin cannot be smaller than its start position"):
                    w.add_pixels(dff)

            df["end1"] = [3]
            with pytest.raises(ValueError, match=r".*end position of a bin cannot be smaller than its start position"):
                w.add_pixels(dff)

            df = pd.DataFrame(
                {
                    "chrom1": [0],
                    "start1": [0],
                    "end1": [2],
                    "chrom2": [0],
                    "start2": [0],
                    "end2": [2],
                    "count": [1],
                }
            )

            with pytest.raises(ValueError, match="DataFrame is not in COO or BG2 format.*"):
                w.add_pixels(df)

            # categorical have incompatible underlying type
            df["chrom1"] = pd.Categorical(df["chrom1"])
            df["chrom2"] = pd.Categorical(df["chrom2"])

            with pytest.raises(ValueError, match="DataFrame is not in COO or BG2 format.*"):
                w.add_pixels(df)
