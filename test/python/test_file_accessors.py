# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import warnings

import pytest

import hictkpy

from .helpers import get_test_dir, numpy_avail, pandas_avail, pyarrow_avail

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        ((get_test_dir() / "data" / "cooler_test_file.mcool::/resolutions/100000").as_posix(), None),
        (get_test_dir() / "data" / "cooler_test_file.mcool", 100_000),
        (get_test_dir() / "data" / "hic_test_file.hic", 100_000),
    ],
)


class TestFileAccessors:
    def test_attributes(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            assert f.resolution() == 100_000
            assert f.nchroms() == 8
            assert f.nbins() == 1380

            assert "chr2L" in f.chromosomes()

            assert len(f.bins()) == 1380
            assert len(f.chromosomes()) == 8

            if f.is_cooler():
                assert f.attributes()["format"] == "HDF5::Cooler"
            else:
                assert f.attributes()["format"] == "HIC"

    @pytest.mark.skipif(not pandas_avail() or not pyarrow_avail(), reason="either pandas or pyarrow are not available")
    def test_normalizations(self, file, resolution):
        cooler_weights = ["KR", "SCALE", "VC", "VC_SQRT", "weight"]
        hic_weights = ["ICE"]

        with hictkpy.File(file, resolution) as f:
            if f.is_cooler():
                assert f.avail_normalizations() == cooler_weights
            else:
                assert f.avail_normalizations() == hic_weights

            assert not f.has_normalization("foo")

            name = "weight" if f.is_cooler() else "ICE"
            weights = f.weights(name)
            assert len(weights) == f.nbins()

            if f.is_cooler():
                df = f.weights(cooler_weights)
                assert len(df.columns) == len(cooler_weights)
            else:
                df = f.weights(hic_weights)
                assert len(df.columns) == len(hic_weights)
            assert len(df) == f.nbins()

    @pytest.mark.skipif(not numpy_avail(), reason="numpy is not available")
    def test_fetch(self, file, resolution):
        import numpy as np

        valid_types = {
            int: np.int64,
            float: np.float64,
            np.uint8: np.uint8,
            np.uint16: np.uint16,
            np.uint32: np.uint32,
            np.uint64: np.uint64,
            np.int8: np.int8,
            np.int16: np.int16,
            np.int32: np.int32,
            np.int64: np.int64,
            "int": np.int32,
            "uint": np.uint32,
            "float": np.float64,
            "uint8": np.uint8,
            "uint16": np.uint16,
            "uint32": np.uint32,
            "uint64": np.uint64,
            "int8": np.int8,
            "int16": np.int16,
            "int32": np.int32,
            "int64": np.int64,
        }

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            invalid_types = [
                object,
                bool,
                str,
                np.array,
                np.bool_,
                np.float16,
                np.complex64,
                "str",
                "foobar",
            ]

            try:
                invalid_types.append(np.float128)
            except:  # noqa
                pass

        with hictkpy.File(file, resolution) as f:
            assert f.fetch().dtype() is np.int32

            for t1, t2 in valid_types.items():
                assert f.fetch(count_type=t1).dtype() is t2

            for t in invalid_types:
                with pytest.raises(TypeError):
                    f.fetch(count_type=t)
