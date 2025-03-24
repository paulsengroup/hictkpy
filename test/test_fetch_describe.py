# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
import tempfile
from math import inf, isnan

import pytest

import hictkpy

from .helpers import numpy_avail, pandas_avail


def isclose(n1, n2) -> bool:
    import math

    return math.isclose(n1, n2, rel_tol=1.0e-4, abs_tol=1.0e-8)


@pytest.mark.skipif(not numpy_avail() or not pandas_avail(), reason="either numpy or pandas are not available")
class TestClass:
    @staticmethod
    def generate_pixels(insert_nan: bool = False, insert_neg_inf: bool = False, insert_pos_inf: bool = False):
        import numpy as np
        import pandas as pd

        chroms = {"chr1": 1500}
        resolution = 10
        num_bins = chroms["chr1"] // resolution

        bins = list(range(0, num_bins))

        bin1_ids = []
        bin2_ids = []

        # generate coordinates
        for i in bins:
            bin1_ids.extend([i] * len(bins))
            bin2_ids.extend(bins)

        bin1_ids = np.array(bin1_ids)
        bin2_ids = np.array(bin2_ids)

        # drop coordinates overlapping with the lower-triangular matrix
        mask = bin1_ids <= bin1_ids
        bin1_ids = bin1_ids[mask]
        bin2_ids = bin2_ids[mask]

        # make matrix sparse
        bin1_ids = bin1_ids[::2]
        bin2_ids = bin2_ids[::2]

        counts = [i + 0.123 for i, _ in enumerate(bin1_ids)]

        if insert_nan:
            counts[-1] = np.nan
        if insert_neg_inf:
            counts[-2] = -np.inf
        if insert_pos_inf:
            counts[-3] = np.inf

        return chroms, resolution, pd.DataFrame({"bin1_id": bin1_ids, "bin2_id": bin2_ids, "count": counts})

    @staticmethod
    def make_cooler_file(chroms, resolution, pixels, tmpdir) -> hictkpy.File:
        with tempfile.NamedTemporaryFile(dir=tmpdir, suffix=".cool") as tmpfile:
            path = pathlib.Path(tmpfile.name)

        w = hictkpy.cooler.FileWriter(path, chroms, resolution)
        w.add_pixels(pixels)
        w.finalize()

        return hictkpy.File(path, resolution)

    def test_describe_all_finite(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=False, insert_neg_inf=False, insert_pos_inf=False), tmpdir
        )

        stats = f.fetch(count_type="float").describe(keep_nans=True, keep_infs=True)

        assert stats.get("nnz", -1) == 11_250
        assert isclose(stats.get("sum", -1), 63277008.749)
        assert isclose(stats.get("min", -1), 0.123)
        assert isclose(stats.get("max", -1), 11249.123)
        assert isclose(stats.get("mean", -1), 5624.6229)
        assert isclose(stats.get("variance", -1), 10547812.5)
        assert isclose(stats.get("skewness", -1), 1.7841628283203783e-15)
        assert isclose(stats.get("kurtosis", -1), -1.2000000189629632)

    def test_describe_subset_finite(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=False, insert_neg_inf=False, insert_pos_inf=False), tmpdir
        )

        stats = f.fetch(count_type="float").describe(["nnz", "kurtosis"])
        assert len(stats) == 2
        assert stats.get("nnz", -1) == 11_250
        assert isclose(stats.get("kurtosis", -1), -1.2000000189629632)

        stats = f.fetch(count_type="float").describe(["max", "variance"])
        assert len(stats) == 2
        assert isclose(stats.get("max", -1), 11249.123)
        assert isclose(stats.get("variance", -1), 10547812.5)

    def test_describe_all_with_nans(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=True, insert_neg_inf=False, insert_pos_inf=False), tmpdir
        )

        stats = f.fetch(count_type="float").describe(keep_nans=True, keep_infs=True)

        assert stats.get("nnz", -1) == 11_250
        assert isnan(stats.get("sum", -1))
        assert isnan(stats.get("min", -1))
        assert isnan(stats.get("max", -1))
        assert isnan(stats.get("mean", -1))
        assert isnan(stats.get("variance", -1))
        assert isnan(stats.get("skewness", -1))
        assert isnan(stats.get("kurtosis", -1))

    def test_describe_subset_with_nans(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=True, insert_neg_inf=False, insert_pos_inf=False), tmpdir
        )

        stats = f.fetch(count_type="float").describe(["nnz", "kurtosis"], keep_nans=True, keep_infs=True)
        assert len(stats) == 2
        assert stats.get("nnz", -1) == 11_250
        assert isnan(stats.get("kurtosis", -1))

        stats = f.fetch(count_type="float").describe(["max", "variance"], keep_nans=True, keep_infs=True)
        assert len(stats) == 2
        assert isnan(stats.get("max", -1))
        assert isnan(stats.get("variance", -1))

    def test_describe_all_with_neg_inf(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=False, insert_neg_inf=True, insert_pos_inf=False), tmpdir
        )

        stats = f.fetch(count_type="float").describe(keep_nans=True, keep_infs=True)

        assert stats.get("nnz", -1) == 11_250
        assert stats.get("sum", -1) == -inf
        assert stats.get("min", -1) == -inf
        assert isclose(stats.get("max", -1), 11249.123)
        assert stats.get("mean", -1) == -inf
        assert isnan(stats.get("variance", -1))
        assert isnan(stats.get("skewness", -1))
        assert isnan(stats.get("kurtosis", -1))

    def test_describe_subset_with_neg_inf(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=False, insert_neg_inf=True, insert_pos_inf=False), tmpdir
        )

        stats = f.fetch(count_type="float").describe(["nnz", "kurtosis"], keep_nans=True, keep_infs=True)
        assert len(stats) == 2
        assert stats.get("nnz", -1) == 11_250
        assert isnan(stats.get("kurtosis", -1))

        stats = f.fetch(count_type="float").describe(["max", "variance"], keep_nans=True, keep_infs=True)
        assert len(stats) == 2
        assert isclose(stats.get("max", -1), 11249.123)
        assert isnan(stats.get("variance", -1))

    def test_describe_all_with_pos_inf(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=False, insert_neg_inf=False, insert_pos_inf=True), tmpdir
        )

        stats = f.fetch(count_type="float").describe(keep_nans=True, keep_infs=True)

        assert stats.get("nnz", -1) == 11_250
        assert stats.get("sum", -1) == inf
        assert isclose(stats.get("min", -1), 0.123)
        assert stats.get("max", -1) == inf
        assert stats.get("mean", -1) == inf
        assert isnan(stats.get("variance", -1))
        assert isnan(stats.get("skewness", -1))
        assert isnan(stats.get("kurtosis", -1))

    def test_describe_subset_with_pos_inf(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=False, insert_neg_inf=False, insert_pos_inf=True), tmpdir
        )

        stats = f.fetch(count_type="float").describe(["nnz", "kurtosis"], keep_nans=True, keep_infs=True)
        assert len(stats) == 2
        assert stats.get("nnz", -1) == 11_250
        assert isnan(stats.get("kurtosis", -1))

        stats = f.fetch(count_type="float").describe(["max", "variance"], keep_nans=True, keep_infs=True)
        assert len(stats) == 2
        assert stats.get("max", -1) == inf
        assert isnan(stats.get("variance", -1))

    def test_describe_all_with_neg_and_pos_inf(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=False, insert_neg_inf=True, insert_pos_inf=True), tmpdir
        )

        stats = f.fetch(count_type="float").describe(keep_nans=True, keep_infs=True)

        assert stats.get("nnz", -1) == 11_250
        assert isnan(stats.get("sum", -1))
        assert stats.get("min", -1) == -inf
        assert stats.get("max", -1) == inf
        assert isnan(stats.get("mean", -1))
        assert isnan(stats.get("variance", -1))
        assert isnan(stats.get("skewness", -1))
        assert isnan(stats.get("kurtosis", -1))

    def test_describe_subset_with_neg_and_pos_inf(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=False, insert_neg_inf=True, insert_pos_inf=True), tmpdir
        )

        stats = f.fetch(count_type="float").describe(["nnz", "kurtosis"], keep_nans=True, keep_infs=True)
        assert len(stats) == 2
        assert stats.get("nnz", -1) == 11_250
        assert isnan(stats.get("kurtosis", -1))

        stats = f.fetch(count_type="float").describe(["max", "variance"], keep_nans=True, keep_infs=True)
        assert len(stats) == 2
        assert stats.get("max", -1) == inf
        assert isnan(stats.get("variance", -1))

    def test_describe_all_with_nan_and_inf(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=True, insert_neg_inf=False, insert_pos_inf=True), tmpdir
        )

        stats = f.fetch(count_type="float").describe(keep_nans=True, keep_infs=True)

        assert stats.get("nnz", -1) == 11_250
        assert isnan(stats.get("sum", -1))
        assert isnan(stats.get("min", -1))
        assert isnan(stats.get("max", -1))
        assert isnan(stats.get("mean", -1))
        assert isnan(stats.get("variance", -1))
        assert isnan(stats.get("skewness", -1))
        assert isnan(stats.get("kurtosis", -1))

    def test_describe_subset_with_nan_and_inf(self, tmpdir):
        f = self.make_cooler_file(
            *self.generate_pixels(insert_nan=True, insert_neg_inf=False, insert_pos_inf=True), tmpdir
        )

        stats = f.fetch(count_type="float").describe(["nnz", "kurtosis"], keep_nans=True, keep_infs=True)
        assert len(stats) == 2
        assert stats.get("nnz", -1) == 11_250
        assert isnan(stats.get("kurtosis", -1))

        stats = f.fetch(count_type="float").describe(["max", "variance"], keep_nans=True, keep_infs=True)
        assert len(stats) == 2
        assert isnan(stats.get("max", -1))
        assert isnan(stats.get("variance", -1))
