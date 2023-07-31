# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os

import numpy as np
import pytest

import hictkpy

testdir = os.path.dirname(os.path.abspath(__file__))


def fetch_and_compare(f):
    m = f.fetch_sparse()
    assert m.shape == (1380, 1380)
    assert m.sum() == 119_208_613

    ### CIS
    m = f.fetch_sparse("chr2R:10,000,000-15,000,000")
    assert m.shape == (50, 50)
    assert m.sum() == 4_519_080

    m = f.fetch_sparse("chr2R:10,000,000-15,000,000", count_type="int")
    assert m.dtype == np.int32

    m = f.fetch_sparse("chr2R:10,000,000-15,000,000", count_type="float")
    assert m.dtype == np.float64

    m = f.fetch_sparse("chr2R\t10000000\t15000000", query_type="BED")
    assert m.shape == (50, 50)

    ### TRANS
    m = f.fetch_sparse("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000")
    assert m.shape == (50, 100)
    assert m.sum() == 83_604

    m = f.fetch_sparse("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="int")
    assert m.dtype == np.int32

    m = f.fetch_sparse("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="float")
    assert m.dtype == np.float64

    m = f.fetch_sparse("chr2R\t10000000\t15000000", "chrX\t0\t10000000", query_type="BED")
    assert m.shape == (50, 100)


@pytest.mark.parametrize(
    "file,resolution",
    [(os.path.join(testdir, "data", "cooler_test_file.cool"), 100_000)],
)
def test_file_fetch_sparse_file(file, resolution):
    f = hictkpy.File(file, resolution)
    fetch_and_compare(f)

    df = f.fetch("chr2R:10,000,000-15,000,000", normalization="weight")
    assert np.isclose(59.349524704033215, df["count"].sum())


@pytest.mark.parametrize(
    "file",
    [(os.path.join(testdir, "data", "cooler_test_file.cool"))],
)
def test_cooler_fetch_sparse_cooler(file):
    f = hictkpy.cooler.File(file)
    fetch_and_compare(f)

    df = f.fetch("chr2R:10,000,000-15,000,000", normalization="weight")
    assert np.isclose(59.349524704033215, df["count"].sum())


@pytest.mark.parametrize(
    "file,resolution",
    [(os.path.join(testdir, "data", "hic_test_file.hic"), 100_000)],
)
def test_hic_fetch_sparse_hic(file, resolution):
    f = hictkpy.hic.File(file, resolution)
    fetch_and_compare(f)

    df = f.fetch("chr2R:10,000,000-15,000,000", normalization="ICE")
    assert np.isclose(59.349524704033215, df["count"].sum())
