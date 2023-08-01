# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os

import pytest

import hictkpy

testdir = os.path.dirname(os.path.abspath(__file__))


def compare_nnz(f):
    assert f.fetch_nnz() == 890_384
    assert f.fetch_nnz("chr2R") == 31_900


@pytest.mark.parametrize(
    "file,resolution",
    [(os.path.join(testdir, "data", "cooler_test_file.cool"), 100_000)],
)
def test_file_fetch_nnz_file(file, resolution):
    f = hictkpy.File(file, resolution)
    compare_nnz(f)


@pytest.mark.parametrize(
    "file",
    [(os.path.join(testdir, "data", "cooler_test_file.cool"))],
)
def test_cooler_fetch_nnz_cooler(file):
    f = hictkpy.cooler.File(file)
    compare_nnz(f)


@pytest.mark.parametrize(
    "file,resolution",
    [(os.path.join(testdir, "data", "hic_test_file.hic"), 100_000)],
)
def test_hic_fetch_nnz_hic(file, resolution):
    f = hictkpy.hic.File(file, resolution)
    compare_nnz(f)
