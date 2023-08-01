# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os

import pytest

import hictkpy

testdir = os.path.dirname(os.path.abspath(__file__))


def compare_sum(f):
    assert f.fetch_sum() == 119_208_613
    assert f.fetch_sum("chr2L") == 19_968_156


@pytest.mark.parametrize(
    "file,resolution",
    [(os.path.join(testdir, "data", "cooler_test_file.cool"), 100_000)],
)
def test_file_fetch_sum_file(file, resolution):
    f = hictkpy.File(file, resolution)
    compare_sum(f)


@pytest.mark.parametrize(
    "file",
    [(os.path.join(testdir, "data", "cooler_test_file.cool"))],
)
def test_cooler_fetch_sum_cooler(file):
    f = hictkpy.cooler.File(file)
    compare_sum(f)


@pytest.mark.parametrize(
    "file,resolution",
    [(os.path.join(testdir, "data", "hic_test_file.hic"), 100_000)],
)
def test_hic_fetch_sum_hic(file, resolution):
    f = hictkpy.hic.File(file, resolution)
    compare_sum(f)
