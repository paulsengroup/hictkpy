import pathlib

import pytest

import hictkpy


def compare_sum(f):
    assert f.fetch_sum() == 119_208_613
    assert f.fetch_sum("chr2L") == 19_968_156


def test_file_fetch_sum_file(file: pathlib.Path = "test/data/cooler_test_file.cool", resolution: int = 100_000):
    f = hictkpy.File(file, resolution)
    compare_sum(f)


def test_cooler_fetch_sum_cooler(file: pathlib.Path = "test/data/cooler_test_file.cool"):
    f = hictkpy.cooler.File(file)
    compare_sum(f)


def test_hic_fetch_sum_hic(file: pathlib.Path = "test/data/hic_test_file.hic", resolution: int = 100_000):
    f = hictkpy.hic.File(file, resolution)
    compare_sum(f)
