import pathlib

import pytest

import hictkpy


def compare_nnz(f):
    assert f.fetch_nnz() == 890_384
    assert f.fetch_nnz("chr2R") == 31_900


def test_file_fetch_nnz_file(file: pathlib.Path = "test/data/cooler_test_file.cool", resolution: int = 100_000):
    f = hictkpy.File(file, resolution)
    compare_nnz(f)


def test_cooler_fetch_nnz_cooler(file: pathlib.Path = "test/data/cooler_test_file.cool"):
    f = hictkpy.cooler.File(file)
    compare_nnz(f)


def test_hic_fetch_nnz_hic(file: pathlib.Path = "test/data/hic_test_file.hic", resolution: int = 100_000):
    f = hictkpy.hic.File(file, resolution)
    compare_nnz(f)
