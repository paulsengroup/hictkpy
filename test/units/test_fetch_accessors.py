import pathlib

import pytest

import hictkpy


def compare_shared_attributes(f):
    assert f.bin_size() == 100_000
    assert f.nbins() == 1380

    assert "chr2L" in f.chromosomes()
    assert len(f.bins()) == 1380
    assert len(f.chromosomes()) == 8


def test_file_fetch_accessors_file(file: pathlib.Path = "test/data/cooler_test_file.cool", resolution: int = 100_000):
    f = hictkpy.File(file, resolution)
    compare_shared_attributes(f)


def test_cooler_fetch_accessors_cooler(file: pathlib.Path = "test/data/cooler_test_file.cool"):
    f = hictkpy.cooler.File(file)
    compare_shared_attributes(f)


def test_hic_fetch_accessors_hic(file: pathlib.Path = "test/data/hic_test_file.hic", resolution: int = 100_000):
    f = hictkpy.hic.File(file, resolution)
    compare_shared_attributes(f)
