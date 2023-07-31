import os

import pytest

import hictkpy

testdir = os.path.dirname(os.path.abspath(__file__))


def compare_shared_attributes(f):
    assert f.bin_size() == 100_000
    assert f.nbins() == 1380

    assert "chr2L" in f.chromosomes()
    assert len(f.bins()) == 1380
    assert len(f.chromosomes()) == 8


@pytest.mark.parametrize(
    "file,resolution",
    [(os.path.join(testdir, "data", "cooler_test_file.cool"), 100_000)],
)
def test_file_fetch_accessors_file(file, resolution):
    f = hictkpy.File(file, resolution)
    compare_shared_attributes(f)


@pytest.mark.parametrize(
    "file",
    [(os.path.join(testdir, "data", "cooler_test_file.cool"))],
)
def test_cooler_fetch_accessors_cooler(file):
    f = hictkpy.cooler.File(file)
    compare_shared_attributes(f)


@pytest.mark.parametrize(
    "file,resolution",
    [(os.path.join(testdir, "data", "hic_test_file.hic"), 100_000)],
)
def test_hic_fetch_accessors_hic(file, resolution):
    f = hictkpy.hic.File(file, resolution)
    compare_shared_attributes(f)
