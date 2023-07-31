import pathlib

import numpy as np
import pytest

import hictkpy


def compare_pixels(f):
    ### GW
    df = f.fetch()
    assert df["count"].sum() == 119_208_613
    assert len(df) == 890_384

    ### CIS
    df = f.fetch("chr2R:10,000,000-15,000,000")
    assert df["count"].sum() == 4_519_080
    assert len(df.columns) == 3

    df = f.fetch("chr2R:10,000,000-15,000,000", join=True)
    assert df["count"].sum() == 4_519_080
    assert len(df.columns) == 7

    df = f.fetch("chr2R:10,000,000-15,000,000", count_type="int")
    assert df["count"].dtype == np.int32

    df = f.fetch("chr2R:10,000,000-15,000,000", count_type="float")
    assert df["count"].dtype == np.float64

    df = f.fetch("chr2R\t10000000\t15000000", query_type="BED")
    assert len(df) == 1275

    ### TRANS
    df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000")
    assert df["count"].sum() == 83_604
    assert len(df.columns) == 3

    df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", join=True)
    assert df["count"].sum() == 83_604
    assert len(df.columns) == 7

    df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="int")
    assert df["count"].dtype == np.int32

    df = f.fetch("chr2R:10,000,000-15,000,000", "chrX:0-10,000,000", count_type="float")
    assert df["count"].dtype == np.float64

    df = f.fetch("chr2R\t10000000\t15000000", "chrX\t0\t10000000", query_type="BED")
    assert len(df) == 4995


def test_file_fetch_file(file: pathlib.Path = "test/data/cooler_test_file.cool", resolution: int = 100_000):
    f = hictkpy.File(file, resolution)
    compare_pixels(f)

    df = f.fetch("chr2R:10,000,000-15,000,000", normalization="weight")
    assert np.isclose(59.349524704033215, df["count"].sum())


def test_cooler_fetch_cooler(file: pathlib.Path = "test/data/cooler_test_file.cool"):
    f = hictkpy.cooler.File(file)
    compare_pixels(f)

    df = f.fetch("chr2R:10,000,000-15,000,000", normalization="weight")
    assert np.isclose(59.349524704033215, df["count"].sum())


def test_hic_fetch_hic(file: pathlib.Path = "test/data/hic_test_file.hic", resolution: int = 100_000):
    f = hictkpy.hic.File(file, resolution)
    compare_pixels(f)

    df = f.fetch("chr2R:10,000,000-15,000,000", normalization="ICE")
    assert np.isclose(59.349524704033215, df["count"].sum())
