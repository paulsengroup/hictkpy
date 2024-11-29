# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

import pytest

import hictkpy

from .helpers import numpy_avail

testdir = pathlib.Path(__file__).resolve().parent

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (testdir / "data" / "cooler_test_file.mcool", 100_000),
        (testdir / "data" / "hic_test_file.hic", 100_000),
    ],
)


@pytest.mark.skipif(not numpy_avail(), reason="numpy is not available")
class TestClass:
    def test_repr(self, file, resolution):
        f = hictkpy.File(file, resolution)

        sel = f.fetch()
        assert str(sel) == "PixelSelector(ALL; COO; int32)"

        sel = f.fetch(join=True)
        assert str(sel) == "PixelSelector(ALL; BG2; int32)"

        sel = f.fetch(count_type="float")
        assert str(sel) == "PixelSelector(ALL; COO; float64)"

        sel = f.fetch("chr2L:0-10,000,000", "chr2L:5,000,000-20,000,000")
        assert str(sel) == "PixelSelector(chr2L:0-10000000; chr2L:5000000-20000000; COO; int32)"

    def test_coords(self, file, resolution):
        f = hictkpy.File(file, resolution)

        sel = f.fetch()
        assert sel.coord1() == ("ALL", 0, len(f.bins()))
        assert sel.coord1() == sel.coord2()

        sel = f.fetch("chr2L:0-10,000,000")
        assert sel.coord1() == ("chr2L", 0, 10_000_000)
        assert sel.coord1() == sel.coord2()

        sel = f.fetch("chr2L:0-10,000,000", "chr2L:5,000,000-20,000,000")
        assert sel.coord1() == ("chr2L", 0, 10_000_000)
        assert sel.coord2() == ("chr2L", 5_000_000, 20_000_000)
