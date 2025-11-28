# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pytest

import hictkpy

from .helpers import get_test_dir, numpy_avail

pytestmark = pytest.mark.parametrize(
    "file,resolution",
    [
        (get_test_dir() / "data" / "cooler_test_file.mcool", 100_000),
        (get_test_dir() / "data" / "hic_test_file.hic", 100_000),
    ],
)


@pytest.mark.skipif(not numpy_avail(), reason="numpy is not available")
class TestPixelSelectorAccessors:
    def test_repr(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            sel = f.fetch()
            assert str(sel) == "PixelSelector(ALL; COO; int32)"

            sel = f.fetch(join=True)
            assert str(sel) == "PixelSelector(ALL; BG2; int32)"

            sel = f.fetch(count_type="float")
            assert str(sel) == "PixelSelector(ALL; COO; float64)"

            sel = f.fetch("chr2L:0-10,000,000", "chr2L:5,000,000-20,000,000")
            assert str(sel) == "PixelSelector(chr2L:0-10000000; chr2L:5000000-20000000; COO; int32)"

    def test_coords(self, file, resolution):
        with hictkpy.File(file, resolution) as f:
            sel = f.fetch()
            assert sel.coord1() is None
            assert sel.coord2() is None

            sel = f.fetch("chr2L:0-10,000,000")
            assert sel.coord1() == ("chr2L", 0, 10_000_000)
            assert sel.coord1() == sel.coord2()

            sel = f.fetch("chr2L:0-10,000,000", "chr2L:5,000,000-20,000,000")
            assert sel.coord1() == ("chr2L", 0, 10_000_000)
            assert sel.coord2() == ("chr2L", 5_000_000, 20_000_000)
