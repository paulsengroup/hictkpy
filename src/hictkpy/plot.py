# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import functools
from typing import Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

import hictkpy


class HeatMap(object):
    def __init__(
        self,
        matrix: NDArray,
        coord1: Tuple[str, int, int] = None,
        coord2: Tuple[str, int, int] = None,
        ax: Optional[plt.Axes] = None,
        bins: Optional[hictkpy.BinTable] = None,
    ):
        if coord1 is None and coord2 is None and bins is None:
            raise ValueError("please provide at least one of the following parameters: coord1, coord2, bins")

        # TODO improve param validation

        self._fig = None
        self._ax = ax
        if self._ax is None:
            self._fig, self._ax = plt.subplots()

        self._matrix = matrix
        self._coord1 = coord1
        self._coord2 = coord1 if coord2 is None else coord2
        self._bins = bins

        if self._coord1 is None:
            assert self._bins is not None
            if len(self._bins) != matrix.shape[0]:
                raise ValueError("TODO")
            if matrix.shape[0] != matrix.shape[1]:
                raise ValueError("TODO")

        self._outdated = True

    def plot(self) -> Tuple[Optional[plt.Figure], plt.Axes]:
        if self._outdated:
            if self._coord1 is None:
                self._plot_gw()
            else:
                self._plot()

        return self._fig, self._ax

    @property
    def coord1(self) -> Union[str, Tuple[str, int, int]]:
        if self._coord1 is None:
            return "ALL"
        return self._coord1

    @property
    def coord2(self) -> Union[str, Tuple[str, int, int]]:
        if self._coord2 is None:
            return "ALL"
        return self._coord2

    def _plot(self):
        from matplotlib.colors import LogNorm

        _, start1, end1 = self._coord1
        _, start2, end2 = self._coord2

        extent = (start1, end1, end2, start2)

        self._ax.imshow(self._matrix, norm=LogNorm(), extent=extent)
        self._format_ticks()

        self._outdated = False

    def _plot_gw(self):
        assert self._bins is not None
        from matplotlib.colors import LogNorm

        genome_size = sum(self._bins.chromosomes(include_ALL=False).values())
        extent = (0, genome_size, genome_size, 0)

        self._ax.imshow(self._matrix, norm=LogNorm(), extent=extent)
        self._format_ticks_gw()

        self._outdated = False

    @staticmethod
    def _get_genomic_distance_formatter(suffix="bp"):
        from matplotlib.ticker import EngFormatter

        fmt = EngFormatter()
        # We do not want to use unit as we do not want to format 0 as "0 bp"
        fmt.ENG_PREFIXES = {
            0: "",
            3: f"k{suffix}",
            6: f"M{suffix}",
            9: f"G{suffix}",
            12: f"T{suffix}",
        }

        return fmt

    @staticmethod
    def _get_chromosome_formatter(chroms: Dict[str, int], resolution: int):
        from matplotlib.ticker import FuncFormatter

        chrom_sizes_bp = list(chroms.values())
        chrom_sizes_bins = [(pos + resolution - 1) // resolution for pos in chrom_sizes_bp]

        prefix_sum = np.cumsum(chrom_sizes_bins, dtype=np.int64)
        labels = list(chroms.keys())

        assert len(prefix_sum) == len(labels)

        def fx(x, pos) -> str:
            i = np.clip(np.searchsorted(prefix_sum, pos, side="right") - 1, 0, len(labels) - 1)
            return labels[i]

        return FuncFormatter(fx)

    @staticmethod
    @functools.cache
    def _get_custom_palettes() -> Dict[str, NDArray]:
        # Source: https://github.com/paulsengroup/StripePy/blob/main/src/stripepy/plot.py
        return {
            "fall": np.array(
                (
                    (255, 255, 255),
                    (255, 255, 204),
                    (255, 237, 160),
                    (254, 217, 118),
                    (254, 178, 76),
                    (253, 141, 60),
                    (252, 78, 42),
                    (227, 26, 28),
                    (189, 0, 38),
                    (128, 0, 38),
                    (0, 0, 0),
                )
            )
            / 255,
            "fruit_punch": np.array(
                (
                    (255, 255, 255),
                    (255, 204, 204),
                    (255, 153, 153),
                    (255, 102, 102),
                    (255, 50, 50),
                    (255, 0, 0),
                )
            )
            / 255,
            "blues": np.array(
                (
                    (255, 255, 255),
                    (180, 204, 225),
                    (116, 169, 207),
                    (54, 144, 192),
                    (5, 112, 176),
                    (4, 87, 135),
                    (3, 65, 100),
                    (2, 40, 66),
                    (1, 20, 30),
                    (0, 0, 0),
                )
            )
            / 255,
            "acidblues": np.array(
                (
                    (255, 255, 255),
                    (162, 192, 222),
                    (140, 137, 187),
                    (140, 87, 167),
                    (140, 45, 143),
                    (120, 20, 120),
                    (90, 15, 90),
                    (60, 10, 60),
                    (30, 5, 30),
                    (0, 0, 0),
                )
            )
            / 255,
            "nmeth": np.array(
                (
                    (236, 250, 255),
                    (148, 189, 217),
                    (118, 169, 68),
                    (131, 111, 43),
                    (122, 47, 25),
                    (41, 0, 20),
                )
            )
            / 255,
        }

    @staticmethod
    def _list_to_colormap(color_list, name=None):
        import matplotlib as mpl

        color_list = np.array(color_list)
        if color_list.min() < 0:
            raise ValueError("Colors should be 0 to 1, or 0 to 255")
        if color_list.max() > 1.0:
            if color_list.max() > 255:
                raise ValueError("Colors should be 0 to 1 or 0 to 255")
            else:
                color_list = color_list / 255.0
        return mpl.colors.LinearSegmentedColormap.from_list(name, color_list, 256)

    def _format_ticks(self, xaxis: bool = True, yaxis: bool = True, rotation: int = 30):
        """
        Function taken from https://cooltools.readthedocs.io/en/latest/notebooks/viz.html
        """
        from matplotlib.ticker import ScalarFormatter

        if xaxis:
            self._ax.xaxis.set_major_formatter(HeatMap._get_genomic_distance_formatter())
        else:
            self._ax.xaxis.set_major_formatter(ScalarFormatter())
        self._ax.xaxis.tick_bottom()

        if yaxis:
            self._ax.yaxis.set_major_formatter(HeatMap._get_genomic_distance_formatter())
        else:
            self._ax.yaxis.set_major_formatter(ScalarFormatter())

        if rotation != 0:
            self._ax.tick_params(axis="x", rotation=rotation)

    def _format_ticks_gw(self, xaxis: bool = True, yaxis: bool = True):
        assert self._bins is not None

        from matplotlib.ticker import ScalarFormatter

        if xaxis or yaxis:
            chroms = self._bins.chromosomes(include_ALL=False)
            prefix_sum = np.cumsum([0] + list(chroms.values()), dtype=np.int64)

            ticks = []
            for x1, x2 in zip(prefix_sum[:-1], prefix_sum[1:]):
                ticks.append((x1 + x2) / 2)

            chrom_ticks = {"ticks": ticks, "labels": list(chroms.keys())}
        else:
            chrom_ticks = {}

        if xaxis:
            self._ax.set_xticks(**chrom_ticks)
        else:
            self._ax.xaxis.set_major_formatter(ScalarFormatter())
        self._ax.xaxis.tick_bottom()

        if yaxis:
            self._ax.set_yticks(**chrom_ticks)
        else:
            self._ax.yaxis.set_major_formatter(ScalarFormatter())

        # TODO fix overlapping labels
