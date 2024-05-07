#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import timeit
from typing import Dict, Tuple

import numpy as np
import numpy.typing as npt
import bioframe as bf

import hictkpy
import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "--bin-size",
        type=int,
        default=1000,
    )

    cli.add_argument(
        "--batch-size",
        type=int,
        default=20_000_000,
    )

    cli.add_argument(
        "--seed",
        type=int,
        default=123456789,
    )

    cli.add_argument(
        "--min-time",
        type=float,
        default=5,
    )

    return cli


def setup(chroms: Dict[str, int], bin_size: int, batch_size: int) -> Tuple[hictkpy.BinTable, npt.NDArray]:
    assert batch_size > 0
    bins = hictkpy.BinTable(chroms, bin_size)
    bin_ids = np.random.randint(0, len(bins), batch_size)

    return bins, bin_ids


def benchmark_hictkpy_int(bins: hictkpy.BinTable, chrom: str, start: int, end: int, iters: int):
    return timeit.timeit(lambda: bins.get_id(chrom, start), number=iters)


def benchmark_hictkpy_vec(
    bins: hictkpy.BinTable, chroms: npt.NDArray, starts: npt.NDArray, ends: npt.NDArray, iters: int
):
    return timeit.timeit(lambda: bins.get_ids(chroms, starts), number=iters)


def benchmark_pandas_int(bins: pd.DataFrame, chrom: str, start: int, end: int, iters: int):
    return timeit.timeit(
        lambda: bins[(bins["chrom"] == chrom) & (bins["start"] <= start) & (bins["end"] > start)].index[0], number=iters
    )


def benchmark_pandas_vec(bins: pd.DataFrame, chroms: npt.NDArray, starts: npt.NDArray, ends: npt.NDArray, iters: int):
    df = pd.DataFrame({"chrom": chroms, "start": starts, "end": ends})
    # fmt: off
    return timeit.timeit(
        lambda: df.apply(
            lambda row: bins[
                (bins["chrom"] == row["chrom"]) &
                (bins["start"] <= row["start"]) &
                (bins["end"] > row["start"])
            ].index[0],
            number=iters,
        )
    )
    # fmt: on


def benchmark(
    bins: hictkpy.BinTable, chroms: npt.NDArray, starts: npt.NDArray, ends: npt.NDArray, mode: str, min_time: float
) -> Tuple[int, float]:
    if len(chroms) == 1:
        fx = benchmark_hictkpy_int if mode == "hictkpy" else benchmark_pandas_int
        chroms = chroms[0]
        starts = starts[0]
        ends = ends[0]
    else:
        fx = benchmark_hictkpy_vec if mode == "hictkpy" else benchmark_pandas_vec

    n = 1
    t = fx(bins, chroms, starts, ends, n)

    if t < min_time:
        n = int(0.5 * (min_time / t))

    while t < min_time:
        n *= 2
        t = fx(bins, chroms, starts, ends, n)

    return n, t / n


def make_time_human_readable(time: float) -> Tuple[float, str]:
    suffix = "s"
    if time < 1.0e-6:
        time *= 1.0e9  # sec -> ns
        suffix = "ns"
    elif time < 1.0e-3:  # sec -> us
        time *= 1.0e6
        suffix = "us"
    elif time < 1:  # sec -> ns
        time *= 1.0e3
        suffix = "ms"

    return time, suffix


def main():
    args = vars(make_cli().parse_args())

    np.random.seed(args["seed"])

    bins, bin_ids = setup(hg38, args["bin_size"], args["batch_size"])

    coords = [bins.get(id) for id in bin_ids]
    chroms = np.array([c.chrom for c in coords])
    starts = np.array([c.start for c in coords])
    ends = np.array([c.end for c in coords])

    iters_htk, time_htk = benchmark(bins, chroms, starts, ends, min_time=args["min_time"], mode="hictkpy")
    iters_pd, time_pd = benchmark(bins.to_df(), chroms, starts, ends, min_time=args["min_time"], mode="pandas")

    time_htk, suffix_htk = make_time_human_readable(time_htk)
    time_pd, suffix_pd = make_time_human_readable(time_pd)

    print(f"hictkpy: {time_htk:.2f} {suffix_htk} {iters_htk} runs")
    print(f"pandas: {time_pd:.2f} {suffix_pd} {iters_pd} runs")


if __name__ == "__main__":
    hg38 = {
        "chr1": 248956422,
        "chr2": 242193529,
        "chr3": 198295559,
        "chr4": 190214555,
        "chr5": 181538259,
        "chr6": 170805979,
        "chr7": 159345973,
        "chr8": 145138636,
        "chr9": 138394717,
        "chr10": 133797422,
        "chr11": 135086622,
        "chr12": 133275309,
        "chr13": 114364328,
        "chr14": 107043718,
        "chr15": 101991189,
        "chr16": 90338345,
        "chr17": 83257441,
        "chr18": 80373285,
        "chr19": 58617616,
        "chr20": 64444167,
        "chr21": 46709983,
        "chr22": 50818468,
        "chrX": 156040895,
        "chrY": 57227415,
    }

    main()
