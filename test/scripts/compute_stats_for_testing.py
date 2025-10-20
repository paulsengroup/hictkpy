#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import warnings
from typing import Any, Dict, List, Sequence

import numpy as np
import scipy.stats as ss

import hictkpy


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "uri",
        type=pathlib.Path,
        help="Path to a .cool, .mcool or hic file.",
    )
    cli.add_argument(
        "resolution",
        type=int,
        help="Resolution to be used.",
    )
    cli.add_argument(
        "--metrics",
        type=str,
        choices={"nnz", "sum", "min", "max", "mean", "variance", "skewness", "kurtosis"},
        help="One or more stats to be computed.",
    )
    cli.add_argument("--range", nargs="+", type=str, help="Coordinates of the genomic regions to be processed")
    cli.add_argument("--range2", nargs="+", type=str, help="Coordinates of the genomic regions to be processed")
    cli.add_argument(
        "--normalization", nargs="+", type=str, help="Balance interactions using the given normalization method."
    )

    return cli


def format_command(
    range1: str | None,
    range2: str | None,
    normalization: str | None,
    keep_nans: bool,
    keep_infs: bool,
    metric: str,
    value: int | float,
) -> str:
    if np.isnan(value):
        prefix = "assert isnan(f.fetch("
        suffix = ")"
    elif np.isinf(value):
        prefix = "assert isinf(f.fetch("
        suffix = ")"
    elif isinstance(value, float):
        prefix = "assert isclose(f.fetch("
        suffix = f", {value:.16g})"
    else:
        prefix = "assert f.fetch("
        suffix = f" == {value:_}"

    args = []
    if range1:
        args.append(f'"{range1}"')
    if range2 and range1 != range2:
        args.append(f'"{range2}"')
    if normalization is not None and normalization != "NONE":
        args.append(f'normalization="{normalization}"')

    msg = prefix + ", ".join(args)
    msg += f").{metric}("
    args = []
    if keep_nans:
        args.append("keep_nans=True")
    if keep_infs:
        args.append("keep_infs=True")
    msg += f"{", ".join(args)})"
    msg += suffix

    return msg


def format_stats(
    range1: str,
    range2: str,
    normalization: str,
    stats: Dict[str, Any],
    keep_nans: bool = False,
    keep_infs: bool = False,
) -> List[str]:
    messages = []
    for metric, value in stats.items():
        messages.append(format_command(range1, range2, normalization, keep_nans, keep_infs, metric, value))
    return messages


def compute_stats(counts: Sequence[int] | Sequence[float], metrics: List[str]) -> Dict[str, Any]:
    with warnings.catch_warnings(action="ignore"):
        stats = ss.describe(counts, nan_policy="propagate")
        stats = {
            "nnz": stats.nobs,
            "sum": sum(counts),
            "min": stats.minmax[0],
            "max": stats.minmax[1],
            "mean": stats.mean,
            "variance": stats.variance,
            "skewness": stats.skewness,
            "kurtosis": stats.kurtosis,
        }

    return {metric: value for metric, value in stats.items() if metric in metrics}


def process_param_combination(f: hictkpy.File, range: str, range2: str, normalization: str, metrics: List[str]):
    counts = f.fetch(range, range2, normalization).to_df()["count"]

    raw = compute_stats(counts, metrics)
    keep_nans = compute_stats(counts[~np.isinf(counts)], metrics)
    keep_infs = compute_stats(counts[~np.isnan(counts)], metrics)
    finite = compute_stats(counts[np.isfinite(counts)], metrics)

    messages = format_stats(range, range2, normalization, finite)
    messages.extend(format_stats(range, range2, normalization, keep_nans, keep_nans=True))
    messages.extend(format_stats(range, range2, normalization, keep_infs, keep_infs=True))
    messages.extend(format_stats(range, range2, normalization, raw, keep_nans=True, keep_infs=True))

    print("\n".join(sorted(messages)))


def main():
    args = vars(make_cli().parse_args())

    if args["metrics"] is None:
        args["metrics"] = ["nnz", "sum", "min", "max", "mean", "variance", "skewness", "kurtosis"]

    if args["range"] is None:
        args["range"] = [""]

    if args["range2"] is None:
        args["range2"] = [""]

    if len(args["range"]) != len(args["range2"]):
        raise RuntimeError("please provide the same number of arguments for --range and --range2")

    if args["normalization"] is None:
        args["normalization"] = ["NONE"]

    with hictkpy.File(args["uri"], args["resolution"]) as f:
        for norm in args["normalization"]:
            for range1, range2 in zip(args["range"], args["range2"]):
                process_param_combination(f, range1, range2, norm, args["metrics"])
                print()


if __name__ == "__main__":
    main()
