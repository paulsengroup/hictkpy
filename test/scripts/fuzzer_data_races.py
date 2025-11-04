#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import json
import logging
import multiprocessing as mp
import pathlib
import random
import sys
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Dict

import hictkpy


def make_cli() -> argparse.ArgumentParser:
    def nproc() -> int:
        return mp.cpu_count()

    def positive_int(arg) -> int:
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    def non_negative_int(arg) -> int:
        if (n := int(arg)) >= 0:
            return n

        raise ValueError("Not a non-negative integer")

    def valid_nproc(arg) -> int:
        if 1 <= (n := int(arg)) <= nproc():
            return n

        raise ValueError(f"Not a number between 1 and {nproc()}")

    cli = argparse.ArgumentParser(description="Fuzzer for hictkpy.")

    cli.add_argument(
        "hic-file",
        type=pathlib.Path,
        help="Path to a .hic file.",
    )
    cli.add_argument(
        "cooler-file",
        type=pathlib.Path,
        help="Path to an .[m]cool file.",
    )
    cli.add_argument(
        "--resolution",
        type=non_negative_int,
        required=True,
        help="Matrix resolution.",
    )
    cli.add_argument(
        "--duration",
        type=positive_int,
        default=60,
        help="Duration in seconds.",
    )
    cli.add_argument(
        "--seed",
        type=int,
        help="Seed used for PRNG.",
    )
    cli.add_argument(
        "--nthreads",
        type=valid_nproc,
        default=nproc(),
        help="Number of threads to use for fuzzing.",
    )
    cli.add_argument(
        "--tmpdir",
        type=pathlib.Path,
        default=None,
        help="Path pointing to an existing folder where to store temporary files.",
    )
    return cli


def call_to_arrow(f: hictkpy.File) -> int:
    return len(f.fetch().to_arrow())


def call_to_pandas(f: hictkpy.File) -> int:
    return len(f.fetch().to_pandas())


def call_to_df(f: hictkpy.File) -> int:
    return len(f.fetch().to_df())


def call_to_numpy(f: hictkpy.File) -> int | float:
    return f.fetch().to_numpy().sum()


def call_to_csr(f: hictkpy.File) -> int | float:
    return f.fetch().to_csr().sum()


def call_to_coo(f: hictkpy.File) -> int | float:
    return f.fetch().to_coo().sum()


def call_describe(f: hictkpy.File) -> int:
    return f.fetch().describe()["nnz"]


def call_nnz(f: hictkpy.File) -> int:
    return f.fetch().nnz()


def call_sum(f: hictkpy.File) -> int | float:
    return f.fetch().sum()


def call_min(f: hictkpy.File) -> int | float:
    return f.fetch().min()


def call_max(f: hictkpy.File) -> int | float:
    return f.fetch().max()


def call_variance(f: hictkpy.File) -> float:
    return f.fetch().variance()


def call_skewness(f: hictkpy.File) -> float:
    return f.fetch().skewness()


def call_kurtosis(f: hictkpy.File) -> float:
    return f.fetch().kurtosis()


@functools.cache
def vtable() -> Dict[str, Callable]:
    names = [
        "call_to_arrow",
        # "call_to_pandas",
        # "call_to_df",
        "call_to_numpy",
        "call_to_csr",
        "call_to_coo",
        "call_describe",
        "call_nnz",
        "call_sum",
        "call_min",
        "call_max",
        "call_variance",
        "call_skewness",
        "call_kurtosis",
    ]

    return {name: globals()[name] for name in names}


def runme(
    path: pathlib.Path,
    resolution: int,
    end_time: float,
    seed: int | None,
) -> Dict[str, int]:
    if seed is None:
        seed = threading.get_ident()
    prng = random.Random(seed)
    counters = {name: 0 for name in vtable()}
    file_type = path.suffix.lstrip(".")
    with hictkpy.File(path, resolution) as f:
        i = 1
        while time.time() < end_time:
            name, fx = prng.choice(list(vtable().items()))
            logging.info(
                "[tid=%d] [iter=%d]: %s[%s]",
                threading.get_ident(),
                i,
                name,
                file_type,
            )
            fx(f)
            counters[name] += 1
            i += 1

    return counters


def run_fuzzer(
    hic_file: pathlib.Path,
    cooler_file: pathlib.Path,
    resolution: int,
    nthreads: int,
    duration: float,
    seed: int | None,
) -> Counter:
    end_time = time.time() + duration

    results = Counter()
    with ThreadPoolExecutor(nthreads) as tpool:
        tasks = []
        for i in range(nthreads):
            if seed is not None:
                seed += 1
            if i % 2 == 0:
                path = hic_file
            else:
                path = cooler_file

            tasks.append(
                tpool.submit(
                    runme,
                    path=path,
                    resolution=resolution,
                    end_time=end_time,
                    seed=seed,
                )
            )

        for t in tasks:
            results += Counter(t.result())

        return results


def main() -> int:
    args = vars(make_cli().parse_args())

    hic_file = args["hic-file"]
    cooler_file = args["cooler-file"]

    assert hictkpy.is_hic(hic_file)
    assert hictkpy.is_cooler(cooler_file) or hictkpy.is_mcool_file(cooler_file)

    results = run_fuzzer(
        hic_file,
        cooler_file,
        resolution=args["resolution"],
        nthreads=args["nthreads"],
        duration=args["duration"],
        seed=args["seed"],
    )
    results["total"] = sum(results.values())

    json.dump(
        results,
        sys.stdout,
        indent=2,
        sort_keys=True,
    )
    sys.stdout.write("\n")

    return 0


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)
    hictkpy.logging.setLevel(logging.WARN)


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
