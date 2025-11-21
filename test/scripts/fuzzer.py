#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import logging
import multiprocessing as mp
import pathlib
import random
import sys
import time
import warnings
from typing import Any, Dict, List, Tuple

import cooler
import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.sparse as ss
from scipy.stats import describe

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

    def valid_fraction(arg) -> float:
        if (n := float(arg)) >= 0 and n <= 1:
            return n

        raise ValueError("Not a number between 0 and 1")

    def valid_nproc(arg) -> int:
        if 1 <= (n := int(arg)) <= nproc():
            return n

        raise ValueError(f"Not a number between 1 and {nproc()}")

    cli = argparse.ArgumentParser(description="Fuzzer for hictkpy.")

    cli.add_argument(
        "test-uri",
        type=pathlib.Path,
        help="Path to a .cool, .mcool or hic file to be used as test file.",
    )
    cli.add_argument(
        "reference-uri",
        type=pathlib.Path,
        help="Path to a .cool or .mcool file to be used as reference.",
    )
    cli.add_argument(
        "--resolution",
        type=non_negative_int,
        help="Matrix resolution.\n"
        "Required when one or both of test-uri and reference-uri are multi-resolution files.",
    )
    cli.add_argument(
        "--1d-to-2d-query-ratio",
        type=valid_fraction,
        default=0.33,
        help="Ratio of 1D to 2D queries. Use 0 or 1 to only test 1D or 2D queries.",
    )
    cli.add_argument(
        "--duration",
        type=positive_int,
        default=60,
        help="Duration in seconds.",
    )
    cli.add_argument(
        "--format",
        type=str,
        choices={"df", "numpy", "csr", "describe"},
        default="df",
        help="Format used to fetch pixels.",
    )
    cli.add_argument(
        "--query-length-avg",
        type=positive_int,
        default=1_000_000,
        help="Average query size.",
    )
    cli.add_argument(
        "--query-length-std",
        type=non_negative_int,
        default=250_000,
        help="Standard deviation for query size.",
    )
    cli.add_argument(
        "--normalization",
        type=str,
        default="NONE",
        help="Name of the dataset to use for balancing.",
    )
    cli.add_argument(
        "--join",
        action="store_true",
        default=False,
        help="Fetch pixels in BG2 format.\n" "Ignored when --format is not df.",
    )
    cli.add_argument(
        "--seed",
        type=int,
        help="Seed used for PRNG.",
    )
    cli.add_argument(
        "--excluded-chroms",
        nargs="+",
        type=str,
        help="One or more chromosomes to be excluded.",
    )
    cli.add_argument(
        "--nproc",
        type=valid_nproc,
        default=nproc(),
        help="Number of test processes to run in parallel.",
    )
    return cli


def postproc_df(df: pd.DataFrame) -> pd.DataFrame:
    if "balanced" in df:
        df["count"] = df["balanced"]

    if "bin1_id" in df:
        return df[["bin1_id", "bin2_id", "count"]]

    df["chrom1"] = df["chrom1"].astype(str)
    df["chrom2"] = df["chrom2"].astype(str)
    return df.set_index(["chrom1", "start1", "end1", "chrom2", "start2", "end2"])[["count"]]


def cooler_dump(selector, query1: str, query2: str) -> pd.DataFrame | ss.coo_matrix | npt.NDArray:
    logging.debug("[cooler] running query for %s, %s...", query1, query2)
    data = selector.fetch(query1, query2)
    if isinstance(data, pd.DataFrame):
        return postproc_df(data)

    return data


def hictk_dump(
    file,
    query1: str,
    query2: str,
    normalization: str = "NONE",
    query_type: str = "df",
) -> pd.DataFrame | npt.NDArray | ss.csr_matrix | List[Any]:
    logging.debug("[hictkpy] running query for %s, %s...", query1, query2)
    if query_type == "df":
        return postproc_df(file.fetch(query1, query2, normalization, join=True).to_df())
    if query_type == "csr":
        return file.fetch(query1, query2, normalization).to_csr(query_span="full")
    if query_type == "numpy":
        return file.fetch(query1, query2, normalization).to_numpy(query_span="full")

    raise NotImplementedError


def read_chrom_sizes_cooler(path_to_cooler_file: str) -> Dict[str, int]:
    return cooler.Cooler(path_to_cooler_file).chromsizes.to_dict()


def generate_query_1d(chroms, weights: npt.NDArray, mean_length: float, stddev_length: float) -> str:
    chrom_name, chrom_size = random.choices(chroms, weights=weights, k=1)[0]

    query_length = max(2.0, random.gauss(mu=mean_length, sigma=stddev_length))

    center_pos = random.randint(0, chrom_size)
    start_pos = max(0.0, center_pos - (query_length / 2))
    end_pos = min(chrom_size, start_pos + query_length)

    return f"{chrom_name}:{start_pos:.0f}-{end_pos:.0f}"


def generate_query_2d(
    chroms,
    weights: npt.NDArray,
    ranks: Dict[str, int],
    mean_length: float,
    stddev_length: float,
) -> Tuple[str, str]:
    q1 = generate_query_1d(chroms, weights, mean_length, stddev_length)
    q2 = generate_query_1d(chroms, weights, mean_length, stddev_length)

    chrom1, _, coord1 = q1.partition(":")
    chrom2, _, coord2 = q2.partition(":")

    if ranks[chrom1] > ranks[chrom2]:
        q1, q2 = q2, q1

    if chrom1 == chrom2:
        start1, _, _ = coord1.partition("-")
        start2, _, _ = coord2.partition("-")
        if int(start1) > int(start2):
            q1, q2 = q2, q1

    return q1, q2


def find_differences_df(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    df = df1.merge(
        df2,
        how="outer",
        left_index=True,
        right_index=True,
        suffixes=("1", "2"),
    )
    # We're mapping False to None so that we can more easily drop identical rows with dropna()
    df["count_close_enough"] = pd.Series(np.isclose(df["count1"], df["count2"])).map({False: None})

    # We're dropping the counts to avoid incorrectly flagging rows with nan as counts
    return df.drop(columns=["count1", "count2"]).dropna()


def compare_dfs(worker_id: int, q1: str, q2: str, expected: pd.DataFrame, found: pd.DataFrame) -> bool:
    if len(expected) != len(found):
        logging.warning(
            "[%d] %s, %s: FAIL! Expected %d nnz, found %d!",
            worker_id,
            q1,
            q2,
            len(expected),
            len(found),
        )
        return False

    if len(expected) != 0:
        diff = find_differences_df(expected, found)
        if len(diff) != 0:
            logging.warning(
                "[%d] %s, %s (%d nnz): FAIL! Found %d differences!",
                worker_id,
                q1,
                q2,
                len(expected),
                len(diff),
            )
            return False

    logging.debug("[%d] %s, %s (%d nnz): OK!", worker_id, q1, q2, len(expected))
    return True


def compare_numpy(
    worker_id: int, q1: str, q2: str, expected: npt.NDArray, found: npt.NDArray, rtol: float = 1.0e-5
) -> bool:
    nnz = np.sum(expected != 0)
    if expected.shape != found.shape:
        logging.warning(
            "[%d] %s, %s (%d nnz): FAIL! Numpy matrices have different shapes! Expected %s, found %s",
            worker_id,
            q1,
            q2,
            nnz,
            expected.shape,
            found.shape,
        )
        return False

    num_differences = (~np.isclose(expected, found, rtol=rtol, equal_nan=True)).sum()
    if num_differences != 0:
        logging.warning(
            "[%d] %s, %s (%d nnz): FAIL! Found %d differences!",
            worker_id,
            q1,
            q2,
            nnz,
            num_differences,
        )
        return False

    logging.debug("[%d] %s, %s (%d nnz): OK!", worker_id, q1, q2, np.sum(expected != 0))
    return True


def compare_csr(
    worker_id: int, q1: str, q2: str, expected: ss.csr_matrix, found: ss.csr_matrix, rtol: float = 1.0e-5
) -> bool:
    if expected.shape != found.shape:
        logging.warning(
            "[%d] %s, %s (%d nnz): FAIL! CSR matrices have different shapes! Expected %s, found %s",
            worker_id,
            q1,
            q2,
            expected.nnz,
            expected.shape,
            found.shape,
        )
        return False

    if expected.nnz != found.nnz:
        logging.warning(
            "[%d] %s, %s (%d nnz): FAIL! CSR matrices have different nnz! Expected %d, found %d",
            worker_id,
            q1,
            q2,
            expected.nnz,
            expected.nnz,
            found.nnz,
        )
        return False

    expected.sort_indices()
    found.sort_indices()

    num_differences = (expected.indices != found.indices).sum()
    if num_differences != 0:
        logging.warning(
            "[%d] %s, %s (%d nnz): FAIL! Found %d differences in CSR matrix.indices!",
            worker_id,
            q1,
            q2,
            expected.nnz,
            num_differences,
        )
        return False

    num_differences = (expected.indptr != found.indptr).sum()
    if num_differences != 0:
        logging.warning(
            "[%d] %s, %s (%d nnz): FAIL! Found %d differences in CSR matrix.indptr!",
            worker_id,
            q1,
            q2,
            expected.nnz,
            num_differences,
        )
        return False

    num_differences = (~np.isclose(expected.data, found.data, rtol=rtol, equal_nan=True)).sum()
    if num_differences != 0:
        logging.warning(
            "[%d] %s, %s (%d nnz): FAIL! Found %d differences!",
            worker_id,
            q1,
            q2,
            expected.nnz,
            num_differences,
        )
        return False

    logging.debug("[%d] %s, %s (%d nnz): OK!", worker_id, q1, q2, expected.nnz)
    return True


def compare_query_results(worker_id: int, q1: str, q2: str, expected, found) -> int:
    if isinstance(found, pd.DataFrame):
        assert isinstance(expected, type(found))
        return int(not compare_dfs(worker_id, q1, q2, expected, found))
    if isinstance(found, (np.ndarray, np.generic)):
        assert isinstance(expected, type(found))
        return int(not compare_numpy(worker_id, q1, q2, expected, found))
    if isinstance(found, ss.csr_matrix):
        return int(not compare_csr(worker_id, q1, q2, expected.tocsr(), found))

    raise NotImplementedError


def compute_stats(df: pd.DataFrame, keep_nans: bool, keep_infs: bool) -> Dict:
    if not keep_nans:
        df = df[~np.isnan(df["count"])]
    if not keep_infs:
        df = df[~np.isinf(df["count"])]

    if len(df) == 0:
        return {
            "nnz": 0,
            "sum": 0.0 if np.issubdtype(df["count"].dtype, np.floating) else 0,
            "min": None,
            "max": None,
            "mean": None,
            "variance": None,
            "skewness": None,
            "kurtosis": None,
        }

    with warnings.catch_warnings(action="ignore"):
        stats = describe(df["count"], nan_policy="propagate")
        return {
            "nnz": stats.nobs,
            "sum": df["count"].sum(skipna=False),
            "min": stats.minmax[0],
            "max": stats.minmax[1],
            "mean": stats.mean,
            "variance": stats.variance if len(df) > 1 else None,
            "skewness": stats.skewness if len(df) > 1 else None,
            "kurtosis": stats.kurtosis if len(df) > 1 else None,
        }


def compare_metric(
    worker_id: int,
    q1: str,
    q2: str,
    metric: str,
    expected,
    found,
    keep_nans: bool,
    keep_infs: bool,
    rtol: float = 1.0e-5,
    atol: float = 1.0e-8,
) -> bool:
    do_numeric_comparison = expected is not None and found is not None
    if do_numeric_comparison:
        if not np.isclose(expected, found, rtol=rtol, atol=atol, equal_nan=True):
            logging.warning(
                "[%d] %s, %s (%s; keep_nans=%s; keep_infs=%s): FAIL! Expected %.16g, found %.16g",
                worker_id,
                q1,
                q2,
                metric,
                keep_nans,
                keep_infs,
                expected,
                found,
            )
            return False
        return True

    if expected is None and found is None:
        return True

    logging.warning(
        "[%d] %s, %s (%s; keep_nans=%s; keep_infs=%s): FAIL! Expected %s, found %s",
        worker_id,
        q1,
        q2,
        metric,
        keep_nans,
        keep_infs,
        expected,
        found,
    )
    return False


def compare_query_stats(worker_id: int, q1: str, q2: str, expected, sel: hictkpy.PixelSelector) -> int:
    exact_metrics = ["nnz", "sum", "min", "max", "min"]
    approx_metrics = ["variance", "skewness", "kurtosis"]
    num_failures = 0
    for keep_nans in [True, False]:
        for keep_infs in [True, False]:
            stats_expected = compute_stats(expected, keep_nans, keep_infs)
            stats_found = sel.describe(keep_nans=keep_nans, keep_infs=keep_infs)

            for metric in exact_metrics:
                n1 = stats_expected[metric]
                n2 = stats_found[metric]
                if not compare_metric(worker_id, q1, q2, metric, n1, n2, keep_nans, keep_infs):
                    num_failures += 1

            for metric in approx_metrics:
                n1 = stats_expected[metric]
                n2 = stats_found[metric]
                if not compare_metric(worker_id, q1, q2, metric, n1, n2, keep_nans, keep_infs, 1.0e-4, 1.0e-6):
                    num_failures += 1

    return int(num_failures != 0)


def seed_prng(worker_id: int, seed):
    seed = hash(tuple([worker_id, seed]))
    logging.info("[%d] seed: %d", worker_id, seed)
    random.seed(seed)


def worker(
    worker_id: int,
    path_to_file: pathlib.Path,
    path_to_reference_file: pathlib.Path,
    resolution: int,
    chroms_flat,
    chrom_ranks,
    query_type: str,
    query_length_mu: float,
    query_length_std: float,
    _1d_to_2d_query_ratio: float,
    balance: str,
    seed: int | None,
    end_time,
    early_return,
) -> Tuple[int, int]:
    setup_logger(logging.INFO)

    if seed is None:
        seed = random.randint(0, 2**64)

    num_failures = 0
    num_queries = 0

    clr_matrix_args = {
        "balance": balance if balance != "NONE" else False,
    }

    try:
        if query_type == "df":
            clr_matrix_args["as_pixels"] = True
            clr_matrix_args["join"] = True
        elif query_type == "csr":
            clr_matrix_args["sparse"] = True
        elif query_type == "numpy":
            pass
        elif query_type == "describe":
            clr_matrix_args["as_pixels"] = True
        else:
            raise NotImplementedError

        seed_prng(worker_id, seed)

        chrom_sizes = np.array([n for _, n in chroms_flat], dtype=int)
        weights = chrom_sizes / chrom_sizes.sum()

        clr = cooler.Cooler(str(path_to_reference_file))
        sel = clr.matrix(**clr_matrix_args)

        with hictkpy.File(path_to_file, resolution) as f:
            while time.time() < end_time:
                if early_return.value:
                    logging.debug(
                        "[%d] early return signal received. Returning immediately!",
                        worker_id,
                    )
                    break

                if _1d_to_2d_query_ratio <= random.random():
                    q1, q2 = generate_query_2d(
                        chroms_flat,
                        weights,
                        chrom_ranks,
                        mean_length=query_length_mu,
                        stddev_length=query_length_std,
                    )
                else:
                    q1 = generate_query_1d(
                        chroms_flat,
                        weights,
                        mean_length=query_length_mu,
                        stddev_length=query_length_std,
                    )
                    q2 = q1

                num_queries += 1

                expected = cooler_dump(sel, q1, q2)

                if query_type == "describe":
                    num_failures += compare_query_stats(
                        worker_id,
                        q1,
                        q2,
                        expected,
                        f.fetch(q1, q2, normalization=balance),
                    )
                else:
                    found = hictk_dump(f, q1, q2, balance, query_type)
                    num_failures += compare_query_results(worker_id, q1, q2, expected, found)

    except:  # noqa
        logging.debug(
            "[%d] exception raised in worker process. Sending early return signal!",
            worker_id,
        )
        early_return.value = True
        raise

    return num_queries, num_failures


def main() -> int:
    args = vars(make_cli().parse_args())

    if cooler.fileops.is_multires_file(str(args["reference-uri"])) and args["resolution"] is None:
        raise RuntimeError("--resolution is required when test-uri or reference-uri are multi-resolution files.")

    reference_uri = str(args["reference-uri"])
    resolution = args["resolution"]

    if cooler.fileops.is_multires_file(reference_uri):
        if resolution is None:
            raise RuntimeError(
                "--resolution is a mandatory option when one or both of test-uri and reference-uri are multi-resolution files."
            )
        reference_uri = f"{reference_uri}::/resolutions/{resolution}"

    if cooler.fileops.is_cooler(reference_uri) and resolution is not None:
        found_res = cooler.Cooler(reference_uri).binsize
        if found_res is None:
            found_res = 0
        if found_res != resolution:
            raise RuntimeError(
                f'Cooler at "{reference_uri}" has an unexpected resolution: expected {resolution}, found {found_res}.'
            )

    chroms = read_chrom_sizes_cooler(reference_uri)

    masked_chroms = args["excluded_chroms"]
    if masked_chroms is not None:
        for chrom in masked_chroms:
            chroms.pop(chrom, None)

        if len(chroms) == 0:
            raise RuntimeError("No chromosomes left after dropping chromosomes specified through --excluded-chroms")

    logging.info(f"Will generate random queries using the following chromosome(s): {', '.join(chroms.keys())}")

    chrom_ranks = {chrom: i for i, chrom in enumerate(chroms.keys())}
    chroms_flat = list(chroms.items())

    end_time = time.time() + args["duration"]

    with mp.Pool(args["nproc"]) as pool, mp.Manager() as manager:
        early_return = manager.Value(bool, False)
        worker_fx = functools.partial(
            worker,
            path_to_file=args["test-uri"],
            path_to_reference_file=reference_uri,
            resolution=resolution,
            chroms_flat=chroms_flat,
            chrom_ranks=chrom_ranks,
            query_type=args["format"],
            query_length_mu=args["query_length_avg"],
            query_length_std=args["query_length_std"],
            _1d_to_2d_query_ratio=args["1d_to_2d_query_ratio"],
            balance=args["normalization"],
            seed=args["seed"],
            end_time=end_time,
            early_return=early_return,
        )
        results = pool.map(
            worker_fx,
            range(1, args["nproc"] + 1),
            chunksize=1,
        )

    num_queries = sum((n for n, _ in results))
    num_failures = sum((n for _, n in results))
    num_passes = num_queries - num_failures
    if num_failures == 0:
        lvl = logging.INFO
    else:
        lvl = logging.WARN

    logging.log(
        lvl,
        "Score: %.4g/100 (%d success and %d failures).",
        100 * num_passes / num_queries,
        num_passes,
        num_failures,
    )

    return num_failures != 0


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)
    hictkpy.logging.setLevel(logging.WARN)


if __name__ == "__main__":
    mp.set_start_method("spawn")
    setup_logger()
    sys.exit(main())
