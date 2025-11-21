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
import tempfile
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Dict, Tuple

import pandas as pd
import pyarrow as pa
import scipy.sparse as ss
from numpy.typing import NDArray

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
        if (n := int(arg)) > 1:
            return n

        raise ValueError("Not a number greater than 1")

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
    cli.add_argument(
        "--range1",
        type=str,
        help="Genomic range in UCSC format used to perform queries.",
    )
    cli.add_argument(
        "--range2",
        type=str,
        help="Genomic range in UCSC format used to perform queries.",
    )
    return cli


def hash_pandas_df(df: pd.DataFrame) -> int:
    data = []
    for name, col in df.items():
        data.extend(
            (
                name,
                col.dtype,
                col.to_numpy().tobytes(),
            )
        )

    # logging.debug("%s", data)
    return hash(tuple(data))


def hash_arrow_table(df: pa.Table) -> int:
    df = df.combine_chunks()

    def collect_column(i):
        field = df.field(i)

        data = [field.name, field.type]
        for chunk in df.column(i).chunks:
            if chunk is not None:
                data.append(chunk.to_numpy().tobytes())

        return tuple(data)

    data = tuple(collect_column(i) for i in range(df.num_columns))
    # logging.debug("%s", data)
    return hash(data)


def hash_numpy_array(v: NDArray) -> int:
    logging.debug("%s %s", v, v.shape)
    return hash(
        (
            v.dtype,
            v.shape,
            v.tobytes(),
        )
    )


def hash_csr_matrix(m: ss.csr_matrix) -> int:
    logging.debug("shape=%s nnz=%s sum=%s dtype=%s", m.shape, m.nnz, m.sum(), m.dtype)
    return hash(
        (
            m.dtype,
            m.shape,
            m.data.tobytes(),
            m.indices.tobytes(),
            m.indptr.tobytes(),
        )
    )


def hash_coo_matrix(m: ss.coo_matrix) -> int:
    logging.debug("shape=%s nnz=%s sum=%s dtype=%s", m.shape, m.nnz, m.sum(), m.dtype)
    return hash(
        (
            m.dtype,
            m.shape,
            m.data.tobytes(),
            m.row.tobytes(),
            m.col.tobytes(),
        )
    )


def hash_hictkpy_file(f: hictkpy.File, range1: str | None = None, range2: str | None = None) -> int:
    attrs = f.attributes()
    attrs.pop("creation-date", None)
    logging.debug(
        "chroms=%s resolution=%d is_hic=%d is_cooler=%s attributes=%s sum=%d",
        json.dumps(f.chromosomes(), sort_keys=True),
        f.resolution(),
        f.is_hic(),
        f.is_cooler(),
        json.dumps(attrs, sort_keys=True),
        f.fetch(range1, range2).sum(),
    )
    return hash(
        (
            json.dumps(f.chromosomes(), sort_keys=True),
            f.resolution(),
            f.is_hic(),
            f.is_cooler(),
            json.dumps(attrs, sort_keys=True),
            hash_csr_matrix(f.fetch(range1, range2).to_csr()),
        )
    )


def call_to_arrow(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    df = f.fetch(range1, range2).to_arrow()
    return hash_arrow_table(df)


def call_to_pandas(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    df = f.fetch(range1, range2).to_pandas()
    return hash_pandas_df(df)


def call_to_df(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    df = f.fetch(range1, range2).to_df()
    return hash_pandas_df(df)


def call_to_numpy(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    m = f.fetch(range1, range2).to_numpy()
    return hash_numpy_array(m)


def call_to_csr(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    m = f.fetch(range1, range2).to_csr()
    return hash_csr_matrix(m)


def call_to_coo(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    m = f.fetch(range1, range2).to_coo()
    return hash_coo_matrix(m)


def call_describe(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    return hash(json.dumps(f.fetch(range1, range2).describe(), sort_keys=True))


def call_nnz(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    return f.fetch(range1, range2).nnz()


def call_sum(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    return hash(f.fetch(range1, range2).sum())


def call_min(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    return hash(f.fetch(range1, range2).min())


def call_max(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    return hash(f.fetch(range1, range2).max())


def call_variance(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    return hash(f.fetch(range1, range2).variance())


def call_skewness(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    return hash(f.fetch(range1, range2).skewness())


def call_kurtosis(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    return hash(f.fetch(range1, range2).kurtosis())


def call_iterator(f: hictkpy.File, range1: str | None, range2: str | None) -> int:
    return hash(tuple(p.count for p in f.fetch(range1, range2)))


def create_file_helper(
    chroms: Dict[str, int],
    resolution: int,
    pixel_chunks: Tuple[pd.DataFrame, ...],
    tmpdir: pathlib.Path,
    format: str,
) -> int:

    if format == "hic":
        Writer = hictkpy.hic.FileWriter
    elif format == "cool":
        Writer = hictkpy.cooler.FileWriter
    else:
        raise NotImplementedError()

    with tempfile.TemporaryDirectory(dir=tmpdir) as tmpdir:
        tmpdir = pathlib.Path(tmpdir)

        dest = tmpdir / f"test.{format}"
        with Writer(
            path=dest,
            chromosomes=chroms,
            resolution=resolution,
            tmpdir=tmpdir,
            compression_lvl=1,
        ) as w:
            for chunk in pixel_chunks:
                w.add_pixels(chunk)

        with hictkpy.File(dest) as f:
            return hash_hictkpy_file(f)


def create_cooler(
    chroms: Dict[str, int],
    resolution: int,
    pixel_chunks: Tuple[pd.DataFrame, ...],
    tmpdir: pathlib.Path,
    format: str = "cool",
) -> int:
    assert format == "cool"
    return create_file_helper(
        chroms,
        resolution,
        pixel_chunks,
        tmpdir,
        format="cool",
    )


def create_hic(
    chroms: Dict[str, int],
    resolution: int,
    pixel_chunks: Tuple[pd.DataFrame, ...],
    tmpdir: pathlib.Path,
    format: str = "hic",
) -> int:
    assert format == "hic"
    return create_file_helper(
        chroms,
        resolution,
        pixel_chunks,
        tmpdir,
        format="hic",
    )


@functools.cache
def vtable() -> Dict[str, Callable]:
    names = [
        "call_to_arrow",
        "call_to_pandas",
        "call_to_df",
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
        "call_iterator",
        "create_cooler",
        "create_hic",
    ]

    return {name: globals()[name] for name in names}


def runner(
    path: pathlib.Path,
    range1: str | None,
    range2: str | None,
    resolution: int,
    end_time: float,
    seed: int | None,
) -> Dict[str, int]:
    if seed is None:
        seed = threading.get_ident()
    prng = random.Random(seed)
    counters = {name: 0 for name in vtable()}
    counters["errors"] = 0

    if range2 is None:
        range2 = range1

    with (
        hictkpy.File(path, resolution) as f,
        tempfile.TemporaryDirectory() as tmpdir,
    ):
        tmpdir = pathlib.Path(tmpdir)

        file_type = "cool" if f.is_cooler() else "hic"
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

            if name.startswith("create_"):
                if file_type == "cool":
                    name = "create_cooler"
                else:
                    name = f"create_{file_type}"
                fx = vtable().get(name)
                chroms, _, chunks = get_create_file_inputs(path, range1, range2, resolution)
                result = fx(chroms, resolution, chunks, tmpdir, file_type)
                expected_result = reference_results[name]
            else:
                result = fx(f, range1, range2)
                expected_result = reference_results[f"{name}_{file_type}"]

            if result != expected_result:
                logging.error(
                    "[tid=%d] [iter=%d]: %s[%s] incorrect result, expected hash=%d, found hash=%d",
                    threading.get_ident(),
                    i,
                    name,
                    file_type,
                    expected_result,
                    result,
                )
                counters["errors"] += 1

            counters[name] += 1
            i += 1

    return counters


def run_fuzzer(
    hic_file: pathlib.Path,
    cooler_file: pathlib.Path,
    range1: str | None,
    range2: str | None,
    resolution: int,
    nthreads: int,
    duration: float,
    seed: int | None,
) -> Dict[str, int]:
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
                    runner,
                    path=path,
                    range1=range1,
                    range2=range2,
                    resolution=resolution,
                    end_time=end_time,
                    seed=seed,
                )
            )

        for t in tasks:
            results += Counter(t.result())

        results = dict(results)
        if "errors" not in results:
            results["errors"] = 0

        results["total"] = sum(v for k, v in results.items() if k != "errors")

        return results


@functools.cache
def get_create_file_inputs(
    path: pathlib.Path,
    range1: str | None,
    range2: str | None,
    resolution: int,
) -> Tuple[Dict[str, int], int, Tuple[pd.DataFrame, ...]]:
    with hictkpy.File(path, resolution) as f:
        df = f.fetch(range1, range2).to_arrow().to_pandas()

        chunk_size = len(df) // 5

        chunks = []
        for start in range(0, len(df), chunk_size):
            end = start + chunk_size
            chunks.append(df[start:end])

        return f.chromosomes(), resolution, tuple(chunks)


def compute_reference_results(
    hic_file: pathlib.Path,
    cooler_file: pathlib.Path,
    range1: str | None,
    range2: str | None,
    resolution: int,
) -> Dict[Tuple[str, str], int]:

    logging.info("generating reference results...")

    results = {}
    with (
        hictkpy.File(hic_file, resolution) as hf,
        hictkpy.File(cooler_file, resolution) as cf,
    ):
        for name, fx in vtable().items():
            if name.startswith("create_"):
                continue

            logging.info("generating reference results for %s()...", name)
            results[f"{name}_hic"] = fx(hf, range1, range2)
            results[f"{name}_cool"] = fx(cf, range1, range2)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)

        chroms, _, pixel_chunks = get_create_file_inputs(hic_file, range1, range2, resolution)

        logging.info("generating reference results for create_hic()...")
        results["create_hic"] = create_hic(chroms, resolution, pixel_chunks, tmpdir)

        chroms, _, pixel_chunks = get_create_file_inputs(cooler_file, range1, range2, resolution)
        logging.info("generating reference results for create_cooler()...")
        results["create_cooler"] = create_cooler(chroms, resolution, pixel_chunks, tmpdir)

    logging.info("### Reference results:\n%s", json.dumps(results, indent=2, sort_keys=True))
    return results


def main() -> int:
    global reference_results

    args = vars(make_cli().parse_args())

    hic_file = args["hic-file"]
    cooler_file = args["cooler-file"]
    duration = args["duration"]

    assert hictkpy.is_hic(hic_file)
    assert hictkpy.is_cooler(cooler_file) or hictkpy.is_mcool_file(cooler_file)

    reference_results = compute_reference_results(
        hic_file,
        cooler_file,
        range1=args["range1"],
        range2=args["range2"],
        resolution=args["resolution"],
    )

    logging.info("fuzzing hictkpy for ~%d seconds", duration)
    results = run_fuzzer(
        hic_file,
        cooler_file,
        range1=args["range1"],
        range2=args["range2"],
        resolution=args["resolution"],
        nthreads=args["nthreads"],
        duration=duration,
        seed=args["seed"],
    )

    json.dump(
        results,
        sys.stdout,
        indent=2,
        sort_keys=True,
    )
    sys.stdout.write("\n")

    return int(results["errors"] != 0)


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)
    hictkpy.logging.setLevel(logging.WARN)


if __name__ == "__main__":
    mp.set_start_method("spawn")
    setup_logger()

    reference_results = {}
    sys.exit(main())
