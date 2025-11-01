import json
import logging
import pathlib
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import hictkpy


def runme(dest):
    t0 = time.time()
    with hictkpy.cooler.FileWriter(dest, chromosomes=chroms, resolution=resolution) as w:
        w.add_pixels(pixels)

    return {
        "tid": threading.current_thread().name,
        "timestamp": time.ctime(),
        "duration": time.time() - t0,
    }


def printstuff(results):
    t0 = time.time()
    while any(r.running() for r in results):
        print(f"{threading.current_thread().name}: {time.time() - t0:.1f}s", file=sys.stderr)
        time.sleep(0.5)

    return [r.result() for r in results]


def bench_mt(nthreads=8):
    assert nthreads > 0
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmpdir:
        with ThreadPoolExecutor(nthreads) as pool:
            tasks = []
            for i in range(nthreads):
                path = pathlib.Path(tmpdir) / f"{i}.cool"
                tasks.append(pool.submit(runme, path))
            results = printstuff(tasks)

    return time.time() - t0, results


def bench_st():
    t0 = time.time()
    t = threading.Thread(target=runme)
    t.start()
    t.join()
    return time.time() - t0


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    setup_logger(0)

    with hictkpy.File("../hictk/test/data/cooler/4DNFIZ1ZVXC8.mcool", 100_000) as f:
        pixels = f.fetch().to_df()
        resolution = f.resolution()
        chroms = f.chromosomes()

    input("Press enter")

    results = []
    for i in range(5):
        duration, res = bench_mt(int(sys.argv[1]))
        print(f"{i}\t{duration}", file=sys.stderr)
        results.append({"iteration": i, "results": res})

    print("DONE!", file=sys.stderr)

    json.dump(
        results,
        sys.stdout,
        indent=2,
    )
