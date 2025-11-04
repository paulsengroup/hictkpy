import logging
import sys
from concurrent.futures import ThreadPoolExecutor

import hictkpy


def fx(i):
    if i % 2 == 0:
        path = "../hictk/test/data/cooler/4DNFIZ1ZVXC8.mcool"
    else:
        path = "../hictk/test/data/hic/4DNFIZ1ZVXC8.hic9"

    print(f"{i} fetching from {path}...", file=sys.stderr)
    with hictkpy.File(path, 50_000) as f:
        res = f.fetch().to_csr().sum()
    print(f"{i} DONE!", file=sys.stderr)
    return res


def test(threads):
    with ThreadPoolExecutor(threads) as pool:
        res = [pool.submit(fx, i + 1) for i in range(threads)]

    return [r.result() for r in res]


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    setup_logger(logging.DEBUG)
    # input("Press enter")
    print(test(16))
