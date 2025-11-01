from concurrent.futures import ThreadPoolExecutor

import hictkpy


def fx():
    f = hictkpy.File("../hictk/test/data/cooler/4DNFIZ1ZVXC8.mcool", 10_000)
    return f.fetch().sum()


def test(threads):
    with ThreadPoolExecutor(threads) as pool:
        res = [pool.submit(fx), pool.submit(fx)]

    return [r.result() for r in res]


if __name__ == "__main__":
    test(2)
