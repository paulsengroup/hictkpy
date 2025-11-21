#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


import argparse
import importlib
import os
import sys

from nanobind.stubgen import StubGen


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()
    cli.add_argument("output-dir", type=str)
    cli.add_argument("--force", default=False, action="store_true")

    return cli


def process_module(mod_name: str, out_file: str, force: bool):
    print(f'Processing module "{mod_name}"...', file=sys.stderr)
    if not force and os.path.exists(out_file):
        raise RuntimeError(f'Refusing to overwrite file "{out_file}". Pass --force to overwrite.')

    mod = importlib.import_module(mod_name)
    sg = StubGen(mod)
    sg.put(mod)

    print(f'Writing stub to file "{out_file}"...', file=sys.stderr)
    with open(out_file, "w") as f:
        print(sg.get(), file=f)


def touch_file(path: str, force: bool):
    if not force and os.path.exists(path):
        raise RuntimeError(f'Refusing to overwrite file "{path}". Pass --force to overwrite.')

    print(f'Touching file "{path}"...', file=sys.stderr)
    with open(path, "w") as f:
        print("", file=f)


def main():
    args = vars(make_cli().parse_args())

    # Ensure that the current directory is on the path
    if "" not in sys.path and "." not in sys.path:
        sys.path.insert(0, "")

    touch_file(os.path.join(args["output-dir"], "py.typed"), args["force"])

    process_module("hictkpy._hictkpy", os.path.join(args["output-dir"], "__init__.pyi"), args["force"])
    for mod in ("cooler", "hic", "logging"):
        process_module(
            mod_name=f"hictkpy._hictkpy.{mod}",
            out_file=os.path.join(args["output-dir"], f"{mod}.pyi"),
            force=args["force"],
        )


if __name__ == "__main__":
    main()
