#!/usr/bin/env python3

# Copyright (c) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import importlib
import inspect
import pathlib
import re
import shutil
import subprocess as sp
import sys
from typing import Dict


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def existing_file(arg):
        if (path := pathlib.Path(arg)).exists():
            return path

        raise FileNotFoundError(f'File "{arg}" does not exists')

    cli.add_argument(
        "conanfile",
        type=existing_file,
        help="Path to the conanfile.py to use as input.",
    )

    return cli


def import_conanfile_from_path(file_path: pathlib.Path, module_name="conanfile"):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def extract_requirements(conanfile: pathlib.Path) -> Dict[str, str]:
    conanfile = import_conanfile_from_path(conanfile, conanfile.stem)
    source = inspect.getsource(conanfile.HictkpyConan.requirements)

    pattern = re.compile(r"^\"([\w\-_]+)/([\w.\-_]+)#(\w+)\"(.*)$")

    requirements = {}
    for line in source.split("\n"):
        line = line.strip()
        if not line.startswith("self.requires("):
            continue

        line = line.removeprefix("self.requires(")

        matches = pattern.search(line).groups()
        if len(matches) < 2:
            raise RuntimeError(f'Failed to parse requirements from line "{line}"')

        name = matches[0]
        version = matches[1]
        requirements[name] = version

    return requirements


@functools.cache
def get_conan() -> pathlib.Path:
    conan = shutil.which("conan")
    if not conan:
        raise RuntimeError("Unable to find Conan in your path")

    return pathlib.Path(conan)


def init_profile():
    sp.check_call([get_conan(), "profile", "detect", "--force"], stdout=sp.DEVNULL, stderr=sp.DEVNULL)


def get_conan_home() -> pathlib.Path:
    return pathlib.Path(sp.check_output([get_conan(), "config", "home"]).decode("utf-8"))


def main():
    args = vars(make_cli().parse_args())
    init_profile()

    profile = get_conan_home() / "profiles" / "default"
    assert profile.is_file()

    requires = extract_requirements(args["conanfile"])
    replace_requires = ["arrow", "boost", "libdeflate", "hdf5", "zstd"]

    with profile.open("a") as f:
        f.write("\n[platform_requires]\n")

        for pkg in replace_requires:
            f.write(f"{pkg}/{requires[pkg]}\n")

        f.write("\n[replace_requires]\n")

        for pkg in replace_requires:
            f.write(f"{pkg}/{requires[pkg]}: {pkg}/[*]\n")

        f.write("\n")


if __name__ == "__main__":
    main()
