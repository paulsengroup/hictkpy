#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import itertools
import os
import pathlib
import platform
import shlex
import shutil
import subprocess as sp
import tarfile
from typing import Any, Dict


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "-o",
        "--output-prefix",
        type=pathlib.Path,
        help="Path where to store the resulting *.cmake files",
    )
    cli.add_argument(
        "--profile",
        nargs="+",
        type=str,
        default=("gcc", "clang"),
        choices={"gcc", "clang", "default"},
        help="Names of the conan profiles to be used.",
    )
    cli.add_argument(
        "--build-type",
        nargs="+",
        type=str,
        default=("Debug", "RelWithDebInfo", "Release"),
        choices={"Debug", "Release", "RelWithDebInfo", "MinSizeRel"},
        help="Conan build types.",
    )
    cli.add_argument(
        "--cppstd",
        type=str,
        default="17",
        choices={"17", "20", "23"},
        help="C++ standard used to compile the dependencies.",
    )
    cli.add_argument(
        "--build-shared-only",
        action="store_true",
        default=False,
        help="Build dependencies as shared libraries only.",
    )
    cli.add_argument(
        "--build-static-only",
        action="store_true",
        default=False,
        help="Build dependencies as static libraries only.",
    )
    cli.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="Print the commands that would be executed if --dry-run was not specified, then exit.",
    )

    return cli


def infer_root_dir() -> pathlib.Path:
    path = pathlib.Path(sp.check_output(["git", "rev-parse", "--show-toplevel"], encoding="utf-8").strip())

    if not path.is_dir():
        raise RuntimeError("Unable to infer repository root!")

    return path


def run_or_print(args, env: Dict[str, str], dry_run: bool):
    if dry_run:
        print(shlex.join(str(x) for x in args))
    else:
        sp.check_call(args, env=env)


def run_conan(
    profile: str,
    build_type: str,
    shared: bool,
    cppstd: str,
    output_prefix: pathlib.Path,
    dry_run: bool,
):
    output_folder = output_prefix / profile / build_type / ("shared" if shared else "static")

    args = [
        "conan",
        "install",
        infer_root_dir() / "conanfile.py",
        "--build=missing",
        "--update",
        "--profile",
        profile,
        "--settings",
        f"build_type={build_type}",
        "--settings",
        f"compiler.cppstd={cppstd}",
        "--output-folder",
        output_folder,
        "--options",
        f"*/*:shared={shared}",
    ]

    env = os.environ.copy()
    env["CC"] = profile
    if profile == "gcc":
        env["CXX"] = "g++"
    elif profile == "clang":
        env["CXX"] = "clang++"
    elif profile != "default":
        raise RuntimeError(
            f'Unrecognized compiler "{profile}". Profiles should be either named "gcc", "clang", or "default".'
        )

    run_or_print(args, env=env, dry_run=dry_run)


def run_local(args: Dict[str, Any]):

    profiles = args["profile"]
    build_types = args["build_type"]
    cppstd = args["cppstd"]
    if args["build_shared_only"]:
        shared_build = [True]
    elif args["build_static_only"]:
        shared_build = [False]
    else:
        shared_build = [True, False]

    dry_run = args["dry_run"]

    output_prefix = infer_root_dir() / "conan-envs"
    if not args["dry_run"]:
        if output_prefix.exists():
            shutil.rmtree(output_prefix)

        output_prefix.mkdir(exist_ok=True)

    for args in itertools.product(profiles, build_types, shared_build):
        run_conan(
            *args,
            cppstd=cppstd,
            output_prefix=output_prefix,
            dry_run=dry_run,
        )


def main():
    args = vars(make_cli().parse_args())

    run_local(args)


if __name__ == "__main__":
    main()
