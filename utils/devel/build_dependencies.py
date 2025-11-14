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
import tempfile
from typing import Any, Dict, List


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

    cli.add_argument(
        "--deploy",
        action="store_true",
        default=False,
        help="Deploy dependencies under the output prefix folder.",
    )

    cli.add_argument(
        "--ci",
        action="store_true",
        default=False,
        help="Install dependencies as required by the CI pipeline.",
    )

    cli.add_argument(
        "--no-update",
        action="store_true",
        default=False,
        help="Do not pass --update to conan install/create.",
    )

    return cli


def get_tempdir() -> pathlib.Path:
    return pathlib.Path(tempfile.gettempdir())


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


def add_output_params(args: List, dest: pathlib.Path, deploy: bool):
    if deploy:
        return args + [
            f"--deployer-folder",
            dest,
            "--deployer=full_deploy",
            "--conf=tools.deployer:symlinks=False",
            "--output-folder",
            dest / "cmake",
        ]

    return args + [
        "--output-folder",
        dest,
    ]


def run_conan(
    profile: str,
    build_type: str,
    shared: bool,
    cppstd: str,
    output_prefix: pathlib.Path,
    dry_run: bool,
    no_update: bool,
    deploy: bool,
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
        "--options",
        f"*/*:shared={shared}",
    ]

    if not no_update:
        args.append("--update")

    args = add_output_params(args, output_folder, deploy)

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

    no_update = args["no_update"]
    deploy = args["deploy"]

    dry_run = args["dry_run"]

    output_prefix = infer_root_dir() / "conan-envs"
    if not args["dry_run"]:
        if output_prefix.exists():
            shutil.rmtree(output_prefix)

        output_prefix.mkdir(exist_ok=True)
        (output_prefix / ".gitignore").write_text("*")

    for args in itertools.product(profiles, build_types, shared_build):
        run_conan(
            *args,
            cppstd=cppstd,
            output_prefix=output_prefix,
            deploy=deploy,
            dry_run=dry_run,
            no_update=no_update,
        )


def run_conan_profile_detect_ci(
    name: str = "default",
    cc: str | None = None,
    cxx: str | None = None,
    libcxx: str | None = None,
) -> pathlib.Path:
    assert "CONAN_HOME" in os.environ

    args = ["conan", "profile", "detect", "--name", name]

    env = os.environ.copy()
    if cc is not None:
        env["CC"] = cc
    if cxx is not None:
        env["CXX"] = cxx

    run_or_print(args, env=env, dry_run=False)

    profile = pathlib.Path(os.getenv("CONAN_HOME")) / "profiles" / name
    if any(x is not None for x in (cc, cxx, libcxx)):
        data = []
        for line in profile.open("r"):
            if line.startswith("compiler.libcxx"):
                continue
            data.append(line)

        data.append("[buildenv]")
        if cc is not None:
            data.append(f"CC={cc}")
        if cxx is not None:
            data.append(f"CXX={cxx}")

        profile.write_text("\n".join(data))

    return profile


def run_conan_install_b2_ci(version: str = "5.3.1"):
    args = ["conan", "install", "--requires", f"b2/{version}", "--build=*"]
    if platform.system() == "Linux":
        args.extend(("--profile:b=gcc", "--profile:h=gcc"))

    run_or_print(args, env=os.environ.copy(), dry_run=False)


def run_ci_linux(args: Dict[str, Any]) -> str:
    run_conan_profile_detect_ci(name="gcc", cc="gcc", cxx="g++", libcxx="libstdc++11")
    run_conan_install_b2_ci()

    cppstd = args["cppstd"]
    conan_args = [
        "conan",
        "install",
        "conanfile.py",
        "--build",
        "*",
        "--profile:b=clang",
        "--profile:h=clang",
        "--settings",
        f"compiler.cppstd={cppstd}",
        "--settings",
        f"build_type=Release",
        "--options=*/*:shared=False",
    ]

    run_conan_profile_detect_ci(name="clang", cc="clang", cxx="clang++", libcxx="libstdc++11")

    dest = get_tempdir() / "deps" / "cmake"
    conan_args = add_output_params(conan_args, dest.parent, True)
    run_or_print(conan_args, os.environ.copy(), False)

    return str(dest)


def run_ci_macos(args: Dict[str, Any]) -> str:
    run_conan_profile_detect_ci(name="default")

    cppstd = args["cppstd"]
    conan_args = [
        "conan",
        "install",
        "conanfile.py",
        "--build",
        "missing",
        "--profile:b",
        "default",
        "--profile:h",
        "default",
        "--settings",
        f"compiler.cppstd={cppstd}",
        "--settings",
        "build_type=Release",
        "--settings:h",
        "os.version=14.0",
    ]

    dest = get_tempdir() / "deps" / "cmake"
    conan_args = add_output_params(conan_args, dest.parent, True)
    run_or_print(conan_args, os.environ.copy(), False)

    return str(dest)


def run_ci_windows(args: Dict[str, Any]) -> str:
    winsdk_version = os.getenv("WINSDK_VERSION")
    if winsdk_version is None:
        raise RuntimeError("Unable to read Windows SDK version from environment variable WINSDK_VERSION")

    run_conan_profile_detect_ci(name="default")
    run_conan_install_b2_ci()

    cppstd = args["cppstd"]
    conan_args = [
        "conan",
        "install",
        "conanfile.py",
        "--build",
        "missing",
        "--profile:b",
        "default",
        "--profile:h",
        "default",
        "--settings",
        f"compiler.cppstd={cppstd}",
        "--settings",
        "compiler.runtime_type=Release",
        "--settings",
        "build_type=Release",
        "--conf:a",
        f"tools.microsoft:winsdk_version={winsdk_version}",
    ]

    dest = get_tempdir() / "deps" / "cmake"
    conan_args = add_output_params(conan_args, dest.parent, True)
    run_or_print(conan_args, os.environ.copy(), False)

    return str(dest)


def run_ci(args: Dict[str, Any]):
    if "CONAN_HOME" not in os.environ:
        raise RuntimeError("Please set the CONAN_HOME environment variable")

    if platform.system() == "Linux":
        run_ci_linux(args)
        return

    if platform.system() == "Darwin":
        run_ci_macos(args)
        return

    if platform.system() == "Windows":
        run_ci_windows(args)
        return

    raise RuntimeError(f"Platform {platform.system()} is not supported!")


def main():
    args = vars(make_cli().parse_args())

    if args["ci"]:
        run_ci(args)
    else:
        run_local(args)


if __name__ == "__main__":
    main()
