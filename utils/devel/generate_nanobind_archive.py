#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


import argparse
import hashlib
import logging
import pathlib
import shlex
import shutil
import subprocess as sp
import tempfile
import textwrap
from tabnanny import check
from typing import Any, List


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def existing_folder(s: str) -> pathlib.Path:
        p = pathlib.Path(s)
        if p.is_dir():
            return p

        raise argparse.ArgumentError(f'"{s}" does not point to an existing folder')

    cli.add_argument(
        "tag",
        type=str,
        help="Git tag used to generate the archive.",
    )
    cli.add_argument(
        "--output-dir",
        type=existing_folder,
        help="Path to an existing folder where to store the archive of nanobind source code.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing file(s).",
    )

    return cli


def get_exec(name: str | pathlib.Path) -> pathlib.Path:
    exec = shutil.which(name)
    if exec:
        return pathlib.Path(exec)

    raise RuntimeError(f"Unable to find {name} in your PATH")


def generate_printable_cmd(cmd: List[Any]) -> str:
    return shlex.join(str(x) for x in cmd)


def guess_repo_root() -> pathlib.Path:
    git = get_exec("git")
    cmd = [git, "rev-parse", "--show-toplevel"]

    logging.info(f"launching {generate_printable_cmd(cmd)}...")
    try:
        return pathlib.Path(sp.check_output(cmd).decode("utf-8").strip())

    except:  # noqa
        pass

    try:
        cwd = pathlib.Path(__file__).parent.resolve()
        return pathlib.Path(sp.check_output(cmd).decode("utf-8").strip(), cwd=cwd)
    except:  # noqa
        raise RuntimeError(
            "Unable to guess the desired output path. Please use --output-dir to specify the desired output folder."
        )


def checkout_nanobind(tag: str, tmpdir: pathlib.Path) -> pathlib.Path:
    git = get_exec("git")
    cloned_repo = tmpdir / f"nanobind-{tag}"

    cmd = [git, "clone", "https://github.com/wjakob/nanobind.git", cloned_repo]
    logging.info(f"launching {generate_printable_cmd(cmd)}...")
    sp.check_call(cmd)

    cmd = [git, "checkout", tag]
    logging.info(f"launching {generate_printable_cmd(cmd)}...")
    sp.check_call(cmd, cwd=cloned_repo)

    cmd = [git, "submodule", "update", "--init", "--recursive"]
    logging.info(f"launching {generate_printable_cmd(cmd)}...")
    sp.check_call(cmd, cwd=cloned_repo)

    shutil.rmtree(cloned_repo / ".git")

    return cloned_repo


def generate_archive(repo: pathlib.Path, tag: str, out_dir: pathlib.Path, force: bool) -> pathlib.Path:
    tar = get_exec("tar")
    xz = get_exec("xz")

    tmpdir = repo.parent

    plain_archive = tmpdir / f"nanobind-{tag}.tar"
    cmd = [tar, "-cf", plain_archive.name, repo.name]
    logging.info(f"launching {generate_printable_cmd(cmd)}...")
    sp.check_call(cmd, cwd=repo.parent)

    dest = out_dir / f"{plain_archive.name}.xz"
    if dest.exists() and not force:
        raise RuntimeError(f"Refusing to overwrite existing file {dest}. Pass --force to overwrite.")

    with dest.open("wb") as f:
        cmd = [xz, "-9", "--extreme", "--stdout", plain_archive]
        logging.info(f"launching {generate_printable_cmd(cmd)}...")
        sp.check_call(cmd, stdout=f)

    return dest


def hash_file(path: pathlib.Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            data = f.read(64 << 10)
            if not data:
                break
            hasher.update(data)

    return hasher.hexdigest()


def main():
    args = vars(make_cli().parse_args())
    git_tag = args["tag"]
    out_dir = args["output_dir"]

    if out_dir is None:
        out_dir = guess_repo_root() / "external"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)

        cloned_repo = checkout_nanobind(git_tag, tmpdir)
        dest = generate_archive(cloned_repo, git_tag, out_dir, args["force"])

    digest = hash_file(dest)

    logging.info(f'Nanobind archive successfully created at "{dest}"!')
    logging.info(f"Archive size: {dest.stat().st_size / 1.0e3}KB")
    logging.info(f"Archive SHA256: {digest}")

    msg = f"""
    FetchContent_Declare(
      nanobind
      URL       "${{CMAKE_CURRENT_SOURCE_DIR}}/external/{dest.name}"
      URL_HASH  "SHA256={digest}"
      EXCLUDE_FROM_ALL
      SYSTEM)
     """

    print(textwrap.dedent(msg))


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    setup_logger()
    main()
