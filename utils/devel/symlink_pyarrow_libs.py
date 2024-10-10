#!/usr/bin/env python3

import glob
import os
import platform
import sys

import pyarrow

"""
This script is largely based on pyarrow.create_library_symlinks()
https://github.com/apache/arrow/blob/c557fe51b3763b9492392f48c5ebcae5a1dd0b42/python/pyarrow/__init__.py#L347-L387

The main difference is that this script detects and fixes broken symlinks.
This script will become redundant once https://github.com/apache/arrow/pull/44228 is merged.
"""


if __name__ == "__main__":
    if platform.system() == "Windows":
        sys.exit(0)

    package_cwd = os.path.dirname(pyarrow.__file__)

    if platform.system() == "Linux":
        bundled_libs = glob.glob(os.path.join(package_cwd, "*.so.*"))

        def get_symlink_path(hard_path):
            return hard_path.rsplit(".", 1)[0]

    else:
        bundled_libs = glob.glob(os.path.join(package_cwd, "*.*.dylib"))

        def get_symlink_path(hard_path):
            return ".".join((hard_path.rsplit(".", 2)[0], "dylib"))

    for lib_hard_path in bundled_libs:
        symlink_path = get_symlink_path(lib_hard_path)
        if os.path.exists(symlink_path):
            continue

        try:
            if os.path.islink(symlink_path):
                # broken symlink
                os.unlink(symlink_path)

            os.symlink(lib_hard_path, symlink_path)
        except PermissionError:
            print(
                "Tried creating symlink {}. If you need to link to "
                "bundled shared libraries, run "
                "pyarrow.create_library_symlinks() as root"
            )
