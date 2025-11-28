# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib


def numpy_avail() -> bool:
    try:
        import numpy
    except ModuleNotFoundError:
        return False

    return True


def pandas_avail() -> bool:
    try:
        import pandas
    except ModuleNotFoundError:
        return False

    return True


def pyarrow_avail() -> bool:
    try:
        import pyarrow
    except ModuleNotFoundError:
        return False

    return True


def scipy_avail() -> bool:
    try:
        import scipy
    except ModuleNotFoundError:
        return False

    return True


def get_test_dir() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parent.parent
