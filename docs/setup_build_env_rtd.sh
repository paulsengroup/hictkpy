#!/usr/bin/env bash
# shellcheck disable=SC2016

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


set -e
set -o pipefail
set -u
set -x


CONDA_PREFIX="$(
  mamba info --envs |
    grep '^base' |
    sed 's/^.*[[:space:]]\+//' |
    tr -d '\n'
)"

mamba install -c conda-forge --file docs/requirements.conda.txt --file docs/requirements.txt

docs/patch_project_for_rtd.py conanfile.py --root-dir "$PWD" --inplace
docs/patch_project_for_rtd.py pyproject.toml --root-dir "$PWD" --inplace
docs/patch_project_for_rtd.py docs/index.rst --root-dir "$PWD" --inplace

# Help CMake find dependencies installed with conda
PREFIX_PATH="$(printf 'list(APPEND CMAKE_PREFIX_PATH "%s")' "$CONDA_PREFIX")"
sed -i "/APPEND CMAKE_MODULE_PATH/a $PREFIX_PATH" CMakeLists.txt

# Override default linker to avoid problems due to missing LLVMgold.so
sed -i '/APPEND CMAKE_PREFIX_PATH/a set(CMAKE_LINKER_TYPE LLD)' CMakeLists.txt

# Don't bother with optimizing the build
sed -i '/APPEND CMAKE_MODULE_PATH/i set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O0")' CMakeLists.txt
sed -i '/APPEND CMAKE_MODULE_PATH/i set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O0")' CMakeLists.txt
