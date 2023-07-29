#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini (roberros@uio.no)
# SPDX-License-Identifier: MIT

set -e
set -u

# readlink -f is not available on macos...
function readlink_py {
  set -eu
  python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
}

install_prefix="$(readlink_py "${1-/usr/local}")"

conan install conanfile.txt \
  -s build_type=Release \
  -s compiler.cppstd=17 \
  --output-folder "$install_prefix/share/cmake/" \
  -o '*/*:shared=True' \
  --build="missing"
