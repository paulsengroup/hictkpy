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

export MACOSX_DEPLOYMENT_TARGET='10.15'
CMAKE_BUILD_PARALLEL_LEVEL="$(python3 -c 'import multiprocessing; multiprocessing.cpu_count()')"
export CMAKE_BUILD_PARALLEL_LEVEL

install_prefix="$(readlink_py "${1-/usr/local}")"
wd="${TMPDIR-/tmp/}/cibw_setup_deps"
mkdir -p "$wd"

# shellcheck disable=SC2064
trap "cd '$PWD'" EXIT

conan install conanfile.txt \
  -s build_type=Release \
  -s compiler.cppstd=17 \
  --output-folder "$install_prefix/share/cmake/" \
  -o '*/*:shared=True' \
  --build="missing"

data_dir="$(readlink_py ../../../external)"
cd "$wd"


# tar -xf "$data_dir/fast_float-v5.2.0.tar.xz"
# cmake -DCMAKE_BUILD_TYPE=Release \
#       -DCMAKE_INSTALL_PREFIX="$install_prefix" \
#       -S fast_float* \
#       -B fast_float_build
# cmake --build fast_float_build
# cmake --install fast_float_build


# tar -xf "$data_dir/fmt-v10.0.0.tar.xz"
# cmake -DCMAKE_BUILD_TYPE=Release \
#       -DCMAKE_INSTALL_PREFIX="$install_prefix" \
#       -DFMT_TEST=OFF \
#       -DFMT_LIB_DIR="lib" \
#       -DFMT_INSTALL=ON \
#       -S fmt* \
#       -B fmt_build
# cmake --build fmt_build
# cmake --install fmt_build


# tar -xf "$data_dir/spdlog-v1.12.0.tar.xz"
# cmake -DCMAKE_BUILD_TYPE=Release \
#       -DSPDLOG_FMT_EXTERNAL_HO=ON \
#       -DSPDLOG_INSTALL=ON \
#       -DSPDLOG_BUILD_SHARED=OFF \
#       -DCMAKE_INSTALL_PREFIX="$install_prefix" \
#       -S spdlog* \
#       -B spdlog_build
# cmake --build spdlog_build
# cmake --install spdlog_build


tar -xf "$data_dir/zlib-v1.2.13.tar.xz"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$install_prefix" \
      -DSKIP_INSTALL_FILES=ON \
      -S zlib* \
      -B zlib_build
cmake --build zlib_build
cmake --install zlib_build


tar -xf "$data_dir/hdf5-v1.14.1.tar.xz"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$install_prefix" \
      -DBUILD_STATIC_LIBS=OFF \
      -DONLY_SHARED_LIBS=ON \
      -DHDF5_ENABLE_THREADSAFE=OFF \
      -DBUILD_TESTING=OFF \
      -DHDF5_BUILD_TOOLS=OFF \
      -DHDF5_BUILD_EXAMPLES=OFF \
      -DHDF5_BUILD_HL_LIB=OFF \
      -DHDF5_BUILD_FORTRAN=OFF \
      -DHDF5_BUILD_CPP_LIB=OFF \
      -DHDF5_ENABLE_Z_LIB_SUPPORT=ON \
      -DHDF5_ENABLE_SZIP_SUPPORT=OFF \
      -DHDF5_PACKAGE_EXTLIBS=OFF \
      -DZLIB_ROOT=staging \
      -S hdf5* \
      -B hdf5_build
cmake --build hdf5_build
cmake --install hdf5_build

tar -xf "$data_dir/highfive-v2.7.1.tar.xz"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$install_prefix" \
      -DHIGHFIVE_PARALLEL_HDF5=OFF \
      -DHIGHFIVE_USE_BOOST=OFF \
      -DHIGHFIVE_USE_EIGEN=OFF \
      -DHIGHFIVE_USE_XTENSOR=OFF \
      -DHIGHFIVE_USE_OPENCV=OFF \
      -DHIGHFIVE_EXAMPLES=OFF \
      -DHIGHFIVE_UNIT_TESTS=OFF \
      -DHIGHFIVE_BUILD_DOCS=OFF \
      -DHIGHFIVE_USE_INSTALL_DEPS=OFF \
      -DHDF5_C_LIBRARIES="$install_prefix/lib/libhdf5.so" \
      -DHDF5_NO_FIND_PACKAGE_CONFIG_FILE=ON \
      -S HighFive* \
      -B HighFive_build
cmake --build HighFive_build
cmake --install HighFive_build


tar -xf "$data_dir/libdeflate-v1.18.tar.xz"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$install_prefix" \
      -DLIBDEFLATE_BUILD_STATIC_LIB=OFF \
      -DLIBDEFLATE_BUILD_SHARED_LIB=ON \
      -DLIBDEFLATE_BUILD_GZIP=OFF \
      -S libdeflate* \
      -B libdeflate_build
cmake --build libdeflate_build
cmake --install libdeflate_build


# tar -xf "$data_dir/parallel-hashmap-v1.3.11.tar.xz"
# cmake -DCMAKE_BUILD_TYPE=Release \
#       -DCMAKE_INSTALL_PREFIX="$install_prefix" \
#       -DPHMAP_INSTALL=ON \
#       -DPHMAP_BUILD_TESTS=OFF \
#       -DPHMAP_BUILD_EXAMPLES=OFF \
#       -S parallel-hashmap* \
#       -B parallel-hashmap_build
# cmake --build parallel-hashmap_build
# cmake --install parallel-hashmap_build
