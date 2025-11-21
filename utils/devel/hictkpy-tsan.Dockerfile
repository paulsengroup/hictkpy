# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# docker build -f utils/devel/hictkpy-tsan.Dockerfile . -t hictkpy-tsan

FROM ghcr.io/paulsengroup/ci-docker-images/ubuntu-python-314 AS base

ARG CC=clang
ARG CXX=clang++

RUN apt-get update \
&&  apt-get install -y -q build-essential clang git \
&&  rm -rf /var/lib/apt/lists/*

COPY . /tmp/hictkpy/

RUN sed -i 's/^cmake\.build-type = .*$/cmake.build-type = "RelWithDebInfo"/' /tmp/hictkpy/pyproject.toml \
&&  sed -i '/tool\.scikit-build\.cmake\.define/a CMAKE_C_FLAGS = "-fsanitize=thread"' /tmp/hictkpy/pyproject.toml \
&&  sed -i '/tool\.scikit-build\.cmake\.define/a CMAKE_CXX_FLAGS = "-fsanitize=thread"' /tmp/hictkpy/pyproject.toml \
&&  sed -i '/tool\.scikit-build\.cmake\.define/a CMAKE_SHARED_LINKER_FLAGS = "-fsanitize=thread"' /tmp/hictkpy/pyproject.toml

RUN /opt/python/tsan/bin/python3 -m venv /opt/python/venv \
&&  /opt/python/venv/bin/pip install '/tmp/hictkpy/[all]' -v \
&&  /opt/python/venv/bin/python -c 'import hictkpy; print(hictkpy.__version__)'

FROM ghcr.io/paulsengroup/ci-docker-images/ubuntu-python-314 AS final

RUN apt-get update \
&&  apt-get install -y -q llvm \
&&  rm -rf /var/lib/apt/lists/*

RUN which llvm-symbolizer

ENV PATH="/opt/python/venv/bin:$PATH"

COPY --from=base /opt/python/venv/ /opt/python/venv

RUN python -c 'import hictkpy; print(hictkpy.__version__)'
