// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

template <typename N, typename OwningNumPyVector = nanobind::ndarray<
                          nanobind::numpy, nanobind::ndim<1>, nanobind::c_contig, N>>
[[nodiscard]] OwningNumPyVector make_owning_numpy(std::vector<N> data);

template <typename N_OUT, typename N_IN,
          typename OwningNumPyVector =
              nanobind::ndarray<nanobind::numpy, nanobind::ndim<1>, nanobind::c_contig, N_OUT>>
[[nodiscard]] OwningNumPyVector make_owning_numpy(const std::vector<N_IN> &data);

}  // namespace hictkpy

#include "../../to_numpy_impl.hpp"
