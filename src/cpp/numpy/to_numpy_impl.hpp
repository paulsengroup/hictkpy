// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"

namespace hictkpy {

template <typename N, typename OwningNumPyVector>
[[nodiscard]] inline OwningNumPyVector make_owning_numpy(std::vector<N> data) {
  auto data_ = std::make_unique<std::vector<N>>(std::move(data));
  auto *data_ptr = data_->data();
  auto size = data_->size();

  HICTKPY_GIL_SCOPED_ACQUIRE
  return OwningNumPyVector{data_ptr, {size}, make_capsule(std::move(data_))};
}

template <typename N_OUT, typename N_IN, typename OwningNumPyVector>
[[nodiscard]] inline OwningNumPyVector make_owning_numpy(const std::vector<N_IN> &data) {
  static_assert(!std::is_same_v<N_IN, N_OUT>);

  auto data_ = std::make_unique<std::vector<N_OUT>>(data.size());

  std::transform(data.begin(), data.end(), data_->begin(),
                 [](const auto res) { return static_cast<N_OUT>(res); });

  auto *data_ptr = data_->data();
  auto size = data_->size();

  HICTKPY_GIL_SCOPED_ACQUIRE
  return OwningNumPyVector{data_ptr, {size}, make_capsule(std::move(data_))};
}

}  // namespace hictkpy
