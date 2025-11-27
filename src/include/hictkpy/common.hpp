// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <hictk/common.hpp>
#include <hictk/type_traits.hpp>

namespace hictkpy {

#define HICTKPY_LIKELY   HICTK_LIKELY
#define HICTKPY_UNLIKELY HICTK_UNLIKELY

[[noreturn]] inline void unreachable_code() { hictk::unreachable_code(); }

template <typename T, typename U>
[[maybe_unused]] [[nodiscard]] constexpr T conditional_static_cast(U value) {
  return hictk::conditional_static_cast<T>(value);
}

template <typename T>
using remove_cvref = hictk::remove_cvref<T>;

template <typename T>
using remove_cvref_t = typename remove_cvref<T>::type;

}  // namespace hictkpy
