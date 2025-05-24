// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cctype>
#include <hictk/common.hpp>
#include <hictk/type_traits.hpp>
#include <string>
#include <string_view>

#include "hictkpy/nanobind.hpp"

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

[[nodiscard]] inline std::string normalize_log_lvl(std::string_view lvl) {
  std::string normalized_lvl{lvl};
  std::transform(normalized_lvl.begin(), normalized_lvl.end(), normalized_lvl.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return normalized_lvl;
}

}  // namespace hictkpy
