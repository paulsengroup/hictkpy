// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <arrow/type_fwd.h>

#include <algorithm>
#include <array>

namespace hictkpy {

constexpr bool is_dictionary_dtype(arrow::Type::type type) noexcept {
  return type == arrow::Type::DICTIONARY;
}

constexpr bool is_string_dtype(arrow::Type::type type) noexcept {
  using T = decltype(type);
  // clang-format off
  constexpr std::array<arrow::Type::type, 3> string_types{
    T::STRING,
    T::STRING_VIEW,
    T::LARGE_STRING
  };
  // clang-format on
  const auto* match = std::find(string_types.begin(), string_types.end(), type);
  return match != string_types.end();
}

constexpr bool is_integral_dtype(arrow::Type::type type) noexcept {
  using T = decltype(type);
  // clang-format off
  constexpr std::array<arrow::Type::type, 8> integral_types{
    T::UINT8, T::UINT16, T::UINT32, T::UINT64,
    T::INT8,  T::INT16,  T::INT32,  T::INT64
  };
  // clang-format on
  const auto* match = std::find(integral_types.begin(), integral_types.end(), type);
  return match != integral_types.end();
}

constexpr bool is_floating_point_dtype(arrow::Type::type type) noexcept {
  using T = decltype(type);
  return type == T::FLOAT || type == T::DOUBLE;
}

constexpr bool is_numeric_dtype(arrow::Type::type type) noexcept {
  return is_integral_dtype(type) || is_floating_point_dtype(type);
}

}  // namespace hictkpy
