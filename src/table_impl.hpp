// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <arrow/type_fwd.h>

namespace hictkpy {

constexpr bool is_dictionary_dtype(arrow::Type::type type) noexcept {
  return type == arrow::Type::DICTIONARY;
}

constexpr bool is_string_dtype(arrow::Type::type type) noexcept {
  using T = decltype(type);
  return type == T::STRING || type == T::STRING_VIEW || type == T::LARGE_STRING;
}

constexpr bool is_integral_dtype(arrow::Type::type type) noexcept {
  switch (type) {
    using T = decltype(type);
    case T::UINT8:
      return true;
    case T::UINT16:
      return true;
    case T::UINT32:
      return true;
    case T::UINT64:
      return true;
    case T::INT8:
      return true;
    case T::INT16:
      return true;
    case T::INT32:
      return true;
    case T::INT64:
      return true;
    default:
      return false;
  }
}

constexpr bool is_floating_point_dtype(arrow::Type::type type) noexcept {
  using T = decltype(type);
  return type == T::FLOAT || type == T::DOUBLE;
}

constexpr bool is_numeric_dtype(arrow::Type::type type) noexcept {
  return is_integral_dtype(type) || is_floating_point_dtype(type);
}

}  // namespace hictkpy
