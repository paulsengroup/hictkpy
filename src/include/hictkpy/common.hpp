// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <hictk/common.hpp>
#include <hictk/numeric_variant.hpp>
#include <hictk/type_traits.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

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

template <typename T>
[[nodiscard]] constexpr std::string_view map_type_to_dtype() {
  if constexpr (std::is_unsigned_v<T>) {
    switch (sizeof(T)) {
      case 1:
        return "uint8";
      case 2:
        return "uint16";
      case 4:
        return "uint32";
      case 8:
        return "uint64";
      default:
        unreachable_code();
    }
  }
  if constexpr (std::is_integral_v<T>) {
    switch (sizeof(T)) {
      case 1:
        return "int8";
      case 2:
        return "int16";
      case 4:
        return "int32";
      case 8:
        return "int64";
      default:
        unreachable_code();
    }
  }
  switch (sizeof(T)) {
    case 2:
      return "float16";
    case 4:
      return "float32";
    case 8:
      return "float64";
    default:
      unreachable_code();
  }

  throw std::runtime_error("Unable to infer numpy dtype.");
}

[[nodiscard]] hictk::internal::NumericVariant map_dtype_to_type(std::string_view dtype);

[[nodiscard]] std::string normalize_log_lvl(std::string_view lvl);

}  // namespace hictkpy
