// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <hictk/numeric_variant.hpp>
#include <stdexcept>
#include <string_view>
#include <type_traits>

#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"

namespace hictkpy {

template <typename T>
[[nodiscard]] constexpr std::string_view type_to_str() {
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

  throw std::runtime_error("Unable to infer dtype.");
}

[[nodiscard]] hictk::internal::NumericVariant map_py_type_to_cpp_type(
    const nanobind::type_object& dtype);
[[nodiscard]] hictk::internal::NumericVariant map_py_type_to_cpp_type(std::string_view dtype);

}  // namespace hictkpy
