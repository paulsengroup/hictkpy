// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <patchlevel.h>

#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/variant.hpp"

namespace hictkpy {

template <typename T>
[[nodiscard]] constexpr std::string_view type_to_str() {
  // NOLINTBEGIN(*-avoid-magic-numbers)
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
  if constexpr (std::is_floating_point_v<T>) {
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
  }
  // NOLINTEND(*-avoid-magic-numbers)

  throw std::runtime_error("Unable to infer dtype.");
}

[[nodiscard]] constexpr bool typing_union_required() noexcept {
#ifdef PY_VERSION_HEX
  constexpr auto python_310_hex = 0x030A00F0;
  return PY_VERSION_HEX < python_310_hex;
#else
  return true;
#endif
}

[[nodiscard]] NumericDtype map_py_numeric_to_cpp_type(const nanobind::type_object& dtype);
[[nodiscard]] NumericDtype map_py_numeric_to_cpp_type(std::string_view dtype);
[[nodiscard]] std::string format_py_type(const nanobind::handle& h);

}  // namespace hictkpy
