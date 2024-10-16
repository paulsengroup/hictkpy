// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/common.hpp"

#include <fmt/format.h>

#include <cstdint>
#include <hictk/numeric_variant.hpp>
#include <stdexcept>
#include <string_view>

namespace hictkpy {

hictk::internal::NumericVariant map_dtype_to_type(std::string_view dtype) {
  static_assert(sizeof(unsigned) == 4);
  static_assert(sizeof(int) == 4);
  static_assert(sizeof(float) == 4);
  static_assert(sizeof(double) == 8);

  if (dtype == "uint8") {
    return std::uint8_t{};
  }
  if (dtype == "uint16") {
    return std::uint16_t{};
  }
  if (dtype == "uint32" || dtype == "uint") {
    return std::uint32_t{};
  }
  if (dtype == "uint64") {
    return std::uint64_t{};
  }

  if (dtype == "int8") {
    return std::int8_t{};
  }
  if (dtype == "int16") {
    return std::int16_t{};
  }
  if (dtype == "int32" || dtype == "int") {
    return std::int32_t{};
  }
  if (dtype == "int64") {
    return std::int64_t{};
  }

  if (dtype == "float16") {
    return float{};
  }
  if (dtype == "float32") {
    return float{};
  }
  if (dtype == "float64" || dtype == "float" || dtype == "double") {
    return double{};
  }

  throw std::runtime_error(fmt::format(
      FMT_STRING(
          "Unable to map dtype \"{}\" to a C++ numeric type. Valid types are: uint, int, float, "
          "double, uint8, uint16, uint32, uint64, int8, int16, int32, int64, float32, and "
          "float64."),
      dtype));
}

}  // namespace hictkpy
