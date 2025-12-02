// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <arrow/compute/cast.h>
#include <arrow/compute/initialize.h>
#include <arrow/type.h>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <memory>
#include <mutex>
#include <stdexcept>
#include <string_view>
#include <tuple>

#include "hictkpy/common.hpp"
#include "hictkpy/type.hpp"

namespace hictkpy {

namespace internal {
inline void init_arrow_compute() {
  static std::once_flag flag;  // NOLINT(*-const-correctness)
  std::call_once(flag, []() {
    const auto status = arrow::compute::Initialize();
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }
  });
}
}  // namespace internal

template <typename... ChunkedArrays>
[[nodiscard]] inline auto normalize_non_uniform_column_types(
    const std::shared_ptr<arrow::DataType> &result_type, const ChunkedArrays &...arrays) {
  const auto first_array = [](auto first, auto &&...) { return first; }(arrays...);

  auto dtype_is_different = [&first_array](const auto &a) {
    return a->type()->id() != first_array->type()->id();
  };

  const bool casting_needed = (dtype_is_different(arrays) || ...);

  if (!casting_needed) {
    return std::make_tuple(arrays...);
  }

  auto cast_array = [&](const auto &array) {
    internal::init_arrow_compute();

    if (array->type()->id() == result_type->id()) {
      return array;
    }
    SPDLOG_DEBUG(FMT_STRING("casting array from {} to {}..."), array->type()->ToString(),
                 result_type->ToString());
    auto res = arrow::compute::Cast(array, result_type);
    if (!res.ok()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("failed to cast array of type {} to type {}: {}"),
                      array->type()->ToString(), result_type->ToString(), res.status().message()));
    }
    return res.ValueUnsafe().chunked_array();
  };

  return std::make_tuple(cast_array(arrays)...);
}

template <typename N>
[[nodiscard]] inline ThinPixelBufferVar allocate_thin_pixel_buffer(std::size_t capacity) {
  ThinPixelBuffer<N> buff;
  buff.reserve(capacity);
  return buff;
}

template <typename N_OUT>
template <typename N_IN,
          typename std::enable_if_t<std::is_integral_v<N_OUT> && std::is_integral_v<N_IN>> *>
constexpr bool SafeNumericConverter<N_OUT>::can_convert(N_IN count) noexcept {
  if constexpr (std::is_same_v<N_IN, N_OUT>) {
    return true;
  }

  constexpr auto same_signedness = std::is_signed_v<N_IN> == std::is_signed_v<N_OUT>;
  // casting to a wider type with the same sign is lossless and safe
  if constexpr (same_signedness && sizeof(N_OUT) >= sizeof(N_IN)) {
    return true;
  }

  // N_IN is wider than N_OUT: need to check at runtime if the conversion is safe
  if constexpr (same_signedness) {
    using NarrowType = N_OUT;
    using WideType = N_IN;
    // NOLINTNEXTLINE(*-signed-char-misuse, *-str34-c)
    constexpr auto lb = static_cast<WideType>(std::numeric_limits<NarrowType>::lowest());
    constexpr auto ub = static_cast<WideType>(std::numeric_limits<NarrowType>::max());

    return count >= lb && count <= ub;  // NOLINT(*-redundant-expression)
  }

  // N_OUT is unsigned while N_IN is not: casting is unsafe is count is negative
  if constexpr (std::is_unsigned_v<N_OUT> && std::is_signed_v<N_IN>) {
    if (count < 0) {
      return false;
    }
  }

  // N_OUT can represent bigger numbers than N_IN: casting is safe
  constexpr auto ub_out = conditional_static_cast<std::uint64_t>(std::numeric_limits<N_OUT>::max());
  constexpr auto ub_in = conditional_static_cast<std::uint64_t>(std::numeric_limits<N_IN>::max());
  if constexpr (ub_out >= ub_in) {
    return true;
  }

  return conditional_static_cast<std::uint64_t>(count) <= ub_out;
}

template <typename N_OUT>
template <typename N_IN, typename std::enable_if_t<std::is_floating_point_v<N_OUT> &&
                                                   std::is_floating_point_v<N_IN>> *>
constexpr bool SafeNumericConverter<N_OUT>::can_convert([[maybe_unused]] N_IN count) noexcept {
  // casting between different fp types is lossy, but safe for our purposes
  return true;
}

template <typename N_OUT>
template <typename N_IN, typename std::enable_if_t<std::is_floating_point_v<N_OUT> !=
                                                   std::is_floating_point_v<N_IN>> *>
constexpr bool SafeNumericConverter<N_OUT>::can_convert([[maybe_unused]] N_IN count) noexcept {
  // casting non-fp to fp is lossy, but safe for our purposes
  return std::is_floating_point_v<N_OUT> && !std::is_floating_point_v<N_IN>;
}

template <typename N_OUT>
template <typename N_IN>
inline N_OUT SafeNumericConverter<N_OUT>::convert(N_IN count) {
  if constexpr (std::is_floating_point_v<N_IN> && !std::is_floating_point_v<N_OUT>) {
    if (!std::isfinite(count)) {
      throw std::invalid_argument("number cannot be converted safely");
    }

    N_IN decimal = std::round(count);
    const auto lb = static_cast<N_IN>(std::numeric_limits<std::int64_t>::lowest());
    const auto ub = static_cast<N_IN>(std::numeric_limits<std::uint64_t>::max());

    if (decimal >= lb && decimal <= ub) {
      if (decimal < 0) {
        return convert(static_cast<std::int64_t>(decimal));
      }
      return convert(static_cast<std::uint64_t>(decimal));
    }
    throw std::invalid_argument("number cannot be converted safely");
  }

  if (can_convert(count)) {
    return conditional_static_cast<N_OUT>(count);
  }
  throw std::invalid_argument("number cannot be converted safely");
}

template <typename N_OUT, typename N_IN>
[[nodiscard]] inline N_OUT safe_numeric_cast(N_IN n) {
  return SafeNumericConverter<N_OUT>::convert(n);
}

template <typename T_OUT, typename T_IN>
[[nodiscard]] inline T_OUT safe_numeric_cast(std::string_view field_name, T_IN n) {
  try {
    return safe_numeric_cast<T_OUT>(n);
  } catch (const std::invalid_argument &) {
    throw std::invalid_argument(fmt::format(FMT_STRING("unable to safely convert {}={} ({}) to {}"),
                                            field_name, n, type_to_str<T_IN>(),
                                            type_to_str<T_OUT>()));
  }
}

}  // namespace hictkpy
