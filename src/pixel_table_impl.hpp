// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <arrow/compute/cast.h>
#include <arrow/type.h>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <memory>
#include <stdexcept>
#include <tuple>

namespace hictkpy {

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
    if (!dtype_is_different(array)) {
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

}  // namespace hictkpy
