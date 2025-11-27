// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <arrow/type_fwd.h>

#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/pixel.hpp>
#include <memory>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictkpy/variant.hpp"

namespace hictkpy {

class PyArrowTable;

template <typename N>
using ThinPixelBuffer = std::vector<hictk::ThinPixel<N>>;

// clang-format off
using ThinPixelBufferVar =
  std::variant<
    ThinPixelBuffer<std::uint8_t>,
    ThinPixelBuffer<std::uint16_t>,
    ThinPixelBuffer<std::uint32_t>,
    ThinPixelBuffer<std::uint64_t>,
    ThinPixelBuffer<std::int8_t>,
    ThinPixelBuffer<std::int16_t>,
    ThinPixelBuffer<std::int32_t>,
    ThinPixelBuffer<std::int64_t>,
    ThinPixelBuffer<float>,
    ThinPixelBuffer<double>
  >;

using ArrowNumericArray =
    std::variant<arrow::UInt8Array*, arrow::UInt16Array*, arrow::UInt32Array*, arrow::UInt64Array*,
                 arrow::Int8Array*,  arrow::Int16Array*,  arrow::Int32Array*,  arrow::Int64Array*,
                 arrow::FloatArray*, arrow::DoubleArray*>;

using ArrowIntegerArray =
    std::variant<arrow::UInt8Array*, arrow::UInt16Array*, arrow::UInt32Array*, arrow::UInt64Array*,
                 arrow::Int8Array*,  arrow::Int16Array*,  arrow::Int32Array*,  arrow::Int64Array*>;

using ArrowStringArray =
    std::variant<arrow::StringArray*, arrow::StringViewArray*, arrow::LargeStringArray*>;
// clang-format on

template <typename... ChunkedArrays>
[[nodiscard]] auto normalize_non_uniform_column_types(
    const std::shared_ptr<arrow::DataType> &result_type, const ChunkedArrays &...arrays);

[[nodiscard]] std::shared_ptr<arrow::Table> ensure_table_has_uniform_chunks(
    std::shared_ptr<arrow::Table> table);

[[nodiscard]] ArrowNumericArray numeric_array_static_pointer_cast_helper(
    const std::shared_ptr<arrow::DataType> &dtype);

[[nodiscard]] ArrowIntegerArray integer_array_static_pointer_cast_helper(
    const std::shared_ptr<arrow::DataType> &dtype);

[[nodiscard]] ArrowStringArray string_array_static_pointer_cast_helper(
    const std::shared_ptr<arrow::DataType> &dtype);

namespace bg2 {

[[nodiscard]] ThinPixelBufferVar convert_table_thin_pixels(const hictk::BinTable &bins,
                                                           std::shared_ptr<arrow::Table> df,
                                                           bool sort,
                                                           const NumericDtype &count_type);
}

namespace coo {

[[nodiscard]] ThinPixelBufferVar convert_table_thin_pixels(std::shared_ptr<arrow::Table> df,
                                                           bool sort,
                                                           const NumericDtype &count_type);
}

[[nodiscard]] ThinPixelBufferVar convert_table_thin_pixels(const hictk::BinTable &bins,
                                                           const PyArrowTable &df, bool sort,
                                                           const NumericDtype &count_type);

template <typename N_OUT>
struct SafeNumericConverter {
  static_assert(std::is_arithmetic_v<N_OUT>);

  SafeNumericConverter() = delete;

  template <typename N_IN, typename std::enable_if_t<std::is_integral_v<N_OUT> &&
                                                     std::is_integral_v<N_IN>> * = nullptr>
  static constexpr bool can_convert(N_IN count) noexcept;
  template <typename N_IN, typename std::enable_if_t<std::is_floating_point_v<N_OUT> &&
                                                     std::is_floating_point_v<N_IN>> * = nullptr>
  static constexpr bool can_convert(N_IN count) noexcept;
  template <typename N_IN, typename std::enable_if_t<std::is_floating_point_v<N_OUT> !=
                                                     std::is_floating_point_v<N_IN>> * = nullptr>
  static constexpr bool can_convert(N_IN count) noexcept;
  template <typename N_IN>
  [[nodiscard]] static N_OUT convert(N_IN count);
};

template <typename N_OUT, typename N_IN>
[[nodiscard]] N_OUT safe_numeric_cast(N_IN n);

template <typename T_OUT, typename T_IN>
[[nodiscard]] T_OUT safe_numeric_cast(std::string_view field_name, T_IN n);

}  // namespace hictkpy

#include "../../pixel_table_impl.hpp"
