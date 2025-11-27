// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/pixel_table.hpp"

#include <arrow/array/array_binary.h>
#include <arrow/array/array_primitive.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <arrow/type_fwd.h>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

#include "hictkpy/common.hpp"
#include "hictkpy/pixel_table_helpers.hpp"
#include "hictkpy/table.hpp"
#include "hictkpy/variant.hpp"

namespace hictkpy {

std::shared_ptr<arrow::Table> ensure_table_has_uniform_chunks(std::shared_ptr<arrow::Table> table) {
  const auto &columns = table->columns();
  if (columns.size() < 2) {
    return table;
  }

  auto combine_table_chunks = [&]() {
    SPDLOG_DEBUG("found uneven chunks while converting arrow::Table to hictk::ThinPixels");
    auto res = table->CombineChunks();
    if (!res.ok()) {
      throw std::runtime_error(fmt::format(FMT_STRING("failed to combine arrow::Table chunks: {}"),
                                           res.status().message()));
    }
    table.reset();
    return res.MoveValueUnsafe();
  };

  auto check_uneven_chunks = [length = columns.front()->length()](const auto &chunk) {
    return chunk->length() != length;
  };

  auto num_chunks = columns.front()->num_chunks();
  auto table_has_uneven_chunks =
      std::any_of(columns.begin(), columns.end(),
                  [num_chunks](const auto &col) { return col->num_chunks() != num_chunks; });

  if (table_has_uneven_chunks) {
    return combine_table_chunks();
  }

  arrow::ArrayVector chunks(columns.size());
  for (std::int32_t i = 0; i < num_chunks; ++i) {
    chunks.clear();
    for (const auto &col : columns) {
      chunks.emplace_back(col->chunk(i));
    }

    if (std::any_of(chunks.begin(), chunks.end(), check_uneven_chunks)) {
      return combine_table_chunks();
    }
  }

  return table;
}

ArrowNumericArray numeric_array_static_pointer_cast_helper(
    const std::shared_ptr<arrow::DataType> &dtype) {
  using ArrowArray = ArrowNumericArray;

  switch (dtype->id()) {
    using T = arrow::Type::type;
    case T::UINT8:
      return ArrowArray{static_cast<arrow::UInt8Array *>(nullptr)};
    case T::UINT16:
      return ArrowArray{static_cast<arrow::UInt16Array *>(nullptr)};
    case T::UINT32:
      return ArrowArray{static_cast<arrow::UInt32Array *>(nullptr)};
    case T::UINT64:
      return ArrowArray{static_cast<arrow::UInt64Array *>(nullptr)};
    case T::INT8:
      return ArrowArray{static_cast<arrow::Int8Array *>(nullptr)};
    case T::INT16:
      return ArrowArray{static_cast<arrow::Int16Array *>(nullptr)};
    case T::INT32:
      return ArrowArray{static_cast<arrow::Int32Array *>(nullptr)};
    case T::INT64:
      return ArrowArray{static_cast<arrow::Int64Array *>(nullptr)};
    case T::FLOAT:
      return ArrowArray{static_cast<arrow::FloatArray *>(nullptr)};
    case T::DOUBLE:
      return ArrowArray{static_cast<arrow::DoubleArray *>(nullptr)};
    default:
      throw std::invalid_argument(
          fmt::format(FMT_STRING("{} is not a valid numeric dtype"), dtype->ToString()));
  }
}

ArrowIntegerArray integer_array_static_pointer_cast_helper(
    const std::shared_ptr<arrow::DataType> &dtype) {
  using ArrowArray = ArrowIntegerArray;

  return std::visit(
      [&]([[maybe_unused]] const auto *a) -> ArrowArray {
        using Array = remove_cvref_t<decltype(*a)>;
        if constexpr (std::is_same_v<arrow::FloatArray, Array> ||
                      std::is_same_v<arrow::DoubleArray, Array>) {
          throw std::invalid_argument(
              fmt::format(FMT_STRING("{} is not a valid integral dtype"), dtype->ToString()));
        } else {
          return {static_cast<Array *>(nullptr)};
        }
      },
      numeric_array_static_pointer_cast_helper(dtype));
}

ArrowStringArray string_array_static_pointer_cast_helper(
    const std::shared_ptr<arrow::DataType> &dtype) {
  using ArrowArray = ArrowStringArray;
  switch (dtype->id()) {
    using T = arrow::Type::type;
    case T::STRING:
      return ArrowArray{static_cast<arrow::StringArray *>(nullptr)};
    case T::STRING_VIEW:
      return ArrowArray{static_cast<arrow::StringViewArray *>(nullptr)};
    case T::LARGE_STRING:
      return ArrowArray{static_cast<arrow::LargeStringArray *>(nullptr)};
    default:
      throw std::invalid_argument(
          fmt::format(FMT_STRING("{} is not a valid string dtype"), dtype->ToString()));
  }
}

ThinPixelBufferVar convert_table_to_thin_pixels(const hictk::BinTable &bin_table,
                                                const PyArrowTable &df, bool sort,
                                                const NumericDtype &count_type) {
  switch (df.type()) {
    case PyArrowTable::Type::COO:
      return coo::convert_table_thin_pixels(df.get(), sort, count_type);
    case PyArrowTable::Type::BG2:
      return bg2::convert_table_thin_pixels(bin_table, df.get(), sort, count_type);
    default:
      unreachable_code();
  }
}

ThinPixelBufferVar allocate_thin_pixel_buffer(std::size_t capacity,
                                              const NumericDtype &count_type) {
  return std::visit(
      [&]([[maybe_unused]] const auto &type) -> ThinPixelBufferVar {
        using N = remove_cvref_t<decltype(type)>;
        return allocate_thin_pixel_buffer<N>(capacity);
      },
      count_type);
}

}  // namespace hictkpy
