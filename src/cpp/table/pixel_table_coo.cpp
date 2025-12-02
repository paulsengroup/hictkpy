// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <arrow/array/array_base.h>
#include <arrow/array/array_primitive.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <hictk/pixel.hpp>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

#include "hictkpy/common.hpp"
#include "hictkpy/pixel_table.hpp"
#include "hictkpy/pixel_table_helpers.hpp"
#include "hictkpy/variant.hpp"

namespace hictkpy::coo {

template <typename ArrowBin, typename ArrowCount, typename Count>
static void process_chunk(const std::shared_ptr<arrow::NumericArray<ArrowBin>> &bin1_ids_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &bin2_ids_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          ThinPixelBuffer<Count> &buff) {
  const auto size = bin1_ids_chunk->data()->length;
  assert(bin2_ids_chunk->data()->length == size);
  assert(counts_chunk->data()->length == size);

  if (size == 0) {
    return;
  }

  for (std::int64_t i = 0; i < size; ++i) {
    const auto b1 = bin1_ids_chunk->Value(i);
    const auto b2 = bin2_ids_chunk->Value(i);
    const auto count = counts_chunk->Value(i);

    if constexpr (std::is_signed_v<decltype(b1)>) {
      if (b1 < 0) {
        throw std::runtime_error("found negative value in bin1_id column");
      }
      if (b2 < 0) {
        throw std::runtime_error("found negative value in bin2_id column");
      }
    }
    buff.emplace_back(hictk::ThinPixel<Count>{safe_numeric_cast<std::uint64_t>("bin1_id", b1),
                                              safe_numeric_cast<std::uint64_t>("bin2_id", b2),
                                              safe_numeric_cast<Count>("count", count)});
  }
}

template <typename ArrowBin, typename ArrowCount>
static void process_chunk(const std::shared_ptr<arrow::NumericArray<ArrowBin>> &bin1_ids_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &bin2_ids_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          ThinPixelBufferVar &buff) {
  std::visit(
      [&](auto &buff_) { process_chunk(bin1_ids_chunk, bin2_ids_chunk, counts_chunk, buff_); },
      buff);
}

template <typename ArrowCount>
static void process_chunk(const std::shared_ptr<arrow::Array> &bin1_ids_chunk,
                          const std::shared_ptr<arrow::Array> &bin2_ids_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          ThinPixelBufferVar &buff) {
  assert(bin1_ids_chunk->type()->id() == bin2_ids_chunk->type()->id());

  const auto array_type_var = [&]() {
    try {
      return integer_array_static_pointer_cast_helper(bin1_ids_chunk->type());
    } catch (const std::exception &e) {
      throw std::invalid_argument(
          fmt::format(FMT_STRING("failed to infer dtype for bin{{1,2}}_id columns: {}"), e.what()));
    } catch (...) {
      throw std::invalid_argument("failed to infer dtype for bin{{1,2}}_id columns: unknown error");
    }
  }();

  return std::visit(
      [&]([[maybe_unused]] const auto &array) {
        using T = remove_cvref_t<decltype(*array)>;
        process_chunk(std::static_pointer_cast<T>(bin1_ids_chunk),
                      std::static_pointer_cast<T>(bin2_ids_chunk), counts_chunk, buff);
      },
      array_type_var);
}

static void process_chunk(const std::shared_ptr<arrow::Array> &bin1_ids_chunk,
                          const std::shared_ptr<arrow::Array> &bin2_ids_chunk,
                          const std::shared_ptr<arrow::Array> &counts_chunk,
                          ThinPixelBufferVar &buff) {
  const auto array_type_var = [&]() {
    try {
      return numeric_array_static_pointer_cast_helper(counts_chunk->type());
    } catch (const std::exception &e) {
      throw std::invalid_argument(
          fmt::format(FMT_STRING("failed to infer dtype for count column: {}"), e.what()));
    } catch (...) {
      throw std::invalid_argument("failed to infer dtype for count column: unknown error");
    }
  }();

  std::visit(
      [&]([[maybe_unused]] const auto &array) {
        using T = remove_cvref_t<decltype(*array)>;
        process_chunk(bin1_ids_chunk, bin2_ids_chunk, std::static_pointer_cast<T>(counts_chunk),
                      buff);
      },
      array_type_var);
}

ThinPixelBufferVar convert_table_thin_pixels(std::shared_ptr<arrow::Table> df, bool sort,
                                             const NumericDtype &count_type) {
  try {
    if (!df) {
      return allocate_thin_pixel_buffer(0, count_type);
    }

    df = ensure_table_has_uniform_chunks(df);

    // we assume that the array types have already been validated
    auto [bin1_ids, bin2_ids] = normalize_non_uniform_column_types(
        arrow::int64(), df->GetColumnByName("bin1_id"), df->GetColumnByName("bin2_id"));
    auto counts = df->GetColumnByName("count");

    auto buffer = allocate_thin_pixel_buffer(static_cast<std::size_t>(df->num_rows()), count_type);

    const auto num_chunks = counts->num_chunks();
    for (std::int32_t i = 0; i < num_chunks; ++i) {
      const auto chunk1 = bin1_ids->chunk(i);
      const auto chunk2 = bin2_ids->chunk(i);
      const auto chunk3 = counts->chunk(i);
      process_chunk(chunk1, chunk2, chunk3, buffer);
    }

    std::visit(
        [&](auto &buffer_) {
          assert(buffer_.size() == static_cast<std::size_t>(df->num_rows()));
          df.reset();

          if (sort) {
            std::sort(buffer_.begin(), buffer_.end());
          }
        },
        buffer);

    return buffer;
  } catch (const std::invalid_argument &e) {
    throw std::invalid_argument(
        fmt::format(FMT_STRING("failed to convert DataFrame to a COO pixels: {}"), e.what()));

  } catch (const std::runtime_error &e) {
    throw std::invalid_argument(
        fmt::format(FMT_STRING("failed to convert DataFrame to a COO pixels: {}"), e.what()));
  } catch (...) {
    throw std::invalid_argument("failed to convert DataFrame to a COO pixels: unknown error");
  }
}

}  // namespace hictkpy::coo
