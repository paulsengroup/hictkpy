// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <arrow/table.h>
#include <fmt/format.h>

#include <memory>
#include <stdexcept>
#include <type_traits>
#include <variant>

#include "hictkpy/common.hpp"
#include "hictkpy/table.hpp"
#include "hictkpy/variant.hpp"

namespace hictkpy::internal {

[[nodiscard]] inline NumericDtype infer_count_type(const std::shared_ptr<arrow::Table>& df) {
  const auto type = infer_column_dtype(df, "count");
  return std::visit(
      [&df](const auto& x) -> NumericDtype {
        using T = remove_cvref_t<decltype(x)>;
        if constexpr (std::is_arithmetic_v<T>) {
          return x;
        } else {
          const auto col = df->GetColumnByName("count");
          if (!col) {
            throw std::runtime_error(
                "unable to infer dtype for column \"count\": column does not exist!");
          }
          throw std::runtime_error(fmt::format(
              FMT_STRING("unable to infer dtype for column \"count\": unable to map type "
                         "\"{}\" to a known numeric type"),
              col->type()->ToString()));
        }
      },
      type);
}

[[noreturn]] inline void raise_invalid_table_format() {
  // clang-format off
  throw std::invalid_argument(
      "DataFrame is not in COO or BG2 format.\n"
      "Please make sure that the DataFrame contains the following columns:\n"
      "- COO: [bin1_id, bin2_id, count]\n"
      "- BG2: [chrom1, start1, end1, chrom2, start2, end2, count]\n"
      "And that columns have appropriate dtypes:\n"
      "[chrom1, chrom2] -> string/categorical[string]\n"
      "[bin1_id, bin2_id, start1, end1, start2, end2] -> integral\n"
      "[count] -> numeric (excluding complex numbers)");
  // clang-format on
}

}  // namespace hictkpy::internal
