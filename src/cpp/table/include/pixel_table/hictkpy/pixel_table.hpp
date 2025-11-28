// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <hictk/bin_table.hpp>
#include <memory>

#include "hictkpy/pixel_table_helpers.hpp"
#include "hictkpy/table.hpp"
#include "hictkpy/variant.hpp"

namespace arrow {
class DataType;
}

namespace hictkpy {

// NOLINTBEGIN(*-redundant-declaration)
template <typename... ChunkedArrays>
[[nodiscard]] auto normalize_non_uniform_column_types(
    const std::shared_ptr<arrow::DataType> &result_type, const ChunkedArrays &...arrays);

[[nodiscard]] ThinPixelBufferVar convert_table_to_thin_pixels(const hictk::BinTable &bin_table,
                                                              const PyArrowTable &df, bool sort,
                                                              const NumericDtype &count_type);

template <typename N>
[[nodiscard]] ThinPixelBufferVar allocate_thin_pixel_buffer(std::size_t capacity);
[[nodiscard]] ThinPixelBufferVar allocate_thin_pixel_buffer(std::size_t capacity,
                                                            const NumericDtype &count_type);
// NOLINTEND(*-redundant-declaration)

}  // namespace hictkpy

#include "../../pixel_table_impl.hpp"
