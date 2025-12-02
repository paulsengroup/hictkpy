// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <arrow/array/array_binary.h>
#include <arrow/array/array_dict.h>
#include <arrow/array/array_primitive.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <arrow/type_fwd.h>
#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <hictk/bin_table.hpp>
#include <hictk/pixel.hpp>
#include <hictk/reference.hpp>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "hictkpy/common.hpp"
#include "hictkpy/pixel_table.hpp"
#include "hictkpy/pixel_table_helpers.hpp"
#include "hictkpy/variant.hpp"

namespace hictkpy::bg2 {

class ChromosomeIDArray {
  // clang-format off
  using ChromIDs = std::variant<
    std::vector<std::uint8_t>,
    std::vector<std::uint16_t>,
    std::vector<std::uint32_t>
  >;
  // clang-format on

  ChromIDs _chrom_ids{};
  std::shared_ptr<arrow::Array> _chroms_arrow{};
  const hictk::Reference *_chroms{};

 public:
  ChromosomeIDArray(const hictk::Reference &chroms, std::shared_ptr<arrow::Array> chroms_arrow)
      : _chrom_ids(allocate_chrom_id_buffer(chroms)),
        _chroms_arrow(std::move(chroms_arrow)),
        _chroms(&chroms) {
    encode();
  }

  [[nodiscard]] const ChromIDs &get() const noexcept { return _chrom_ids; }

  template <typename I>
  [[nodiscard]] const std::vector<I> &get() const {
    return std::get<std::vector<I>>(_chrom_ids);
  }

  // NOLINTNEXTLINE(*-exception-escape)
  [[nodiscard]] std::size_t size() const noexcept {
    return std::visit([](const auto &v) { return v.size(); }, _chrom_ids);
  }

  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

 private:
  [[nodiscard]] static ChromIDs allocate_chrom_id_buffer(const hictk::Reference &chroms) {
    if (chroms.size() <= std::numeric_limits<std::uint8_t>::max()) {
      return {std::vector<std::uint8_t>{}};
    }
    if (chroms.size() <= std::numeric_limits<std::uint16_t>::max()) {
      return {std::vector<std::uint16_t>{}};
    }
    assert(chroms.size() <= std::numeric_limits<std::uint32_t>::max());
    return {std::vector<std::uint32_t>{}};
  }

  [[nodiscard]] std::uint32_t encode(std::string_view chrom) { return _chroms->at(chrom).id(); }

  template <typename StringLikeArray>
  void encode(const std::shared_ptr<StringLikeArray> &chroms) {
    using T = remove_cvref_t<StringLikeArray>;
    static_assert(std::is_same_v<T, arrow::StringArray> ||
                  std::is_same_v<T, arrow::LargeStringArray> ||
                  std::is_same_v<T, arrow::StringViewArray>);
    const auto num_records = chroms->length();
    try {
      std::visit(
          [&](auto &chrom_ids) {
            using N = remove_cvref_t<decltype(chrom_ids.front())>;
            chrom_ids.clear();
            chrom_ids.reserve(static_cast<std::size_t>(num_records));
            for (std::int64_t i = 0; i < num_records; ++i) {
              const auto chrom_id = conditional_static_cast<N>(encode(chroms->GetView(i)));
              chrom_ids.push_back(chrom_id);
            }
          },
          _chrom_ids);
    } catch (const std::exception &e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("failed to encode chromosomes: {}"), e.what()));
    } catch (...) {
      throw std::runtime_error("failed to encode chromosomes: unknown error");
    }
  }

  template <typename IndicesArray, typename StringArray>
  void encode(const std::shared_ptr<IndicesArray> &indices,
              const std::shared_ptr<StringArray> &dict) {
    const auto num_records = indices->length();
    try {
      std::visit(
          [&](auto &chrom_ids) {
            using N = remove_cvref_t<decltype(chrom_ids.front())>;
            chrom_ids.clear();
            chrom_ids.reserve(static_cast<std::size_t>(num_records));
            for (std::int64_t i = 0; i < num_records; ++i) {
              const auto dict_idx = conditional_static_cast<std::int64_t>(indices->Value(i));
              const auto chrom_name = dict->GetView(dict_idx);
              const auto chrom_id = conditional_static_cast<N>(encode(chrom_name));
              chrom_ids.push_back(chrom_id);
            }
          },
          _chrom_ids);
    } catch (const std::exception &e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("failed to encode chromosomes: {}"), e.what()));
    } catch (...) {
      throw std::runtime_error("failed to encode chromosomes: unknown error");
    }
  }

  // NOLINTNEXTLINE(*-convert-member-functions-to-static)
  void encode(const std::shared_ptr<arrow::DictionaryArray> &chroms) {
    const auto idx_dtype = chroms->indices()->type();
    const auto dict_dtype = chroms->dictionary()->type();

    if (!is_integral_dtype(idx_dtype->id())) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("unable to decode chromosomes: expected dictionary with index of "
                                 "type integer, found {}"),
                      idx_dtype->ToString()));
    }
    if (!is_string_dtype(dict_dtype->id())) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("unable to decode chromosomes: expected dictionary with value of "
                                 "type string, found {}"),
                      dict_dtype->ToString()));
    }

    const auto idx_dtype_var = integer_array_static_pointer_cast_helper(idx_dtype);
    const auto dict_dtype_var = string_array_static_pointer_cast_helper(dict_dtype);

    std::visit(
        [&]([[maybe_unused]] const auto *dict_dtype_) {
          using DictArray = remove_cvref_t<decltype(*dict_dtype_)>;
          std::visit(
              [&]([[maybe_unused]] const auto *idx_dtype_) {
                using IndicesArray = remove_cvref_t<decltype(*idx_dtype_)>;
                encode(std::static_pointer_cast<IndicesArray>(chroms->indices()),
                       std::make_shared<DictArray>(chroms->dictionary()->data()));
              },
              idx_dtype_var);
        },
        dict_dtype_var);
  }

  void encode() {
    // NOLINTBEGIN(*-avoid-return-with-void-value)
    switch (_chroms_arrow->type()->id()) {
      using T = arrow::Type::type;
      case T::STRING:
        return encode(std::make_shared<arrow::StringArray>(_chroms_arrow->data()));
      case T::STRING_VIEW:
        return encode(std::make_shared<arrow::StringViewArray>(_chroms_arrow->data()));
      case T::LARGE_STRING:
        return encode(std::make_shared<arrow::LargeStringArray>(_chroms_arrow->data()));
      case T::DICTIONARY:
        return encode(std::make_shared<arrow::DictionaryArray>(_chroms_arrow->data()));
      default:
        unreachable_code();
    }
    // NOLINTEND(*-avoid-return-with-void-value)
  }
};

template <typename ArrowBin, typename ArrowCount, typename ChromID, typename Count>
static void process_chunk(const hictk::BinTable &bins, const std::vector<ChromID> &chrom1_id_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end1_chunk,
                          const std::vector<ChromID> &chrom2_id_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          std::vector<hictk::ThinPixel<Count>> &buff) {
  const auto size = static_cast<std::int64_t>(chrom1_id_chunk.size());
  assert(start1_chunk->length() == size);
  assert(end1_chunk->length() == size);
  assert(chrom2_id_chunk.size() == static_cast<std::size_t>(size));
  assert(start2_chunk->length() == size);
  assert(end2_chunk->length() == size);
  assert(counts_chunk->length() == size);

  auto get_bin_checked = [&](std::uint_fast8_t idx, auto chrom_id, auto start, auto end) {
    auto bin = bins.at(chrom_id, start);
    if (bin.end() == end) {
      return bin;
    }
    const auto res = bins.resolution();
    throw std::runtime_error(fmt::format(
        FMT_STRING("invalid end{}: expected {}, found {}, (start{}={}; bin_size={})"), idx,
        bin.end(), end, idx, bin.start(), res == 0 ? "variable" : fmt::to_string(res)));
  };

  try {
    for (std::int64_t i = 0; i < size; ++i) {
      const auto ii = static_cast<std::size_t>(i);
      const auto chrom1_id = conditional_static_cast<std::uint32_t>(chrom1_id_chunk[ii]);
      const auto start1 = start1_chunk->Value(i);
      const auto end1 = end1_chunk->Value(i);
      const auto chrom2_id = conditional_static_cast<std::uint32_t>(chrom2_id_chunk[ii]);
      const auto start2 = start2_chunk->Value(i);
      const auto end2 = end2_chunk->Value(i);
      const auto count = counts_chunk->Value(i);
      try {
        if (start1 < 0 || end1 < 0 || start2 < 0 || end2 < 0) {
          throw std::runtime_error("genomic coordinates cannot be negative");
        }

        if (end1 < start1 || end2 < start2) {
          throw std::runtime_error(
              "end position of a bin cannot be smaller than its start position");
        }

        // clang-format off
        const hictk::Pixel p{
            get_bin_checked(
              1,
              chrom1_id,
              safe_numeric_cast<std::uint32_t>("start1", start1),
              safe_numeric_cast<std::uint32_t>("end1", end1)
            ),
            get_bin_checked(
              2,
              chrom2_id,
              safe_numeric_cast<std::uint32_t>("start2", start2),
              safe_numeric_cast<std::uint32_t>("end2", end2)
            ),
            safe_numeric_cast<Count>("count", count)
        };
        // clang-format on
        buff.emplace_back(p.to_thin());
      } catch (const std::exception &e) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("failed to map {}:{}-{}; {}:{}-{} to a valid pixel: {}"),
                        bins.chromosomes().at(chrom1_id).name(), start1, end1,
                        bins.chromosomes().at(chrom2_id).name(), start1, end1, e.what()));
      } catch (...) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("failed to map {}:{}-{}; {}:{}-{} to a valid pixel: unknown error"),
            bins.chromosomes().at(chrom1_id).name(), start1, end1,
            bins.chromosomes().at(chrom2_id).name(), start1, end1));
      }
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to process BG2 pixel: {}"), e.what()));
  } catch (...) {
    throw std::runtime_error("failed to process BG2 pixel: unknown error");
  }
}

template <typename ArrowBin, typename ArrowCount, typename ChromID>
static void process_chunk(const hictk::BinTable &bins, const std::vector<ChromID> &chrom1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end1_chunk,
                          const std::vector<ChromID> &chrom2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          ThinPixelBufferVar &buff) {
  std::visit(
      [&](auto &buff_) {
        process_chunk(bins, chrom1_chunk, start1_chunk, end1_chunk, chrom2_chunk, start2_chunk,
                      end2_chunk, counts_chunk, buff_);
      },
      buff);
}

template <typename ArrowBin, typename ArrowCount>
static void process_chunk(const hictk::BinTable &bins,
                          const std::shared_ptr<arrow::Array> &chrom1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end1_chunk,
                          const std::shared_ptr<arrow::Array> &chrom2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          ThinPixelBufferVar &buff) {
  const ChromosomeIDArray chrom1_ids(bins.chromosomes(), chrom1_chunk);
  const ChromosomeIDArray chrom2_ids(bins.chromosomes(), chrom2_chunk);

  std::visit(
      [&](const auto &chrom1_ids_buff) {
        using T = remove_cvref_t<decltype(chrom1_ids_buff)>;
        assert(std::holds_alternative<T>(chrom2_ids.get()));

        const auto &chrom2_ids_buff = std::get<T>(chrom2_ids.get());
        process_chunk(bins, chrom1_ids_buff, start1_chunk, end1_chunk, chrom2_ids_buff,
                      start2_chunk, end2_chunk, counts_chunk, buff);
      },
      chrom1_ids.get());
}

template <typename ArrowCount>
static void process_chunk(const hictk::BinTable &bins,
                          const std::shared_ptr<arrow::Array> &chrom1_chunk,
                          const std::shared_ptr<arrow::Array> &start1_chunk,
                          const std::shared_ptr<arrow::Array> &end1_chunk,
                          const std::shared_ptr<arrow::Array> &chrom2_chunk,
                          const std::shared_ptr<arrow::Array> &start2_chunk,
                          const std::shared_ptr<arrow::Array> &end2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          ThinPixelBufferVar &buff) {
  assert(start1_chunk->type()->id() == end1_chunk->type()->id());
  assert(start1_chunk->type()->id() == start2_chunk->type()->id());
  assert(start1_chunk->type()->id() == end2_chunk->type()->id());

  const auto array_type_var = [&]() {
    try {
      return numeric_array_static_pointer_cast_helper(start1_chunk->type());
    } catch (const std::exception &e) {
      throw std::invalid_argument(fmt::format(
          FMT_STRING("failed to infer dtype for start{{1,2}} and end{{1,2}} columns: {}"),
          e.what()));
    } catch (...) {
      throw std::invalid_argument(
          "failed to infer dtype for start{{1,2}} and end{{1,2}} columns: unknown error");
    }
  }();

  return std::visit(
      [&]([[maybe_unused]] const auto &array) {
        using T = remove_cvref_t<decltype(*array)>;
        // clang-format off
        process_chunk(
            bins,
            chrom1_chunk,
            std::static_pointer_cast<T>(start1_chunk),
            std::static_pointer_cast<T>(end1_chunk),
            chrom2_chunk,
            std::static_pointer_cast<T>(start2_chunk),
            std::static_pointer_cast<T>(end2_chunk),
            counts_chunk,
            buff
        );
        // clang-format on
      },
      array_type_var);
}

static void process_chunk(const hictk::BinTable &bins,
                          const std::shared_ptr<arrow::Array> &chrom1_chunk,
                          const std::shared_ptr<arrow::Array> &start1_chunk,
                          const std::shared_ptr<arrow::Array> &end1_chunk,
                          const std::shared_ptr<arrow::Array> &chrom2_chunk,
                          const std::shared_ptr<arrow::Array> &start2_chunk,
                          const std::shared_ptr<arrow::Array> &end2_chunk,
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
        process_chunk(bins, chrom1_chunk, start1_chunk, end1_chunk, chrom2_chunk, start2_chunk,
                      end2_chunk, std::static_pointer_cast<T>(counts_chunk), buff);
      },
      array_type_var);
}

ThinPixelBufferVar convert_table_thin_pixels(const hictk::BinTable &bins,
                                             std::shared_ptr<arrow::Table> df, bool sort,
                                             const NumericDtype &count_type) {
  try {
    if (!df) {
      return allocate_thin_pixel_buffer(0, count_type);
    }
    df = ensure_table_has_uniform_chunks(df);
    auto chrom1 = df->GetColumnByName("chrom1");
    auto chrom2 = df->GetColumnByName("chrom2");
    auto counts = df->GetColumnByName("count");
    auto [start1, end1, start2, end2] = normalize_non_uniform_column_types(
        arrow::int64(), df->GetColumnByName("start1"), df->GetColumnByName("end1"),
        df->GetColumnByName("start2"), df->GetColumnByName("end2"));

    auto buffer = allocate_thin_pixel_buffer(static_cast<std::size_t>(df->num_rows()), count_type);

    const auto num_chunks = counts->num_chunks();
    for (std::int32_t i = 0; i < num_chunks; ++i) {
      auto chunk1 = chrom1->chunk(i);
      auto chunk2 = start1->chunk(i);
      auto chunk3 = end1->chunk(i);
      auto chunk4 = chrom2->chunk(i);
      auto chunk5 = start2->chunk(i);
      auto chunk6 = end2->chunk(i);
      auto chunk7 = counts->chunk(i);
      process_chunk(bins, chunk1, chunk2, chunk3, chunk4, chunk5, chunk6, chunk7, buffer);
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

}  // namespace hictkpy::bg2
