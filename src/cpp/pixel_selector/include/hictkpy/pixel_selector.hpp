// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/cooler/pixel_selector.hpp>
#include <hictk/hic/pixel_selector.hpp>
#include <hictk/pixel.hpp>
#include <hictk/transformers/common.hpp>
#include <hictk/transformers/to_dataframe.hpp>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <variant>
#include <vector>

#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/variant.hpp"

namespace hictkpy {

struct PixelSelector {
  // clang-format off
  using SelectorVar =
    std::variant<std::shared_ptr<const hictk::cooler::PixelSelector>,
                 std::shared_ptr<const hictk::hic::PixelSelector>,
                 std::shared_ptr<const hictk::hic::PixelSelectorAll>>;
  // clang-format on
  using PixelVar = NumericDtype;
  using QuerySpan = hictk::transformers::QuerySpan;
  using PixelFormat = hictk::transformers::DataFrameFormat;
  using DenseMatrix = nanobind::ndarray<nanobind::numpy, nanobind::ndim<2>, nanobind::c_contig>;

  SelectorVar selector{};
  PixelVar pixel_count{std::int32_t{0}};
  hictk::transformers::DataFrameFormat pixel_format{PixelFormat::COO};
  std::optional<std::uint64_t> _diagonal_band_width{};

  PixelSelector() = default;

  PixelSelector(std::shared_ptr<const hictk::cooler::PixelSelector> sel_, PixelVar count_type,
                bool join, std::optional<std::int64_t> diagonal_band_width);
  PixelSelector(std::shared_ptr<const hictk::hic::PixelSelector> sel_, PixelVar count_type,
                bool join, std::optional<std::int64_t> diagonal_band_width);
  PixelSelector(std::shared_ptr<const hictk::hic::PixelSelectorAll> sel_, PixelVar count_type,
                bool join, std::optional<std::int64_t> diagonal_band_width);

  PixelSelector(std::shared_ptr<const hictk::cooler::PixelSelector> sel_,
                const nanobind::type_object& type, bool join,
                std::optional<std::int64_t> diagonal_band_width);
  PixelSelector(std::shared_ptr<const hictk::hic::PixelSelector> sel_,
                const nanobind::type_object& type, bool join,
                std::optional<std::int64_t> diagonal_band_width);
  PixelSelector(std::shared_ptr<const hictk::hic::PixelSelectorAll> sel_,
                const nanobind::type_object& type, bool join,
                std::optional<std::int64_t> diagonal_band_width);

  [[nodiscard]] std::string repr() const;

  using GenomicCoordTuple = std::tuple<std::string, std::int64_t, std::int64_t>;

  [[nodiscard]] auto get_coord1() const -> std::optional<GenomicCoordTuple>;
  [[nodiscard]] auto get_coord2() const -> std::optional<GenomicCoordTuple>;
  [[nodiscard]] std::int64_t size(bool upper_triangular) const;
  [[nodiscard]] nanobind::type_object dtype() const;

  [[nodiscard]] nanobind::iterator make_iterable() const;
  [[nodiscard]] std::optional<nanobind::object> to_arrow(std::string_view span) const;
  [[nodiscard]] nanobind::object to_pandas(std::string_view span) const;
  [[nodiscard]] nanobind::object to_coo(std::string_view span, bool low_memory) const;
  [[nodiscard]] std::optional<nanobind::object> to_csr(std::string_view span,
                                                       bool low_memory) const;
  [[nodiscard]] auto to_numpy(std::string_view span) const -> DenseMatrix;

  [[nodiscard]] nanobind::dict describe(const std::vector<std::string>& metrics, bool keep_nans,
                                        bool keep_infs, bool keep_zeros, bool exact) const;
  [[nodiscard]] std::int64_t nnz(bool keep_nans, bool keep_infs) const;
  [[nodiscard]] std::optional<std::variant<std::int64_t, double>> sum(bool keep_nans,
                                                                      bool keep_infs) const;
  [[nodiscard]] std::optional<std::variant<std::int64_t, double>> min(bool keep_nans,
                                                                      bool keep_infs,
                                                                      bool keep_zeros) const;
  [[nodiscard]] std::optional<std::variant<std::int64_t, double>> max(bool keep_nans,
                                                                      bool keep_infs,
                                                                      bool keep_zeros) const;
  [[nodiscard]] double mean(bool keep_nans, bool keep_infs, bool keep_zeros) const;
  [[nodiscard]] double variance(bool keep_nans, bool keep_infs, bool keep_zeros, bool exact) const;
  [[nodiscard]] double skewness(bool keep_nans, bool keep_infs, bool keep_zeros, bool exact) const;
  [[nodiscard]] double kurtosis(bool keep_nans, bool keep_infs, bool keep_zeros, bool exact) const;

  [[nodiscard]] static auto parse_span(std::string_view span) -> QuerySpan;
  [[nodiscard]] static std::string_view count_type_to_str(const PixelVar& var);

  [[nodiscard]] CoolerGlobalLock::UniqueLock lock() const noexcept;

  static void bind(nanobind::module_& m);

 private:
  // NOLINTBEGIN(bugprone-exception-escape)
  [[nodiscard]] hictk::PixelCoordinates coord1() const noexcept;
  [[nodiscard]] hictk::PixelCoordinates coord2() const noexcept;

  [[nodiscard]] const hictk::BinTable& bins() const noexcept;
  // NOLINTEND(bugprone-exception-escape)
};

}  // namespace hictkpy
