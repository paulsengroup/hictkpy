// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/cooler/pixel_selector.hpp>
#include <hictk/hic/pixel_selector.hpp>
#include <hictk/numeric_variant.hpp>
#include <hictk/pixel.hpp>
#include <hictk/transformers/common.hpp>
#include <hictk/transformers/to_dataframe.hpp>
#include <memory>
#include <string>
#include <string_view>
#include <tuple>
#include <variant>
#include <vector>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

struct PixelSelector {
  // clang-format off
  using SelectorVar =
    std::variant<std::shared_ptr<const hictk::cooler::PixelSelector>,
                 std::shared_ptr<const hictk::hic::PixelSelector>,
                 std::shared_ptr<const hictk::hic::PixelSelectorAll>>;
  // clang-format on
  using PixelVar = hictk::internal::NumericVariant;
  using QuerySpan = hictk::transformers::QuerySpan;
  using PixelFormat = hictk::transformers::DataFrameFormat;

  SelectorVar selector{};
  PixelVar pixel_count{std::int32_t{0}};
  hictk::transformers::DataFrameFormat pixel_format{hictk::transformers::DataFrameFormat::COO};

  PixelSelector() = default;

  PixelSelector(std::shared_ptr<const hictk::cooler::PixelSelector> sel_, std::string_view type,
                bool join);
  PixelSelector(std::shared_ptr<const hictk::hic::PixelSelector> sel_, std::string_view type,
                bool join);
  PixelSelector(std::shared_ptr<const hictk::hic::PixelSelectorAll> sel_, std::string_view type,
                bool join);

  [[nodiscard]] std::string repr() const;

  using PixelCoordTuple =
      std::tuple<std::string, std::int32_t, std::int32_t, std::string, std::int32_t, std::int32_t>;

  [[nodiscard]] auto get_coord1() const -> PixelCoordTuple;
  [[nodiscard]] auto get_coord2() const -> PixelCoordTuple;

  [[nodiscard]] nanobind::iterator make_iterable() const;
  [[nodiscard]] nanobind::object to_arrow(std::string_view span) const;
  [[nodiscard]] nanobind::object to_df(std::string_view span) const;
  [[nodiscard]] nanobind::object to_pandas(std::string_view span) const;
  [[nodiscard]] nanobind::object to_coo(std::string_view span) const;
  [[nodiscard]] nanobind::object to_csr(std::string_view span) const;
  [[nodiscard]] nanobind::object to_numpy(std::string_view span) const;

  [[nodiscard]] nanobind::dict describe(const std::vector<std::string>& metrics, bool keep_nans,
                                        bool keep_infs) const;
  [[nodiscard]] std::int64_t nnz(bool keep_nans, bool keep_infs) const;
  [[nodiscard]] nanobind::object sum(bool keep_nans, bool keep_infs) const;
  [[nodiscard]] nanobind::object min(bool keep_nans = false, bool keep_infs = false) const;
  [[nodiscard]] nanobind::object max(bool keep_nans = false, bool keep_infs = false) const;
  [[nodiscard]] double mean(bool keep_nans = false, bool keep_infs = false) const;
  [[nodiscard]] double variance(bool keep_nans = false, bool keep_infs = false) const;
  [[nodiscard]] double skewness(bool keep_nans = false, bool keep_infs = false) const;
  [[nodiscard]] double kurtosis(bool keep_nans = false, bool keep_infs = false) const;

  [[nodiscard]] static auto parse_span(std::string_view span) -> QuerySpan;
  [[nodiscard]] static auto parse_count_type(std::string_view type) -> PixelVar;
  [[nodiscard]] static std::string_view count_type_to_str(const PixelVar& var);

  static void bind(nanobind::module_& m);

 private:
  // NOLINTBEGIN(bugprone-exception-escape)
  [[nodiscard]] hictk::PixelCoordinates coord1() const noexcept;
  [[nodiscard]] hictk::PixelCoordinates coord2() const noexcept;

  [[nodiscard]] const hictk::BinTable& bins() const noexcept;
  // NOLINTEND(bugprone-exception-escape)
};

}  // namespace hictkpy
