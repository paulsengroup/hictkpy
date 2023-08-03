// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <memory>
#include <tuple>
#include <variant>

#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/hic/pixel_selector.hpp"
#include "hictk/transformers/join_genomic_coords.hpp"

namespace hictkpy {

struct PixelSelector {
  // clang-format off
  using SelectorVar =
    std::variant<std::shared_ptr<const hictk::cooler::PixelSelector>,
                 std::shared_ptr<const hictk::hic::PixelSelector>,
                 std::shared_ptr<const hictk::hic::PixelSelectorAll>>;

  using PixelVar = std::variant<std::int32_t, double>;
  // clang-format on

  SelectorVar selector{};
  PixelVar pixel_count{std::int32_t(0)};
  bool join{};

  PixelSelector() = default;

  PixelSelector(std::shared_ptr<const hictk::cooler::PixelSelector> sel_, std::string_view type,
                bool join_);
  PixelSelector(std::shared_ptr<const hictk::hic::PixelSelector> sel_, std::string_view type,
                bool join_);
  PixelSelector(std::shared_ptr<const hictk::hic::PixelSelectorAll> sel_, std::string_view type,
                bool join_);

  using PixelCoordTuple =
      std::tuple<std::string, std::int32_t, std::int32_t, std::string, std::int32_t, std::int32_t>;

  [[nodiscard]] auto get_coord1() const -> PixelCoordTuple;
  [[nodiscard]] auto get_coord2() const -> PixelCoordTuple;

  [[nodiscard]] pybind11::iterator make_iterable() const;
  [[nodiscard]] pybind11::object to_df() const;
  [[nodiscard]] pybind11::object to_coo() const;
  [[nodiscard]] pybind11::object to_numpy() const;
  [[nodiscard]] pybind11::object sum() const;
  [[nodiscard]] std::int64_t nnz() const;

 private:
  [[nodiscard]] constexpr bool int_pixels() const noexcept;
  [[nodiscard]] constexpr bool float_pixels() const noexcept;

  [[nodiscard]] hictk::PixelCoordinates coord1() const noexcept;
  [[nodiscard]] hictk::PixelCoordinates coord2() const noexcept;

  [[nodiscard]] const hictk::BinTable& bins() const noexcept;
};

}  // namespace hictkpy
