// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <hictk/bin.hpp>
#include <hictk/pixel.hpp>
#include <optional>
#include <string>
#include <string_view>
#include <variant>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

class Pixel {
  std::optional<hictk::PixelCoordinates> _coords{};
  std::int64_t _bin1_id{};
  std::int64_t _bin2_id{};
  std::variant<std::int64_t, double> _count{std::int64_t{}};

 public:
  Pixel() = default;
  template <typename N>
  Pixel(hictk::Pixel<N> p) noexcept;  // NOLINT(*-explicit-conversions)

  template <typename N>
  Pixel(const hictk::ThinPixel<N> &p) noexcept;  // NOLINT(*-explicit-conversions)

  template <typename N>
  Pixel(hictk::Bin bin1, hictk::Bin bin2, N count) noexcept;

  template <typename N>
  Pixel(std::int64_t bin1_id, std::int64_t bin2_id, N count_) noexcept;

  template <typename N>
  Pixel &operator=(hictk::Pixel<N> p) noexcept;

  template <typename N>
  Pixel &operator=(const hictk::ThinPixel<N> &p) noexcept;

  [[nodiscard]] constexpr std::int64_t bin1_id() const noexcept;
  [[nodiscard]] constexpr std::int64_t bin2_id() const noexcept;
  [[nodiscard]] std::variant<std::int64_t, double> count() const noexcept;

  [[nodiscard]] const hictk::PixelCoordinates &coords() const;

  [[nodiscard]] const hictk::Bin &bin1() const;
  [[nodiscard]] const hictk::Bin &bin2() const;

  [[nodiscard]] std::string_view chrom1() const;
  [[nodiscard]] std::int64_t start1() const;
  [[nodiscard]] std::int64_t end1() const;
  [[nodiscard]] std::string_view chrom2() const;
  [[nodiscard]] std::int64_t start2() const;
  [[nodiscard]] std::int64_t end2() const;

  [[nodiscard]] std::string repr() const;
  [[nodiscard]] std::string str() const;

  static void bind(nanobind::module_ &m);

 private:
  template <typename N>
  [[nodiscard]] static constexpr std::variant<std::int64_t, double> cast_count(N n) noexcept;
};

}  // namespace hictkpy

#include "../../pixel_impl.hpp"
