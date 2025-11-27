// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <hictk/bin.hpp>
#include <hictk/pixel.hpp>
#include <type_traits>
#include <utility>

#include "hictkpy/common.hpp"

namespace hictkpy {
template <typename N>
inline Pixel::Pixel(hictk::Pixel<N> p) noexcept
    : _coords(std::move(p.coords)),
      _bin1_id(static_cast<std::int64_t>(_coords->bin1.id())),
      _bin2_id(static_cast<std::int64_t>(_coords->bin2.id())),
      _count(cast_count(p.count)) {}

template <typename N>
inline Pixel::Pixel(const hictk::ThinPixel<N> &p) noexcept
    : Pixel(static_cast<std::int64_t>(p.bin1_id), static_cast<std::int64_t>(p.bin2_id), p.count) {}

template <typename N>
inline Pixel::Pixel(hictk::Bin bin1, hictk::Bin bin2, N count) noexcept
    : Pixel(hictk::Pixel{{std::move(bin1), std::move(bin2)}, count}) {}

template <typename N>
inline Pixel::Pixel(std::int64_t bin1_id, std::int64_t bin2_id, N count_) noexcept
    : _bin1_id(bin1_id), _bin2_id(bin2_id), _count(cast_count(count_)) {}

template <typename N>
inline Pixel &Pixel::operator=(hictk::Pixel<N> p) noexcept {
  _coords = std::move(p.coords);
  _bin1_id = static_cast<std::int64_t>(_coords->bin1.id());
  _bin2_id = static_cast<std::int64_t>(_coords->bin2.id());
  _count = cast_count(p.count);

  return *this;
}

template <typename N>
inline Pixel &Pixel::operator=(const hictk::ThinPixel<N> &p) noexcept {
  _coords = hictk::PixelCoordinates{};
  _bin1_id = static_cast<std::int64_t>(p.bin1_id);
  _bin2_id = static_cast<std::int64_t>(p.bin2_id);
  _count = cast_count(p.count);

  return *this;
}

constexpr std::int64_t Pixel::bin1_id() const noexcept { return _bin1_id; }
constexpr std::int64_t Pixel::bin2_id() const noexcept { return _bin2_id; }

template <typename N>
constexpr std::variant<std::int64_t, double> Pixel::cast_count(N n) noexcept {
  static_assert(std::is_arithmetic_v<N>);
  if constexpr (std::is_floating_point_v<N>) {
    return {conditional_static_cast<double>(n)};
  } else {
    return {conditional_static_cast<std::int64_t>(n)};
  }
}

}  // namespace hictkpy
