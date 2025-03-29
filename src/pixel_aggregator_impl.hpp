// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "hictkpy/common.hpp"

namespace hictkpy {
namespace internal {
template <bool keep_nans, bool keep_infs, typename N>
[[nodiscard]] inline bool drop_value([[maybe_unused]] N n) noexcept {
  static_assert(std::is_arithmetic_v<N>);
  if constexpr (!std::is_floating_point_v<N>) {
    return false;
  } else {
    // MSVC gets confused if this chunk of code is not in an else branch...
    if constexpr (!keep_nans && !keep_infs) {
      return !std::isfinite(n);
    }
    if constexpr (!keep_nans) {
      return std::isnan(n);
    }
    if constexpr (!keep_infs) {
      return std::isinf(n);
    }
  }

  return false;
}

}  // namespace internal

template <typename PixelIt>
template <bool keep_nans, bool keep_infs>
inline Stats PixelAggregator<PixelIt>::compute(PixelIt first, PixelIt last, std::uint64_t size,
                                               const phmap::flat_hash_set<std::string>& metrics,
                                               bool keep_zeros, [[maybe_unused]] bool exact) {
  validate_metrics(metrics);
  reset();

  auto break_on_non_finite = [this]() constexpr noexcept {
    return _nan_found || _neg_inf_found || _pos_inf_found;
  };
  auto break_on_nans = [this]() constexpr noexcept { return _nan_found; };

  process_pixels<keep_nans, keep_infs>(first, last, break_on_non_finite);
  if (first == last) {
    if (keep_zeros) {
      update_with_zeros(size - _nnz);
    }
    return extract(metrics);
  }

  if (metrics.contains("nnz")) {
    // if we need to compute the nnz, then we have no tricks up our sleeves
    process_pixels<keep_nans, keep_infs>(first, last, []() constexpr noexcept { return false; });
  } else {
    process_pixels<keep_nans, keep_infs>(first, last, break_on_nans);
  }

  if (keep_zeros) {
    update_with_zeros(size - _nnz);
  }

  return extract(metrics);
  // TODO deal with exact computations
}

template <typename PixelIt>
inline void PixelAggregator<PixelIt>::validate_metrics(
    const phmap::flat_hash_set<std::string>& metrics) {
  for (const auto& metric : metrics) {
    if (std::find(valid_metrics.begin(), valid_metrics.end(), metric) == valid_metrics.end()) {
      throw std::out_of_range(
          fmt::format(FMT_STRING("unknown metric \"{}\". Valid metrics are: {}"), metric,
                      fmt::join(valid_metrics, ", ")));
    }
  }
}

template <typename PixelIt>
template <bool keep_nans, bool keep_infs>
inline void PixelAggregator<PixelIt>::update_finiteness_counters(N n) noexcept {
  static_assert(std::is_arithmetic_v<N>);
  if constexpr (std::is_floating_point_v<N>) {
    if constexpr (keep_nans) {
      if (HICTKPY_UNLIKELY(std::isnan(n))) {
        _nan_found = true;
        return;
      }
    }

    if constexpr (keep_infs) {
      if (HICTKPY_UNLIKELY(std::isinf(n))) {
        if (n > 0) {
          _pos_inf_found = true;
        } else {
          _neg_inf_found = true;
        }
        return;
      }
    }
  }
  _finite_found = true;
}

template <typename PixelIt>
template <bool keep_nans, bool keep_infs, typename StopCondition>
inline void PixelAggregator<PixelIt>::process_pixels(PixelIt& first, const PixelIt& last,
                                                     StopCondition break_condition) {
  // This is just a workaround to allow wrapping drop_value and early_return with HICTKPY_LIKELY
  auto drop_pixel = [](const auto n) noexcept {
    return internal::drop_value<keep_nans, keep_infs>(n);
  };

  for (; first != last && !break_condition(); ++first) {
    if (HICTKPY_LIKELY(!drop_pixel(first->count))) {
      update<keep_nans, keep_infs>(first->count);
    }
  }
}

template <typename PixelIt>
Stats PixelAggregator<PixelIt>::extract(const phmap::flat_hash_set<std::string>& metrics) {
  assert(!metrics.empty());

  Stats stats{};

  if (metrics.contains("nnz")) {
    stats.nnz = static_cast<std::int64_t>(_nnz);
  }

  if (metrics.contains("sum")) {
    stats.sum = _sum;
  }

  if (metrics.contains("min")) {
    stats.min = compute_min();
  }

  if (metrics.contains("max")) {
    stats.max = compute_max();
  }

  if (metrics.contains("mean")) {
    stats.mean = compute_mean();
  }

  if (metrics.contains("variance")) {
    stats.variance = compute_variance();
  }

  if (metrics.contains("skewness")) {
    stats.skewness = compute_skewness();
  }

  if (metrics.contains("kurtosis")) {
    stats.kurtosis = compute_kurtosis();
  }

  return stats;
}

template <typename PixelIt>
template <bool keep_nans, bool keep_infs>
inline void PixelAggregator<PixelIt>::update(N n) noexcept {
  assert(n != 0);
  update_finiteness_counters<keep_nans, keep_infs>(n);
  ++_nnz;
  _min = std::min(_min, conditional_static_cast<CountT>(n));
  _max = std::max(_max, conditional_static_cast<CountT>(n));
  _sum += conditional_static_cast<CountT>(n);

  // clang-format off
  // the formulas to compute the mean and 2nd, 3rd, and 4th moments are taken from the following two
  // publications:
  // Chan TF, Golub GH, LeVeque RJ. Updating formulae and a pairwise algorithm for computing sample variances.
  // In COMPSTAT 1982 5th Symposium held at Toulouse 1982
  // https://link.springer.com/chapter/10.1007/978-3-642-51461-6_3
  // Terriberry TB. Computing higher-order moments online [Internet]. 2008 Jul 13
  // https://web.archive.org/web/20140423031833/http://people.xiph.org/~tterribe/notes/homs.html  // codespell:ignore
  // clang-format on

  const auto n_fp = conditional_static_cast<double>(n);
  const auto count_fp = static_cast<double>(_nnz);

  const auto delta = n_fp - _online_mean;
  const auto delta_scaled = delta / static_cast<double>(_nnz);
  const auto delta_scaled_squared = delta_scaled * delta_scaled;
  const auto term1 = delta * delta_scaled * static_cast<double>(_nnz - 1);

  _online_mean += delta_scaled;
  _online_m4 += term1 * delta_scaled_squared * ((count_fp * count_fp) - (3 * count_fp) + 3) +
                6 * delta_scaled_squared * _online_m2 - 4 * delta_scaled * _online_m3;
  _online_m3 += term1 * delta_scaled * (count_fp - 2) - 3 * delta_scaled * _online_m2;
  _online_m2 += term1;
}

template <typename PixelIt>
inline void PixelAggregator<PixelIt>::update_with_zeros(
    [[maybe_unused]] std::uint64_t num_zeros) noexcept {}

template <typename PixelIt>
inline void PixelAggregator<PixelIt>::reset() noexcept {
  _nnz = 0;
  _num_zeros = 0;
  _min = std::numeric_limits<CountT>::max();
  _max = std::numeric_limits<CountT>::lowest();
  _sum = 0;
  _online_mean = 0;
  _online_m2 = 0;
  _online_m3 = 0;
  _online_m4 = 0;

  _finite_found = false;
  _nan_found = false;
  _neg_inf_found = false;
  _pos_inf_found = false;
}

template <typename PixelIt>
inline std::uint64_t PixelAggregator<PixelIt>::count() const noexcept {
  return _num_zeros + _nnz;
}

template <typename PixelIt>
inline auto PixelAggregator<PixelIt>::compute_min() const noexcept -> std::optional<CountT> {
  if constexpr (std::is_floating_point_v<CountT>) {
    if (_nan_found) {
      return std::numeric_limits<CountT>::quiet_NaN();
    }
    if (_neg_inf_found) {
      return -std::numeric_limits<CountT>::infinity();
    }
    if (!_finite_found && _pos_inf_found) {
      return std::numeric_limits<CountT>::infinity();
    }
  }
  if (count() == 0) {
    return {};
  }
  return _min;
}

template <typename PixelIt>
inline auto PixelAggregator<PixelIt>::compute_max() const noexcept -> std::optional<CountT> {
  if constexpr (std::is_floating_point_v<CountT>) {
    if (_nan_found) {
      return std::numeric_limits<CountT>::quiet_NaN();
    }
    if (_pos_inf_found) {
      return std::numeric_limits<CountT>::infinity();
    }
    if (!_finite_found && _neg_inf_found) {
      return -std::numeric_limits<CountT>::infinity();
    }
  }
  if (count() == 0) {
    return {};
  }
  return _max;
}

template <typename PixelIt>
inline std::optional<double> PixelAggregator<PixelIt>::compute_mean() const noexcept {
  if (count() == 0) {
    return {};
  }
  if (_nan_found || (_neg_inf_found && _pos_inf_found)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return conditional_static_cast<double>(_sum) / static_cast<double>(count());
}

template <typename PixelIt>
inline std::optional<double> PixelAggregator<PixelIt>::compute_variance() const noexcept {
  if (count() < 2) {
    return {};
  }
  if (_nan_found || _pos_inf_found || _neg_inf_found) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return _online_m2 / static_cast<double>(count() - 1);
}

template <typename PixelIt>
inline std::optional<double> PixelAggregator<PixelIt>::compute_skewness() const noexcept {
  if (count() < 2) {
    return {};
  }
  if (_nan_found || _pos_inf_found || _neg_inf_found) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const auto m2 = _online_m2 / static_cast<double>(count());
  const auto m3 = _online_m3 / static_cast<double>(count());
  return m3 / std::pow(m2, 1.5);  // NOLINT(*-avoid-magic-numbers)
}

template <typename PixelIt>
inline std::optional<double> PixelAggregator<PixelIt>::compute_kurtosis() const noexcept {
  if (count() < 2) {
    return {};
  }
  if (_nan_found || _pos_inf_found || _neg_inf_found) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const auto m2 = _online_m2 / static_cast<double>(count());
  const auto m4 = _online_m4 / static_cast<double>(count());
  return (m4 / (m2 * m2)) - 3.0;  // NOLINT(*-avoid-magic-numbers)
}

}  // namespace hictkpy
