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
                                               bool keep_zeros, bool exact) {
  validate_metrics(metrics);
  reset();

  if (exact) {
    const auto stats =
        compute_online<keep_nans, keep_infs>(first, last, size, {"nnz", "mean"}, keep_zeros);
    assert(stats.nnz.has_value());
    reset();
    return compute_exact<keep_nans, keep_infs>(
        std::move(first), std::move(last), size, metrics, static_cast<std::uint64_t>(*stats.nnz),
        stats.mean.value_or(std::numeric_limits<double>::quiet_NaN()), keep_zeros);
  }

  return compute_online<keep_nans, keep_infs>(std::move(first), std::move(last), size, metrics,
                                              keep_zeros);
}

template <typename PixelIt>
template <bool keep_nans, bool keep_infs>
inline Stats PixelAggregator<PixelIt>::compute_online(
    PixelIt first, PixelIt last, std::uint64_t size,
    const phmap::flat_hash_set<std::string>& metrics, bool keep_zeros) {
  auto break_on_non_finite = [this]() constexpr noexcept {
    return _nan_found || _neg_inf_found || _pos_inf_found;
  };
  auto break_on_nans = [this]() constexpr noexcept { return _nan_found; };

  process_pixels<keep_nans, keep_infs>(first, last, break_on_non_finite);
  if (first != last) {
    if (metrics.contains("nnz")) {
      // if we need to compute the nnz, then we have no tricks up our sleeves
      process_pixels<keep_nans, keep_infs>(first, last, []() constexpr noexcept { return false; });
    } else {
      process_pixels<keep_nans, keep_infs>(first, last, break_on_nans);
    }
  }

  if (keep_zeros) {
    update_with_zeros(size - _nnz);
  }

  return extract(metrics);
}

template <typename PixelIt>
template <bool keep_nans, bool keep_infs>
inline Stats PixelAggregator<PixelIt>::compute_exact(
    PixelIt first, PixelIt last, std::uint64_t size,
    const phmap::flat_hash_set<std::string>& metrics, std::uint64_t nnz, double mean,
    bool keep_zeros) {
  if (size == 0) {
    return {};
  }

  if (size == nnz) {
    keep_zeros = false;
  }

  _nnz = nnz;
  _num_zeros = keep_zeros ? size - nnz : 0;
  _online_mean = mean;

  const auto count_fp = static_cast<double>(count());
  double variance = 0;
  const auto variance_denom = count() < 2 ? std::numeric_limits<double>::quiet_NaN() : count_fp - 1;

  while (first != last) {
    const auto n = first->count;
    std::ignore = ++first;
    if (internal::drop_value<keep_nans, keep_infs>(n)) {
      continue;
    }
    update_finiteness_counters<keep_nans, keep_infs>(n);
    _min = std::min(_min, conditional_static_cast<CountT>(n));
    _max = std::max(_max, conditional_static_cast<CountT>(n));
    _sum += conditional_static_cast<CountT>(n);
    const auto delta = conditional_static_cast<double>(n) - mean;
    variance += delta * delta / variance_denom;
    _online_m2 += delta * delta / count_fp;
    _online_m3 += delta * delta * delta / count_fp;
    _online_m4 += delta * delta * delta * delta / count_fp;
  }

  if (keep_zeros) {
    _min = std::min(_min, CountT{0});
    _max = std::max(_max, CountT{0});
    _online_mean = conditional_static_cast<double>(_sum) / count_fp;
    const auto delta = -mean;

    variance += static_cast<double>(size - nnz) * delta * delta / variance_denom;
    _online_m2 += static_cast<double>(size - nnz) * delta * delta / count_fp;
    _online_m3 += static_cast<double>(size - nnz) * delta * delta * delta / count_fp;
    _online_m4 += static_cast<double>(size - nnz) * delta * delta * delta * delta / count_fp;
  }

  auto stats = extract(metrics, true);
  if (metrics.contains("variance") && count() > 1) {
    stats.variance = variance;
  }

  return stats;
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
Stats PixelAggregator<PixelIt>::extract(const phmap::flat_hash_set<std::string>& metrics,
                                        bool no_divide) {
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
    stats.variance = compute_variance(no_divide);
  }

  if (metrics.contains("skewness")) {
    stats.skewness = compute_skewness(no_divide);
  }

  if (metrics.contains("kurtosis")) {
    stats.kurtosis = compute_kurtosis(no_divide);
  }

  return stats;
}

template <typename PixelIt>
template <bool keep_nans, bool keep_infs>
inline void PixelAggregator<PixelIt>::update(N n) noexcept {
  update_finiteness_counters<keep_nans, keep_infs>(n);
  if (n != 0) {
    ++_nnz;
  } else {
    ++_num_zeros;
  }
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
  //
  // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
  // clang-format on

  const auto n_fp = conditional_static_cast<double>(n);
  const auto count_fp = static_cast<double>(count());

  const auto delta = n_fp - _online_mean;
  const auto delta_scaled = delta / static_cast<double>(count());
  const auto delta_scaled_squared = delta_scaled * delta_scaled;
  const auto term1 = delta * delta_scaled * static_cast<double>(count() - 1);

  // NOLINTBEGIN(*-math-missing-parentheses, *-avoid-magic-numbers)
  _online_mean += delta_scaled;
  _online_m4 += term1 * delta_scaled_squared * ((count_fp * count_fp) - (3 * count_fp) + 3) +
                6 * delta_scaled_squared * _online_m2 - 4 * delta_scaled * _online_m3;
  _online_m3 += term1 * delta_scaled * (count_fp - 2) - 3 * delta_scaled * _online_m2;
  _online_m2 += term1;
  // NOLINTEND(*-math-missing-parentheses, *-avoid-magic-numbers)
}

template <typename PixelIt>
inline void PixelAggregator<PixelIt>::update_with_zeros(std::uint64_t num_zeros) noexcept {
  if (num_zeros == 0) {
    return;
  }

  if (_nan_found || _neg_inf_found || _pos_inf_found) {
    _min = std::min(_min, CountT{0});
    _max = std::max(_max, CountT{0});
    _num_zeros += num_zeros;
    return;
  }

  // TODO optimize
  for (std::uint64_t i = 0; i < num_zeros; ++i) {
    update<true, true>(0);
  }
}

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
inline std::optional<double> PixelAggregator<PixelIt>::compute_variance(
    bool no_divide) const noexcept {
  if (count() < 2) {
    return {};
  }
  if (_nan_found || _pos_inf_found || _neg_inf_found) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (no_divide) {
    return _online_m2;
  }
  return _online_m2 / static_cast<double>(count() - 1);
}

template <typename PixelIt>
inline std::optional<double> PixelAggregator<PixelIt>::compute_skewness(
    bool no_divide) const noexcept {
  if (count() < 2) {
    return {};
  }
  if (_nan_found || _pos_inf_found || _neg_inf_found) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const auto m2 = no_divide ? _online_m2 : _online_m2 / static_cast<double>(count());
  const auto m3 = no_divide ? _online_m3 : _online_m3 / static_cast<double>(count());
  return m3 / std::pow(m2, 1.5);  // NOLINT(*-avoid-magic-numbers)
}

template <typename PixelIt>
inline std::optional<double> PixelAggregator<PixelIt>::compute_kurtosis(
    bool no_divide) const noexcept {
  if (count() < 2) {
    return {};
  }
  if (_nan_found || _pos_inf_found || _neg_inf_found) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const auto m2 = no_divide ? _online_m2 : _online_m2 / static_cast<double>(count());
  const auto m4 = no_divide ? _online_m4 : _online_m4 / static_cast<double>(count());
  return (m4 / (m2 * m2)) - 3.0;  // NOLINT(*-avoid-magic-numbers)
}

}  // namespace hictkpy
