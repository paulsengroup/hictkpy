// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>

#include "hictkpy/common.hpp"

namespace hictkpy {

template <typename It>
[[nodiscard]] inline double compute_variance_exact(It first, It last, double mean,
                                                   std::size_t size) {
  assert(size != 0);
  if (HICTKPY_UNLIKELY(size == 1)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return std::accumulate(
      std::move(first), std::move(last), 0.0, [&](const auto accumulator, const auto& pixel) {
        const auto n = conditional_static_cast<double>(pixel.count);
        return accumulator + (((n - mean) * (n - mean)) / static_cast<double>(size - 1));
      });
}

template <typename N, bool keep_nans, bool keep_infs, typename PixelSelector>
inline Stats PixelAggregator::compute(const PixelSelector& sel,
                                      const phmap::flat_hash_set<std::string>& metrics,
                                      bool exact) {
  static_assert(std::is_same_v<N, std::int64_t> || std::is_same_v<N, double>);

  validate_metrics(metrics);

  reset<N>(metrics);
  auto first = sel.template begin<N>();
  auto last = sel.template end<N>();

  if (auto stats = handle_edge_cases<N, keep_nans, keep_infs>(first, last, metrics);
      stats.has_value()) {
    return *stats;
  }

  auto break_on_non_finite = [this]() constexpr noexcept {
    return _nan_found || _neg_inf_found || _pos_inf_found;
  };
  auto break_on_nans = [this]() constexpr noexcept { return _nan_found; };

  const auto nnz = std::visit(
      [&](auto& accumulator) -> std::optional<std::size_t> {
        std::size_t nnz_ = 0;
        if (_compute_kurtosis) {
          constexpr bool skip_kurtosis = false;
          auto result = process_pixels_until_true<keep_nans, keep_infs, skip_kurtosis>(
              accumulator, std::move(first), last, break_on_non_finite);
          first = std::move(result.first);
          nnz_ = result.second;
        } else {
          constexpr bool skip_kurtosis = true;
          auto result = process_pixels_until_true<keep_nans, keep_infs, skip_kurtosis>(
              accumulator, std::move(first), last, break_on_non_finite);
          first = std::move(result.first);
          nnz_ = result.second;
        }
        if (first == last) {
          return nnz_;
        }

        const auto skip_sum = sum_can_be_skipped(accumulator);
        const auto skip_min = min_can_be_skipped(accumulator);
        const auto skip_max = max_can_be_skipped(accumulator);
        const auto skip_mean = mean_can_be_skipped(accumulator);
        disable_redundant_accumulators(accumulator);

        if (_compute_count) {
          // if we need to compute the nnz, then we have no tricks up our sleeves
          process_all_remaining_pixels<keep_nans, keep_infs>(accumulator, std::move(first),
                                                             std::move(last));
          return {};
        }

        if (skip_sum && skip_min && skip_max && skip_mean) {
          return {};
        }

        constexpr bool skip_kurtosis = true;
        std::ignore = process_pixels_until_true<keep_nans, keep_infs, skip_kurtosis>(
            accumulator, std::move(first), std::move(last), break_on_nans);

        return {};
      },
      _accumulator);

  auto stats = extract<N>(metrics);
  const auto variance_may_be_inaccurate = nnz.value_or(10'000) < 10'000;
  if ((variance_may_be_inaccurate || exact) && stats.variance.has_value()) {
    const auto mean = std::visit([&](const auto& accumulator) { return extract_mean(accumulator); },
                                 _accumulator);
    stats.variance =
        compute_variance_exact(sel.template begin<N>(), sel.template end<N>(), mean, *nnz);
  }
  return stats;
}

inline void PixelAggregator::validate_metrics(const phmap::flat_hash_set<std::string>& metrics) {
  for (const auto& metric : metrics) {
    if (std::find(valid_metrics.begin(), valid_metrics.end(), metric) == valid_metrics.end()) {
      throw std::out_of_range(
          fmt::format(FMT_STRING("unknown metric \"{}\". Valid metrics are: {}"), metric,
                      fmt::join(valid_metrics, ", ")));
    }
  }
}

template <bool keep_nans, bool keep_infs, typename N>
inline void PixelAggregator::process_non_finite(N n) noexcept {
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
      }
    }
  }
}

template <bool keep_nans, bool keep_infs, bool skip_kurtosis, typename N, typename It,
          typename StopCondition>
inline std::pair<It, std::size_t> PixelAggregator::process_pixels_until_true(
    Accumulator<N>& accumulator, It first, It last, StopCondition break_loop) {
  // This is just a workaround to allow wrapping drop_value and early_return with HICTKPY_UNLIKELY
  auto drop_pixel = [this](const auto n) noexcept { return drop_value<keep_nans, keep_infs>(n); };

  std::size_t nnz = 0;
  while (first != last) {
    const auto n = first->count;
    if (HICTKPY_UNLIKELY(drop_pixel(n))) {
      std::ignore = ++first;
      continue;
    }
    process_non_finite<keep_nans, keep_infs>(n);
    if (HICTKPY_UNLIKELY(break_loop())) {
      break;
    }
    accumulator(n);
    ++nnz;

    if constexpr (!skip_kurtosis && std::is_integral_v<N>) {
      assert(_kurtosis_accumulator.has_value());
      (*_kurtosis_accumulator)(static_cast<double>(n));
    }
    std::ignore = ++first;
  }
  return std::make_pair(first, nnz);
}

template <typename N>
inline bool PixelAggregator::sum_can_be_skipped(Accumulator<N>& accumulator) const {
  if (_compute_sum && (_nan_found || (_neg_inf_found && _pos_inf_found))) {
    accumulator.template drop<boost::accumulators::tag::sum>();
    return true;
  }
  return false;
}

template <typename N>
inline bool PixelAggregator::min_can_be_skipped(Accumulator<N>& accumulator) const {
  if (_compute_min && _nan_found) {
    accumulator.template drop<boost::accumulators::tag::min>();
    return true;
  }
  return false;
}

template <typename N>
inline bool PixelAggregator::max_can_be_skipped(Accumulator<N>& accumulator) const {
  if (_compute_max && _nan_found) {
    accumulator.template drop<boost::accumulators::tag::max>();
    return true;
  }
  return false;
}

template <typename N>
inline bool PixelAggregator::mean_can_be_skipped(Accumulator<N>& accumulator) const {
  if (_compute_mean && (_nan_found || (_neg_inf_found && _pos_inf_found))) {
    accumulator.template drop<boost::accumulators::tag::mean>();
    return true;
  }
  return false;
}

template <typename N>
inline void PixelAggregator::disable_redundant_accumulators(Accumulator<N>& accumulator) {
  if (_compute_variance) {
    accumulator.template drop<boost::accumulators::tag::variance>();
  }
  if (_compute_skewness) {
    accumulator.template drop<boost::accumulators::tag::skewness>();
  }
  if (_compute_kurtosis) {
    accumulator.template drop<boost::accumulators::tag::kurtosis>();
    _kurtosis_accumulator.reset();
  }
}

template <bool keep_nans, bool keep_infs, typename N, typename It>
inline void PixelAggregator::process_all_remaining_pixels(Accumulator<N>& accumulator, It&& first,
                                                          It&& last) {
  assert(_compute_count);
  // This is just a workaround to allow wrapping drop_value and early_return with HICTKPY_UNLIKELY
  auto drop_pixel = [this](const auto n) noexcept { return drop_value<keep_nans, keep_infs>(n); };

  std::for_each(std::forward<It>(first), std::forward<It>(last), [&](const auto& pixel) {
    if (HICTKPY_UNLIKELY(drop_pixel(pixel.count))) {
      return;
    }
    process_non_finite<keep_nans, keep_infs>(pixel.count);
    accumulator(pixel.count);
  });
}

template <bool keep_nans, bool keep_infs, typename N>
inline bool PixelAggregator::drop_value([[maybe_unused]] N n) noexcept {
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

template <typename N>
inline void PixelAggregator::reset(const phmap::flat_hash_set<std::string>& metrics) {
  if (metrics.empty()) {
    throw std::runtime_error("provide one or more statistics to be computed");
  }

  _compute_count = true;
  _compute_sum = true;
  _compute_min = true;
  _compute_max = true;
  _compute_mean = true;
  _compute_variance = true;
  _compute_skewness = true;
  _compute_kurtosis = true;

  _nan_found = false;
  _neg_inf_found = false;
  _pos_inf_found = false;

  if constexpr (std::is_floating_point_v<N>) {
    _accumulator = Accumulator<double>{};
    _kurtosis_accumulator.reset();
  } else {
    _accumulator = Accumulator<std::int64_t>{};
    if (metrics.contains("kurtosis")) {
      _kurtosis_accumulator = KurtosisAccumulator{};
    } else {
      _kurtosis_accumulator.reset();
      _compute_kurtosis = false;
    }
  }

  std::visit(
      [&](auto& accumulator) {
        if (!metrics.contains("nnz")) {
          accumulator.template drop<boost::accumulators::tag::count>();
          _compute_count = false;
        }

        if (!metrics.contains("sum")) {
          accumulator.template drop<boost::accumulators::tag::sum>();
          _compute_sum = false;
        }

        if (!metrics.contains("min")) {
          accumulator.template drop<boost::accumulators::tag::min>();
          _compute_min = false;
        }

        if (!metrics.contains("max")) {
          accumulator.template drop<boost::accumulators::tag::max>();
          _compute_max = false;
        }

        if (!metrics.contains("mean")) {
          accumulator.template drop<boost::accumulators::tag::mean>();
          _compute_mean = false;
        }

        if (!metrics.contains("variance")) {
          accumulator.template drop<boost::accumulators::tag::variance>();
          _compute_variance = false;
        }

        if (!metrics.contains("skewness")) {
          accumulator.template drop<boost::accumulators::tag::skewness>();
          _compute_skewness = false;
        }

        if (!metrics.contains("kurtosis") || std::is_integral_v<N>) {
          accumulator.template drop<boost::accumulators::tag::kurtosis>();
          _compute_kurtosis = metrics.contains("kurtosis");
        }
      },
      _accumulator);
}

template <typename N>
inline Stats PixelAggregator::extract(const phmap::flat_hash_set<std::string>& metrics) {
  assert(!metrics.empty());

  Stats stats{};

  std::visit(
      [&](const auto& accumulator) {
        if (metrics.contains("nnz")) {
          stats.nnz = extract_nnz(accumulator);
        }

        if (metrics.contains("sum")) {
          stats.sum = extract_sum(accumulator);
        }

        if (metrics.contains("min")) {
          stats.min = extract_min(accumulator);
        }

        if (metrics.contains("max")) {
          stats.max = extract_max(accumulator);
        }

        if (metrics.contains("mean")) {
          stats.mean = extract_mean(accumulator);
        }

        if (metrics.contains("variance")) {
          stats.variance = extract_variance(accumulator);
        }

        if (metrics.contains("skewness")) {
          stats.skewness = extract_skewness(accumulator);
        }

        if (metrics.contains("kurtosis")) {
          stats.kurtosis = extract_kurtosis(accumulator);
        }
      },
      _accumulator);

  return stats;
}

template <typename N, std::size_t keep_nans, std::size_t keep_infs, typename It>
inline std::optional<Stats> PixelAggregator::handle_edge_cases(
    It first, It last, const phmap::flat_hash_set<std::string>& metrics) {
  // Handle selectors with zero or one (valid) pixels
  std::size_t nnz = 0;
  N value{};
  while (first != last && nnz < 2) {
    if (!drop_value<keep_nans, keep_infs>(first->count)) {
      ++nnz;
      value = first->count;
    }
    ++first;
  }

  if (HICTKPY_UNLIKELY(nnz == 0)) {
    Stats stats{};
    if (metrics.contains("nnz")) {
      stats.nnz = 0;
    }
    if (metrics.contains("sum")) {
      stats.sum = N{0};
    }
    return stats;
  }

  if (HICTKPY_UNLIKELY(nnz == 1)) {
    std::visit([&](auto& accumulator) { accumulator(value); }, _accumulator);
    auto stats = extract<N>(metrics);
    if (stats.variance.has_value()) {
      stats.variance = std::numeric_limits<double>::quiet_NaN();
    }
    if (stats.skewness.has_value()) {
      stats.skewness = std::numeric_limits<double>::quiet_NaN();
    }
    if (stats.kurtosis.has_value()) {
      stats.kurtosis = std::numeric_limits<double>::quiet_NaN();
    }

    return stats;
  }

  return {};
}

template <typename N>
inline std::int64_t PixelAggregator::extract_nnz(const Accumulator<N>& accumulator) {
  return static_cast<std::int64_t>(boost::accumulators::count(accumulator));
}

template <typename N>
inline N PixelAggregator::extract_sum(const Accumulator<N>& accumulator) const {
  if constexpr (std::is_floating_point_v<N>) {
    if (_nan_found || (_neg_inf_found && _pos_inf_found)) {
      return std::numeric_limits<N>::quiet_NaN();
    }
    if (_pos_inf_found) {
      return std::numeric_limits<N>::infinity();
    }
    if (_neg_inf_found) {
      return -std::numeric_limits<N>::infinity();
    }
  }
  return boost::accumulators::sum(accumulator);
}

template <typename N>
inline N PixelAggregator::extract_min(const Accumulator<N>& accumulator) const {
  if constexpr (std::is_floating_point_v<N>) {
    if (_nan_found) {
      return std::numeric_limits<N>::quiet_NaN();
    }
    if (_neg_inf_found) {
      return -std::numeric_limits<N>::infinity();
    }
  }
  return boost::accumulators::min(accumulator);
}

template <typename N>
inline N PixelAggregator::extract_max(const Accumulator<N>& accumulator) const {
  if constexpr (std::is_floating_point_v<N>) {
    if (_nan_found) {
      return std::numeric_limits<N>::quiet_NaN();
    }
    if (_pos_inf_found) {
      return std::numeric_limits<N>::infinity();
    }
  }
  return boost::accumulators::max(accumulator);
}

template <typename N>
inline double PixelAggregator::extract_mean(const Accumulator<N>& accumulator) const {
  if (_nan_found || (_neg_inf_found && _pos_inf_found)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (_pos_inf_found) {
    return std::numeric_limits<double>::infinity();
  }
  if (_neg_inf_found) {
    return -std::numeric_limits<double>::infinity();
  }
  return boost::accumulators::mean(accumulator);
}

template <typename N>
inline double PixelAggregator::extract_variance(const Accumulator<N>& accumulator) const {
  if (_nan_found || _pos_inf_found || _neg_inf_found) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return boost::accumulators::variance(accumulator);
}

template <typename N>
inline double PixelAggregator::extract_skewness(const Accumulator<N>& accumulator) const {
  if (_nan_found || _pos_inf_found || _neg_inf_found) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return boost::accumulators::skewness(accumulator);
}

template <typename N>
inline double PixelAggregator::extract_kurtosis(const Accumulator<N>& accumulator) const {
  if (_nan_found || _pos_inf_found || _neg_inf_found) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (_kurtosis_accumulator.has_value()) {
    return boost::accumulators::kurtosis(*_kurtosis_accumulator);
  }
  return boost::accumulators::kurtosis(accumulator);
}

}  // namespace hictkpy
