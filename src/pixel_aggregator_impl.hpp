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

template <bool keep_nans, bool keep_infs, typename PixelIt>
[[nodiscard]] inline auto compute_descriptors_exact(PixelIt first, PixelIt last, std::size_t nnz,
                                                    double mean) {
  assert(first != last);

  struct Result {
    double variance{};
    double m2{};
    double m3{};
    double m4{};
  };

  if (std::isnan(mean)) {
    return Result{
        std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  }

  Result res{};

  const auto nnz_fp = static_cast<double>(nnz);
  const auto variance_denom = nnz < 2 ? std::numeric_limits<double>::quiet_NaN() : nnz_fp - 1;

  while (first != last) {
    const auto n = conditional_static_cast<double>(first->count);
    std::ignore = ++first;
    if (internal::drop_value<keep_nans, keep_infs>(n)) {
      continue;
    }
    const auto delta = n - mean;
    res.variance += delta * delta / variance_denom;
    res.m2 += delta * delta;
    res.m3 += delta * delta * delta;
    res.m4 += delta * delta * delta * delta;
  }

  res.m2 /= nnz_fp;
  res.m3 /= nnz_fp;
  res.m4 /= nnz_fp;

  if (nnz < 2) {
    assert(std::isnan(res.variance));
  }

  return res;
}

[[nodiscard]] inline double skewness(std::size_t size, double mean, double m2, double m3) {
  if (size < 2 || std::isnan(mean) || m2 == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return m3 / std::pow(m2, 1.5);  // NOLINT(*-avoid-magic-numbers)
}

[[nodiscard]] inline double kurtosis(std::size_t size, double mean, double m2, double m4) {
  if (size < 2 || std::isnan(mean) || m2 == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return (m4 / (m2 * m2)) - 3.0;  // NOLINT(*-avoid-magic-numbers)
}

}  // namespace internal

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

  auto nnz = std::visit(
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
          if (_compute_kurtosis) {
            constexpr auto skip_kurtosis = false;
            process_all_remaining_pixels<keep_nans, keep_infs, skip_kurtosis>(
                accumulator, std::move(first), std::move(last));
          } else {
            constexpr auto skip_kurtosis = true;
            process_all_remaining_pixels<keep_nans, keep_infs, skip_kurtosis>(
                accumulator, std::move(first), std::move(last));
          }
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

  auto stats = extract(metrics);
  if (_neg_inf_found || _pos_inf_found || _nan_found) {
    return stats;
  }

  if (!exact && nnz.value_or(10'000) >= 10'000) {  // NOLINT(*-avoid-magic-numbers)
    // exact computation is not required and sample size is big enough
    return stats;
  }

  if (!stats.variance.has_value() && !stats.skewness.has_value() && !stats.kurtosis.has_value()) {
    // all requested metrics have already been computed exactly
    return stats;
  }

  const auto mean =
      _kurtosis_accumulator.has_value()
          ? extract_mean(*_kurtosis_accumulator)
          : std::visit([&](const auto& accumulator) { return extract_mean(accumulator); },
                       _accumulator);

  nnz = static_cast<std::size_t>(
      _kurtosis_accumulator.has_value()
          ? extract_nnz(*_kurtosis_accumulator)
          : std::visit([&](const auto& accumulator) { return extract_nnz(accumulator); },
                       _accumulator));

  const auto [variance, m2, m3, m4] = internal::compute_descriptors_exact<keep_nans, keep_infs>(
      sel.template begin<N>(), sel.template end<N>(), *nnz, mean);

  if (stats.variance.has_value()) {
    stats.variance = variance;
  }

  if (stats.skewness.has_value()) {
    stats.skewness = internal::skewness(*nnz, mean, m2, m3);
  }

  if (stats.kurtosis.has_value()) {
    stats.kurtosis = internal::kurtosis(*nnz, mean, m2, m4);
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
inline void PixelAggregator::update_finiteness_counters(N n) noexcept {
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

template <bool keep_nans, bool keep_infs, bool skip_kurtosis, typename N, typename It,
          typename StopCondition>
inline std::pair<It, std::size_t> PixelAggregator::process_pixels_until_true(
    Accumulator<N>& accumulator, It first, It last, StopCondition break_loop) {
  // This is just a workaround to allow wrapping drop_value and early_return with HICTKPY_UNLIKELY
  auto drop_pixel = [](const auto n) noexcept {
    return internal::drop_value<keep_nans, keep_infs>(n);
  };

  std::size_t nnz = 0;
  while (first != last) {
    const auto n = first->count;
    if (HICTKPY_UNLIKELY(drop_pixel(n))) {
      std::ignore = ++first;
      continue;
    }
    update_finiteness_counters<keep_nans, keep_infs>(n);
    if (HICTKPY_UNLIKELY(break_loop())) {
      break;
    }
    accumulator(n);
    ++nnz;

    if constexpr (!skip_kurtosis) {
      assert(_kurtosis_accumulator.has_value());
      (*_kurtosis_accumulator)(conditional_static_cast<double>(n));
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
}

template <bool keep_nans, bool keep_infs, bool skip_kurtosis, typename N, typename It>
inline void PixelAggregator::process_all_remaining_pixels(Accumulator<N>& accumulator, It&& first,
                                                          It&& last) {
  assert(_compute_count);
  // This is just a workaround to allow wrapping drop_value and early_return with HICTKPY_UNLIKELY
  auto drop_pixel = [](const auto n) noexcept {
    return internal::drop_value<keep_nans, keep_infs>(n);
  };

  std::for_each(std::forward<It>(first), std::forward<It>(last), [&](const auto& pixel) {
    if (HICTKPY_UNLIKELY(drop_pixel(pixel.count))) {
      return;
    }
    update_finiteness_counters<keep_nans, keep_infs>(pixel.count);
    accumulator(pixel.count);
    if constexpr (!skip_kurtosis) {
      assert(_kurtosis_accumulator);
      (*_kurtosis_accumulator)(conditional_static_cast<double>(pixel.count));
    }
  });
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

  _finite_found = false;
  _nan_found = false;
  _neg_inf_found = false;
  _pos_inf_found = false;

  if constexpr (std::is_floating_point_v<N>) {
    _accumulator = Accumulator<double>{};
  } else {
    _accumulator = Accumulator<std::int64_t>{};
  }

  if (metrics.contains("kurtosis")) {
    _kurtosis_accumulator = KurtosisAccumulator{};
  } else {
    _kurtosis_accumulator.reset();
    _compute_kurtosis = false;
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
      },
      _accumulator);
}

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
      },
      _accumulator);

  if (metrics.contains("kurtosis")) {
    stats.kurtosis = extract_kurtosis();
  }

  return stats;
}

template <typename N, std::size_t keep_nans, std::size_t keep_infs, typename It>
inline std::optional<Stats> PixelAggregator::handle_edge_cases(
    It first, It last, const phmap::flat_hash_set<std::string>& metrics) {
  // Handle selectors with zero or one (valid) pixel
  std::size_t nnz = 0;
  N value{};
  while (first != last && nnz < 2) {
    if (!internal::drop_value<keep_nans, keep_infs>(first->count)) {
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
    update_finiteness_counters<keep_nans, keep_infs>(value);
    std::visit([&](auto& accumulator) { accumulator(value); }, _accumulator);
    auto stats = extract(metrics);
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

inline std::int64_t PixelAggregator::extract_nnz(const KurtosisAccumulator& accumulator) {
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
    if (!_finite_found && _pos_inf_found) {
      return std::numeric_limits<N>::infinity();
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
    if (!_finite_found && _neg_inf_found) {
      return -std::numeric_limits<N>::infinity();
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

inline double PixelAggregator::extract_mean(const KurtosisAccumulator& accumulator) const {
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

inline double PixelAggregator::extract_kurtosis() const {
  if (_nan_found || _pos_inf_found || _neg_inf_found) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  assert(_kurtosis_accumulator.has_value());
  return boost::accumulators::kurtosis(*_kurtosis_accumulator);
}

}  // namespace hictkpy
