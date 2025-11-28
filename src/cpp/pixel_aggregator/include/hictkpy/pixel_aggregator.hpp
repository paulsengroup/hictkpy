// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/phmap.h>

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>

#include "hictkpy/common.hpp"

namespace hictkpy {

struct Stats {
  std::optional<std::int64_t> nnz{};
  std::optional<std::variant<std::int64_t, double>> sum{};
  std::optional<std::variant<std::int64_t, double>> min{};
  std::optional<std::variant<std::int64_t, double>> max{};
  std::optional<double> mean{};
  std::optional<double> variance{};
  std::optional<double> skewness{};
  std::optional<double> kurtosis{};
};

template <typename PixelIt>
class PixelAggregator {
  using N = remove_cvref_t<decltype(std::declval<PixelIt>()->count)>;
  static_assert(std::is_arithmetic_v<N>);
  using CountT = std::conditional_t<std::is_floating_point_v<N>, double, std::int64_t>;
  std::uint64_t _nnz{};
  std::uint64_t _num_zeros{};
  CountT _min{};
  CountT _max{};
  CountT _sum{};
  double _online_mean{};
  double _online_m2{};
  double _online_m3{};
  double _online_m4{};

  bool _finite_found{false};
  bool _nan_found{false};
  bool _neg_inf_found{false};
  bool _pos_inf_found{false};

 public:
  static constexpr std::array<std::string_view, 8> valid_metrics{
      "nnz", "sum", "min", "max", "mean", "variance", "skewness", "kurtosis"};

  PixelAggregator() = default;

  // NOLINTBEGIN(*-unnecessary-value-param)
  template <bool keep_nans, bool keep_infs>
  Stats compute(PixelIt first, PixelIt last, std::uint64_t size,
                const phmap::flat_hash_set<std::string>& metrics, bool keep_zeros, bool exact);

 private:
  template <bool keep_nans, bool keep_infs>
  Stats compute_online(PixelIt first, PixelIt last, std::uint64_t size,
                       const phmap::flat_hash_set<std::string>& metrics, bool keep_zeros);

  template <bool keep_nans, bool keep_infs>
  Stats compute_exact(PixelIt first, PixelIt last, std::uint64_t size,
                      const phmap::flat_hash_set<std::string>& metrics, std::uint64_t nnz,
                      double mean, bool keep_zeros);
  // NOLINTEND(*-unnecessary-value-param)
  static void validate_metrics(const phmap::flat_hash_set<std::string>& metrics);
  template <bool keep_nans, bool keep_infs>
  void update_finiteness_counters(N n) noexcept;
  template <bool keep_nans, bool keep_infs, typename StopCondition>
  void process_pixels(PixelIt& first, const PixelIt& last, StopCondition break_condition);
  [[nodiscard]] Stats extract(const phmap::flat_hash_set<std::string>& metrics,
                              bool no_divide = false);
  void reset() noexcept;

  template <bool keep_nans, bool keep_infs>
  void update(N n) noexcept;
  void update_with_zeros(std::uint64_t num_zeros) noexcept;

  [[nodiscard]] std::uint64_t count() const noexcept;
  [[nodiscard]] auto compute_min() const noexcept -> std::optional<CountT>;
  [[nodiscard]] auto compute_max() const noexcept -> std::optional<CountT>;
  [[nodiscard]] std::optional<double> compute_mean() const noexcept;
  [[nodiscard]] std::optional<double> compute_variance(bool no_divide) const noexcept;
  [[nodiscard]] std::optional<double> compute_skewness(bool no_divide) const noexcept;
  [[nodiscard]] std::optional<double> compute_kurtosis(bool no_divide) const noexcept;
};

}  // namespace hictkpy

#include "../../pixel_aggregator_impl.hpp"
