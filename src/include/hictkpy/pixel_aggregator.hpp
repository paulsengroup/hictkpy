// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/phmap.h>

#include <array>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

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

class PixelAggregator {
  template <typename N>
  using Accumulator = boost::accumulators::accumulator_set<
      N, boost::accumulators::stats<
             boost::accumulators::droppable<boost::accumulators::tag::count>,
             boost::accumulators::droppable<boost::accumulators::tag::sum>,
             boost::accumulators::droppable<boost::accumulators::tag::min>,
             boost::accumulators::droppable<boost::accumulators::tag::max>,
             boost::accumulators::droppable<boost::accumulators::tag::mean>,
             boost::accumulators::droppable<boost::accumulators::tag::variance>,
             boost::accumulators::droppable<boost::accumulators::tag::skewness>>>;
  // When using ints there is a high chance that computing the 3rd or 4th moments causes an int
  // overflow, leading to completely incorrect results.
  // At the same time we cannot use doubles everywhere, because for large numbers the value of
  // e.g. sum, mean etc. will no longer be exact
  using KurtosisAccumulator = boost::accumulators::accumulator_set<
      double, boost::accumulators::stats<boost::accumulators::tag::kurtosis>>;

  std::variant<Accumulator<std::int64_t>, Accumulator<double>> _accumulator{
      Accumulator<std::int64_t>{}};
  std::optional<KurtosisAccumulator> _kurtosis_accumulator{};

  bool _compute_count{true};
  bool _compute_sum{true};
  bool _compute_min{true};
  bool _compute_max{true};
  bool _compute_mean{true};
  bool _compute_variance{true};
  bool _compute_skewness{true};
  bool _compute_kurtosis{true};

  bool _finite_found{false};
  bool _nan_found{false};
  bool _neg_inf_found{false};
  bool _pos_inf_found{false};

 public:
  static constexpr std::array<std::string_view, 8> valid_metrics{
      "nnz", "sum", "min", "max", "mean", "variance", "skewness", "kurtosis"};

  PixelAggregator() = default;

  template <typename N, bool keep_nans, bool keep_infs, typename PixelSelector>
  [[nodiscard]] Stats compute(const PixelSelector& sel,
                              const phmap::flat_hash_set<std::string>& metrics, bool exact);

 private:
  static void validate_metrics(const phmap::flat_hash_set<std::string>& metrics);

  template <bool keep_nans, bool keep_infs, typename N>
  void update_finiteness_counters(N n) noexcept;
  template <bool keep_nans, bool keep_infs, bool skip_kurtosis, typename N, typename It,
            typename StopCondition>
  [[nodiscard]] std::pair<It, std::size_t> process_pixels_until_true(Accumulator<N>& accumulator,
                                                                     It first, It last,
                                                                     StopCondition break_loop);

  template <typename N>
  [[nodiscard]] bool sum_can_be_skipped(Accumulator<N>& accumulator) const;
  template <typename N>
  [[nodiscard]] bool min_can_be_skipped(Accumulator<N>& accumulator) const;
  template <typename N>
  [[nodiscard]] bool max_can_be_skipped(Accumulator<N>& accumulator) const;
  template <typename N>
  [[nodiscard]] bool mean_can_be_skipped(Accumulator<N>& accumulator) const;
  template <typename N>
  void disable_redundant_accumulators(Accumulator<N>& accumulator);

  template <bool keep_nans, bool keep_infs, bool skip_kurtosis, typename N, typename It>
  void process_all_remaining_pixels(Accumulator<N>& accumulator, It&& first, It&& last);

  template <typename N>
  void reset(const phmap::flat_hash_set<std::string>& metrics);

  [[nodiscard]] Stats extract(const phmap::flat_hash_set<std::string>& metrics);
  template <typename N, std::size_t keep_nans, std::size_t keep_infs, typename It>
  [[nodiscard]] std::optional<Stats> handle_edge_cases(
      It first, It last, const phmap::flat_hash_set<std::string>& metrics);

  template <typename N>
  [[nodiscard]] static std::int64_t extract_nnz(const Accumulator<N>& accumulator);
  template <typename N>
  [[nodiscard]] N extract_sum(const Accumulator<N>& accumulator) const;
  template <typename N>
  [[nodiscard]] N extract_min(const Accumulator<N>& accumulator) const;
  template <typename N>
  [[nodiscard]] N extract_max(const Accumulator<N>& accumulator) const;
  template <typename N>
  [[nodiscard]] double extract_mean(const Accumulator<N>& accumulator) const;
  [[nodiscard]] double extract_mean(const KurtosisAccumulator& accumulator) const;
  template <typename N>
  [[nodiscard]] double extract_variance(const Accumulator<N>& accumulator) const;
  template <typename N>
  [[nodiscard]] double extract_skewness(const Accumulator<N>& accumulator) const;
  [[nodiscard]] double extract_kurtosis() const;
};

}  // namespace hictkpy

#include "../../pixel_aggregator_impl.hpp"
