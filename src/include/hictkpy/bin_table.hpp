// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <cstdint>
#include <hictk/bin.hpp>
#include <hictk/bin_table.hpp>
#include <string_view>

#include "common.hpp"

namespace hictkpy {

class BinTable {
  hictk::BinTable _bins{};
  std::vector<std::string> _chrom_names{};

 public:
  using BinIDsT = nanobind::ndarray<nanobind::numpy, nanobind::shape<nanobind::any>, std::int64_t>;
  explicit BinTable(const hictk::BinTable& bins);
  BinTable(nanobind::dict chromosomes, std::uint32_t resolution);
  // BinTable(nanobind::object bins);

  [[nodiscard]] const hictk::Reference& chromosomes() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;

  [[nodiscard]] std::string repr() const;

  [[nodiscard]] hictk::Bin bin_id_to_coord(std::uint64_t bin_id) const;
  [[nodiscard]] nanobind::object bin_ids_to_coords(BinIDsT bin_ids) const;

  [[nodiscard]] nanobind::object coord_to_bin(std::string_view chrom, std::uint32_t pos) const;
  nanobind::object coords_to_bins(const std::vector<std::string>& chroms,
                                  const std::vector<std::uint32_t>& positions) const;

  [[nodiscard]] std::int64_t coord_to_bin_id(std::string_view chrom, std::uint32_t pos) const;
  auto coords_to_bin_ids(const std::vector<std::string>& chroms,
                         const std::vector<std::uint32_t>& positions) const -> BinIDsT;

  [[nodiscard]] nanobind::object merge_coords(nanobind::object df) const;

  [[nodiscard]] nanobind::iterator make_iterable() const;

  [[nodiscard]] nanobind::object to_df() const;
};

}  // namespace hictkpy
