// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <hictk/bin.hpp>
#include <hictk/bin_table.hpp>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "hictkpy/nanobind.hpp"
#include "hictkpy/reference.hpp"

namespace hictkpy {

class BinTable {
  std::shared_ptr<const hictk::BinTable> _bins = std::make_shared<const hictk::BinTable>();

 public:
  using BinIDsT = nanobind::ndarray<nanobind::numpy, nanobind::shape<-1>, std::int64_t>;
  BinTable() = default;
  explicit BinTable(std::shared_ptr<const hictk::BinTable> bins) noexcept;
  explicit BinTable(hictk::BinTable bins);
  BinTable(const ChromosomeDict& chromosomes, std::uint32_t resolution);
  // BinTable(nanobind::object bins);

  [[nodiscard]] const hictk::Reference& chromosomes() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;
  [[nodiscard]] std::string_view type() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;

  [[nodiscard]] std::string repr() const;

  [[nodiscard]] hictk::Bin bin_id_to_coord(std::uint64_t bin_id) const;
  [[nodiscard]] nanobind::object bin_ids_to_coords(std::vector<std::uint64_t> bin_ids) const;

  [[nodiscard]] hictk::Bin coord_to_bin(std::string_view chrom, std::uint32_t pos) const;
  [[nodiscard]] nanobind::object coords_to_bins(const std::vector<std::string>& chroms,
                                                const std::vector<std::uint32_t>& positions) const;

  [[nodiscard]] std::int64_t coord_to_bin_id(std::string_view chrom, std::uint32_t pos) const;
  [[nodiscard]] auto coords_to_bin_ids(const std::vector<std::string>& chroms,
                                       const std::vector<std::uint32_t>& positions) const
      -> BinIDsT;

  [[nodiscard]] nanobind::object merge_coords(nanobind::object df) const;

  [[nodiscard]] nanobind::iterator make_iterable() const;

  [[nodiscard]] nanobind::object to_df(std::optional<std::string_view> range,
                                       std::string_view query_type) const;

  [[nodiscard]] std::shared_ptr<const hictk::BinTable> get() const noexcept;

  static void bind(nanobind::module_& m);

 private:
  [[nodiscard]] std::vector<std::string> chrom_names(bool include_ALL = false) const;
};

// TODO remove after merging https://github.com/paulsengroup/hictk/pull/298
template <typename T>
inline hictkpy::BinTable get_bins_from_sc_file(const T& obj) {
  return hictkpy::BinTable(obj.bins());
}

template <typename T>
inline hictkpy::BinTable get_bins_from_object(const T& obj) {
  return hictkpy::BinTable(obj.bins_ptr());
}

}  // namespace hictkpy
