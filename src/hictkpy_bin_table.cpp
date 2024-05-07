// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <algorithm>
#include <hictkpy/common.hpp>
#include <string>
#include <string_view>
#include <vector>

#include "hictkpy/bin_table.hpp"

namespace nb = nanobind;

namespace hictkpy {

BinTable::BinTable(const hictk::BinTable& bins) : _bins(bins) {}

BinTable::BinTable(nb::dict chromosomes, std::uint32_t resolution) {
  std::vector<std::uint32_t> chrom_sizes{};

  std::for_each(chromosomes.begin(), chromosomes.end(), [&](const auto& kv) {
    _chrom_names.push_back(nb::cast<std::string>(kv.first));
    chrom_sizes.push_back(nb::cast<std::uint32_t>(kv.second));
  });

  _bins = hictk::BinTable(
      hictk::Reference{_chrom_names.begin(), _chrom_names.end(), chrom_sizes.begin()}, resolution);
}

const hictk::Reference& BinTable::chromosomes() const noexcept { return _bins.chromosomes(); }

std::uint32_t BinTable::resolution() const noexcept { return _bins.resolution(); }

std::size_t BinTable::size() const noexcept { return _bins.size(); }

std::string BinTable::repr() const {
  return fmt::format(FMT_STRING("BinTable(num_chroms={}; bin_size={};)"),
                     _bins.chromosomes().size(), _bins.resolution());
}

nb::object BinTable::bin_ids_to_coords(BinIDsT bin_ids) const {
  auto pd = nb::module_::import_("pandas");

  Dynamic1DA<std::uint32_t> chrom_ids(bin_ids.size());
  Dynamic1DA<std::int32_t> starts(bin_ids.size());
  Dynamic1DA<std::int32_t> ends(bin_ids.size());

  std::visit(
      [&](const auto& bins) {
        for (std::size_t i = 0; i < bin_ids.size(); ++i) {
          const auto bin = bins.at(bin_ids(i));
          chrom_ids.push_back(bin.id());
          starts.push_back(static_cast<std::int32_t>(bin.start()));
          ends.push_back(static_cast<std::int32_t>(bin.end()));
        }
      },
      _bins.get());

  nb::dict py_bins_dict{};  // NOLINT

  py_bins_dict["chrom"] =
      pd.attr("Categorical")
          .attr("from_codes")(chrom_ids(), "categories"_a = _chrom_names, "validate"_a = false);

  py_bins_dict["chrom"] = chrom_ids();
  py_bins_dict["start"] = starts();
  py_bins_dict["end"] = ends();

  return pd.attr("DataFrame")(py_bins_dict, "copy"_a = false);
}

hictk::Bin BinTable::bin_id_to_coord(std::uint64_t bin_id) const { return _bins.at(bin_id); }

nb::object BinTable::coord_to_bin(std::string_view chrom, std::uint32_t pos) const {
  return nb::cast(_bins.at(chrom, pos));
}

nanobind::object BinTable::coords_to_bins(const std::vector<std::string>& chroms,
                                          const std::vector<std::uint32_t>& positions) const {
  auto pd = nb::module_::import_("pandas");

  if (chroms.size() != positions.size()) {
    throw std::runtime_error("chroms and positions should have the same size");
  }

  Dynamic1DA<std::uint32_t> chrom_ids(chroms.size());
  Dynamic1DA<std::int32_t> starts(chroms.size());
  Dynamic1DA<std::int32_t> ends(chroms.size());

  std::visit(
      [&](const auto& bins) {
        for (std::size_t i = 0; i < chroms.size(); ++i) {
          const auto bin = bins.at(chroms[i], positions[i]);
          chrom_ids.push_back(bin.id());
          starts.push_back(static_cast<std::int32_t>(bin.start()));
          ends.push_back(static_cast<std::int32_t>(bin.end()));
        }
      },
      _bins.get());

  nb::dict py_bins_dict{};  // NOLINT

  py_bins_dict["chrom"] =
      pd.attr("Categorical")
          .attr("from_codes")(chrom_ids(), "categories"_a = _chrom_names, "validate"_a = false);

  py_bins_dict["chrom"] = chrom_ids();
  py_bins_dict["start"] = starts();
  py_bins_dict["end"] = ends();

  return pd.attr("DataFrame")(py_bins_dict, "copy"_a = false);
}

std::int64_t BinTable::coord_to_bin_id(std::string_view chrom, std::uint32_t pos) const {
  return static_cast<std::int64_t>(_bins.at(chrom, pos).id());
}

auto BinTable::coords_to_bin_ids(const std::vector<std::string>& chroms,
                                 const std::vector<std::uint32_t>& positions) const -> BinIDsT {
  if (chroms.size() != positions.size()) {
    throw std::runtime_error("chroms and positions should have the same size");
  }

  Dynamic1DA<std::int64_t> bin_ids{chroms.size()};

  std::visit(
      [&](const auto& bins) {
        for (std::size_t i = 0; i < chroms.size(); ++i) {
          bin_ids.push_back(static_cast<std::int64_t>(bins.at(chroms[i], positions[i]).id()));
        }
      },
      _bins.get());

  return bin_ids();
}

nb::object BinTable::merge_coords(nb::object df) const {
  auto pd = nb::module_::import_("pandas");

  using Buffer64T = nb::ndarray<nb::numpy, nb::shape<nb::any>, std::int64_t>;
  auto bin1_ids = nb::cast<Buffer64T>(df.attr("__getitem__")("bin1_id").attr("to_numpy")());
  auto bin2_ids = nb::cast<Buffer64T>(df.attr("__getitem__")("bin2_id").attr("to_numpy")());

  const auto n = bin1_ids.size();

  Dynamic1DA<std::uint32_t> chrom1_ids(n);
  Dynamic1DA<std::int32_t> starts1(n);
  Dynamic1DA<std::int32_t> ends1(n);
  Dynamic1DA<std::uint32_t> chrom2_ids(n);
  Dynamic1DA<std::int32_t> starts2(n);
  Dynamic1DA<std::int32_t> ends2(n);

  std::visit(
      [&](const auto& bins) {
        for (std::size_t i = 0; i < n; ++i) {
          const auto bin1 = bins.at(bin1_ids(i));
          chrom1_ids.push_back(bin1.chrom().id());
          starts1.push_back(static_cast<std::int32_t>(bin1.start()));
          ends1.push_back(static_cast<std::int32_t>(bin1.end()));

          const auto bin2 = bins.at(bin2_ids(i));
          chrom2_ids.push_back(bin2.chrom().id());
          starts2.push_back(static_cast<std::int32_t>(bin2.start()));
          ends2.push_back(static_cast<std::int32_t>(bin2.end()));
        }
      },
      _bins.get());

  df = df.attr("copy")();

  df.attr("__setitem__")(
      "chrom1",
      pd.attr("Categorical")
          .attr("from_codes")(chrom1_ids(), "categories"_a = _chrom_names, "validate"_a = false));
  df.attr("__setitem__")("start1", starts1());
  df.attr("__setitem__")("end1", ends1());

  df.attr("__setitem__")(
      "chrom2",
      pd.attr("Categorical")
          .attr("from_codes")(chrom2_ids(), "categories"_a = _chrom_names, "validate"_a = false));
  df.attr("__setitem__")("start2", starts1());
  df.attr("__setitem__")("end2", ends1());

  return df;
}

nb::iterator BinTable::make_iterable() const {
  return std::visit(
      [](const auto& bins) {
        return nb::make_iterator(nb::type<hictk::BinTable>(), "BinTableIterator", bins.begin(),
                                 bins.end());
      },
      _bins.get());
}

nb::object BinTable::to_df() const {
  auto pd = nb::module_::import_("pandas");

  using Buffer64T = nb::ndarray<nb::numpy, nb::shape<nb::any>, std::int64_t>;

  const auto n = _bins.size();

  Dynamic1DA<std::uint32_t> chrom_ids(n);
  Dynamic1DA<std::int32_t> starts(n);
  Dynamic1DA<std::int32_t> ends(n);

  std::visit(
      [&](const auto& bins) {
        for (const auto& bin : bins) {
          chrom_ids.push_back(bin.chrom().id());
          starts.push_back(static_cast<std::int32_t>(bin.start()));
          ends.push_back(static_cast<std::int32_t>(bin.end()));
        }
      },
      _bins.get());

  nb::dict py_bins_dict{};  // NOLINT

  py_bins_dict["chrom"] =
      pd.attr("Categorical")
          .attr("from_codes")(chrom_ids(), "categories"_a = _chrom_names, "validate"_a = false);
  py_bins_dict["start"] = starts();
  py_bins_dict["end"] = ends();

  auto df = pd.attr("DataFrame")(py_bins_dict, "copy"_a = false);
  return df;
}

}  // namespace hictkpy
