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

namespace hictkpy::bin_table {

void ctor(hictk::BinTable* bins, nanobind::dict chromosomes, std::uint32_t resolution) {
  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> chrom_sizes{};

  std::for_each(chromosomes.begin(), chromosomes.end(), [&](const auto& kv) {
    chrom_names.push_back(nb::cast<std::string>(kv.first));
    chrom_sizes.push_back(nb::cast<std::uint32_t>(kv.second));
  });

  new (bins) hictk::BinTable(
      hictk::Reference{chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()}, resolution);
}

std::string repr(const hictk::BinTable& bins) {
  return fmt::format(FMT_STRING("BinTable(num_chroms={}; bin_size={};)"), bins.chromosomes().size(),
                     bins.resolution());
}

hictk::Bin bin_id_to_coords(const hictk::BinTable& bins, std::uint64_t bin_id) {
  return bins.at(bin_id);
}

nb::object try_bin_id_to_coords(const hictk::BinTable& bins, std::uint64_t bin_id) {
  try {
    return nb::cast(bins.at(bin_id));
  } catch (const std::out_of_range&) {
    return nb::none();
  }
}

nb::object try_coords_to_bin(const hictk::BinTable& bins, std::string_view chrom,
                                std::uint32_t pos) {
  try {
    return nb::cast(bins.at(chrom, pos));
  } catch (const std::out_of_range&) {
    return nb::none();
  }
}

nb::object merge_coords(const hictk::BinTable& bins, nb::object df) {
  auto pd = nb::module_::import_("pandas");

  using Buffer64T = nb::ndarray<nb::numpy, nb::shape<nb::any>, std::int64_t>;
  auto bin1_ids = nb::cast<Buffer64T>(df.attr("__getitem__")("bin1_id").attr("to_numpy")());
  auto bin2_ids = nb::cast<Buffer64T>(df.attr("__getitem__")("bin2_id").attr("to_numpy")());

  const auto n = bin1_ids.size();

  Dynamic1DA<uint32_t> chrom1_ids(n);
  Dynamic1DA<std::int32_t> starts1(n);
  Dynamic1DA<std::int32_t> ends1(n);
  Dynamic1DA<uint32_t> chrom2_ids(n);
  Dynamic1DA<std::int32_t> starts2(n);
  Dynamic1DA<std::int32_t> ends2(n);

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

  std::vector<std::string> chrom_names{};
  std::transform(bins.chromosomes().begin(), bins.chromosomes().end(),
                 std::back_inserter(chrom_names),
                 [&](const hictk::Chromosome& chrom) { return std::string{chrom.name()}; });

  df = df.attr("copy")();

  df.attr("__setitem__")(
      "chrom1",
      pd.attr("Categorical")
          .attr("from_codes")(chrom1_ids(), "categories"_a = chrom_names, "validate"_a = false));
  df.attr("__setitem__")("start1", pd.attr("Series")(starts1(), "copy"_a = false));
  df.attr("__setitem__")("end1", pd.attr("Series")(ends1(), "copy"_a = false));

  df.attr("__setitem__")(
      "chrom2",
      pd.attr("Categorical")
          .attr("from_codes")(chrom2_ids(), "categories"_a = chrom_names, "validate"_a = false));
  df.attr("__setitem__")("start2", pd.attr("Series")(starts1(), "copy"_a = false));
  df.attr("__setitem__")("end2", pd.attr("Series")(ends1(), "copy"_a = false));

  return df;
}

nb::iterator make_iterable(const hictk::BinTable& bins) {
  return nb::make_iterator(nb::type<hictk::BinTable>(), "BinTableIterator", bins.begin(),
                           bins.end());
}

nb::object to_df(const hictk::BinTable& bins) {
  auto pd = nb::module_::import_("pandas");

  using Buffer64T = nb::ndarray<nb::numpy, nb::shape<nb::any>, std::int64_t>;

  const auto n = bins.size();

  Dynamic1DA<uint32_t> chrom_ids(n);
  Dynamic1DA<std::int32_t> starts(n);
  Dynamic1DA<std::int32_t> ends(n);

  for (const auto& bin : bins) {
    chrom_ids.push_back(bin.chrom().id());
    starts.push_back(static_cast<std::int32_t>(bin.start()));
    ends.push_back(static_cast<std::int32_t>(bin.end()));
  }

  std::vector<std::string> chrom_names{};
  std::transform(bins.chromosomes().begin(), bins.chromosomes().end(),
                 std::back_inserter(chrom_names),
                 [&](const hictk::Chromosome& chrom) { return std::string{chrom.name()}; });

  nb::dict py_bins_dict{};  // NOLINT

  py_bins_dict["chrom"] =
      pd.attr("Categorical")
          .attr("from_codes")(chrom_ids(), "categories"_a = chrom_names, "validate"_a = false);
  py_bins_dict["start"] = pd.attr("Series")(starts(), "copy"_a = false);
  py_bins_dict["end"] = pd.attr("Series")(ends(), "copy"_a = false);

  auto df = pd.attr("DataFrame")(py_bins_dict, "copy"_a = false);
  return df;
}

}  // namespace hictkpy::bin_table
