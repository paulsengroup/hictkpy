// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstdint>
#include <hictk/chromosome.hpp>
#include <iterator>
#include <string>
#include <vector>

#include "hictkpy/dynamic_1d_array.hpp"
#include "hictkpy/nanobind.hpp"

namespace hictkpy {

template <typename File>
inline nanobind::object get_bins_from_file(const File &f) {
  auto pd = nanobind::module_::import_("pandas");

  Dynamic1DA<std::int32_t> chrom_ids{};
  Dynamic1DA<std::int32_t> starts{};
  Dynamic1DA<std::int32_t> ends{};
  for (const auto &bin : f.bins()) {
    chrom_ids.push_back(static_cast<std::int32_t>(bin.chrom().id()));
    starts.push_back(static_cast<std::int32_t>(bin.start()));
    ends.push_back(static_cast<std::int32_t>(bin.end()));
  }

  std::vector<std::string> chrom_names{};
  std::transform(f.chromosomes().begin(), f.chromosomes().end(), std::back_inserter(chrom_names),
                 [&](const hictk::Chromosome &chrom) { return std::string{chrom.name()}; });

  nanobind::dict py_bins_dict{};  // NOLINT

  py_bins_dict["chrom"] =
      pd.attr("Categorical")
          .attr("from_codes")(chrom_ids(), nanobind::arg("categories") = chrom_names,
                              nanobind::arg("validate") = false);
  py_bins_dict["start"] = pd.attr("Series")(starts(), nanobind::arg("copy") = false);
  py_bins_dict["end"] = pd.attr("Series")(ends(), nanobind::arg("copy") = false);

  auto df = pd.attr("DataFrame")(py_bins_dict, nanobind::arg("copy") = false);
  return df;
}

}  // namespace hictkpy
