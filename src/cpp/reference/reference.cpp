// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/reference.hpp"

#include <cstddef>
#include <cstdint>
#include <hictk/reference.hpp>
#include <string>
#include <utility>
#include <vector>

#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"

namespace nb = nanobind;

namespace hictkpy {

hictk::Reference chromosome_dict_to_reference(const ChromosomeDict& chromosomes) {
  std::vector<std::string> chrom_names;
  std::vector<std::uint32_t> chrom_sizes;

  [&]() {
    HICTKPY_GIL_SCOPED_ACQUIRE
    const auto chrom_list = chromosomes.items();
    const auto num_chroms = chrom_list.size();

    chrom_names.reserve(num_chroms);
    chrom_sizes.reserve(num_chroms);

    // NOLINTNEXTLINE(*-loop-convert)
    for (std::size_t i = 0; i < chrom_list.size(); ++i) {
      const auto kv = chrom_list[i];
      chrom_names.emplace_back(nb::cast<std::string>(kv[0]));
      chrom_sizes.push_back(nb::cast<std::uint32_t>(kv[1]));
    }

    return std::make_pair(chrom_names, chrom_sizes);
  }();

  return {chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()};
}

}  // namespace hictkpy
