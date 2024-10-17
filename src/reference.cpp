// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/reference.hpp"

#include <cstddef>
#include <cstdint>
#include <hictk/reference.hpp>
#include <string>
#include <vector>

#include "hictkpy/nanobind.hpp"

namespace nb = nanobind;

namespace hictkpy {

hictk::Reference chromosome_dict_to_reference(const ChromosomeDict& chromosomes) {
  const auto chrom_list = chromosomes.items();

  std::vector<std::string> chrom_names(chrom_list.size());
  std::vector<std::uint32_t> chrom_sizes(chrom_list.size());

  for (std::size_t i = 0; i < chrom_list.size(); ++i) {
    const auto kv = chrom_list[i];
    chrom_names[i] = nb::cast<std::string>(kv[0]);
    chrom_sizes[i] = nb::cast<std::uint32_t>(kv[1]);
  }

  return {chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()};
}

}  // namespace hictkpy
