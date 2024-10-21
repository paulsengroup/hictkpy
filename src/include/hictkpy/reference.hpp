// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <hictk/reference.hpp>
#include <string>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

using ChromosomeDict = nanobind::typed<nanobind::dict, std::string, std::uint32_t>;

[[nodiscard]] hictk::Reference chromosome_dict_to_reference(const ChromosomeDict &chromosomes);

template <typename T>
inline nanobind::typed<nanobind::dict, std::string, std::uint32_t> get_chromosomes_from_object(
    const T &obj, bool include_all = false) {
  nanobind::typed<nanobind::dict, std::string, std::uint32_t> py_chroms{};  // NOLINT
  for (const auto &chrom : obj.chromosomes()) {
    if (!include_all && chrom.is_all()) {
      continue;
    }
    const std::string name{chrom.name()};
    py_chroms[name.c_str()] = chrom.size();
  }

  return py_chroms;
}

}  // namespace hictkpy
