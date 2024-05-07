// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <nanobind/nanobind.h>

#include <cstdint>
#include <string_view>

#include "hictk/bin_table.hpp"

namespace hictkpy::bin_table {

void ctor(hictk::BinTable* bins, nanobind::dict chromosomes, std::uint32_t resolution);

[[nodiscard]] std::string repr(const hictk::BinTable& bins);

[[nodiscard]] nanobind::tuple bin_id_to_coords(const hictk::BinTable& bins, std::uint64_t bin_id);

[[nodiscard]] nanobind::object try_bin_id_to_coords(const hictk::BinTable& bins,
                                                    std::uint64_t bin_id);
[[nodiscard]] nanobind::object try_coords_to_bin_id(const hictk::BinTable& bins,
                                                    std::string_view chrom, std::uint32_t pos);

[[nodiscard]] nanobind::object merge_coords(const hictk::BinTable& bins, nanobind::object df);

[[nodiscard]] nanobind::iterator make_iterable(const hictk::BinTable& bins);

[[nodiscard]] nanobind::object to_df(const hictk::BinTable& bins);

}  // namespace hictkpy::bin_table
