// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/string_view.h>

#include <cstdint>
#include <string>
#include <string_view>

#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/file.hpp"

namespace hictkpy::singlecell_file {
void ctor(hictk::cooler::SingleCellFile* fp, std::string_view path);

[[nodiscard]] std::string repr(const hictk::cooler::SingleCellFile& sclr);

[[nodiscard]] nanobind::dict get_attrs(const hictk::cooler::SingleCellFile& sclr);

[[nodiscard]] hictk::File getitem(const hictk::cooler::SingleCellFile& sclr,
                                  std::string_view cell_id);

[[nodiscard]] std::vector<std::string> get_cells(const hictk::cooler::SingleCellFile& sclr);
}  // namespace hictkpy::singlecell_file
