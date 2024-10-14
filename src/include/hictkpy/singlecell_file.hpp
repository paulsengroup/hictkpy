// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// clang-format off
#include "hictkpy/suppress_warnings.hpp"
HICTKPY_DISABLE_WARNING_PUSH
HICTKPY_DISABLE_WARNING_OLD_STYLE_CAST
HICTKPY_DISABLE_WARNING_PEDANTIC
HICTKPY_DISABLE_WARNING_SHADOW
HICTKPY_DISABLE_WARNING_SIGN_CONVERSION
HICTKPY_DISABLE_WARNING_USELESS_CAST
#include <nanobind/nanobind.h>
#include <nanobind/stl/string_view.h>
HICTKPY_DISABLE_WARNING_POP
// clang-format on

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
