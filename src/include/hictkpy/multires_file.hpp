// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <string_view>

#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/file.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/pixel_selector.hpp"

namespace hictkpy::multires_file {
[[nodiscard]] hictk::cooler::MultiResFile ctor(std::string_view path);

[[nodiscard]] std::string repr(const hictk::cooler::MultiResFile& mclr);

[[nodiscard]] py::dict get_attrs(const hictk::cooler::MultiResFile& mclr);

[[nodiscard]] hictk::File getitem(const hictk::cooler::MultiResFile& mclr,
                                  std::uint32_t resolution);
}  // namespace hictkpy::multires_file
