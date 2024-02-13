// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/string.h>

#include <cstdint>
#include <string>
#include <string_view>

#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/file.hpp"

namespace hictkpy::multires_file {
void ctor(hictk::cooler::MultiResFile* fp, std::string_view path);

[[nodiscard]] std::string repr(const hictk::cooler::MultiResFile& mclr);

[[nodiscard]] nanobind::dict get_attrs(const hictk::cooler::MultiResFile& mclr);

[[nodiscard]] hictk::File getitem(const hictk::cooler::MultiResFile& mclr,
                                  std::uint32_t resolution);
}  // namespace hictkpy::multires_file
