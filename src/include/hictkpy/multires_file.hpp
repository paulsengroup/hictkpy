// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <pybind11/pybind11.h>

#include <cstdint>
#include <string>
#include <string_view>

#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/file.hpp"

namespace hictkpy::multires_file {
[[nodiscard]] hictk::cooler::MultiResFile ctor(std::string_view path);

[[nodiscard]] std::string repr(const hictk::cooler::MultiResFile& mclr);

[[nodiscard]] pybind11::dict get_attrs(const hictk::cooler::MultiResFile& mclr);

[[nodiscard]] hictk::File getitem(const hictk::cooler::MultiResFile& mclr,
                                  std::uint32_t resolution);
}  // namespace hictkpy::multires_file
