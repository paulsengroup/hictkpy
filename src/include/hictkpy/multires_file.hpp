// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <string_view>

#include "hictk/file.hpp"
#include "hictk/multires_file.hpp"

namespace hictkpy::multires_file {
void ctor(hictk::MultiResFile* fp, std::string_view path);

[[nodiscard]] std::string repr(const hictk::MultiResFile& mrf);
}  // namespace hictkpy::multires_file
