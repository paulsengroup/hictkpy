// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>

#include "hictkpy/nanobind.hpp"

namespace hictkpy::file {

[[nodiscard]] bool is_cooler(const std::filesystem::path &uri);
[[nodiscard]] bool is_hic(const std::filesystem::path &uri);

void declare_file_class(nanobind::module_ &m);

}  // namespace hictkpy::file
