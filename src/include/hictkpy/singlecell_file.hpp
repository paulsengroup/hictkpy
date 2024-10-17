// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>

#include "hictkpy/nanobind.hpp"

namespace hictkpy::singlecell_file {

[[nodiscard]] bool is_scool_file(const std::filesystem::path& path);

void declare_singlecell_file_class(nanobind::module_& m);

}  // namespace hictkpy::singlecell_file
