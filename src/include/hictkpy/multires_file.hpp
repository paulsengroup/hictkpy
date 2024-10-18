// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>

#include "hictkpy/nanobind.hpp"

namespace hictkpy::multires_file {

[[nodiscard]] bool is_mcool_file(const std::filesystem::path& path);

void declare_multires_file_class(nanobind::module_& m);

}  // namespace hictkpy::multires_file
