// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string_view>

#include "hictkpy/nanobind.hpp"

namespace hictkpy::multires_file {

[[nodiscard]] bool is_mcool_file(std::string_view uri);

void declare_multires_file_class(nanobind::module_& m);

}  // namespace hictkpy::multires_file
