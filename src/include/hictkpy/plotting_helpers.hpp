// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <hictk/file.hpp>
#include <optional>
#include <string_view>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

[[nodiscard]] nanobind::object coarsen(const hictk::File& f, std::uint32_t target_resolution,
                                       std::optional<std::string_view> range1 = {},
                                       std::optional<std::string_view> range2 = {},
                                       std::optional<std::string_view> normalization = {},
                                       bool floating_point = false, bool exact = false);

}  // namespace hictkpy
