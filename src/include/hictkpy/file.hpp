// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <string_view>

#include "hictk/file.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/pixel_selector.hpp"

namespace hictkpy::file {
[[nodiscard]] hictk::File ctor(std::string_view path, std::int32_t resolution,
                               std::string_view matrix_type, std::string_view matrix_unit);

[[nodiscard]] bool is_cooler(std::string_view uri);

[[nodiscard]] hictkpy::PixelSelector fetch(const hictk::File &f, std::string_view range1,
                                           std::string_view range2, std::string_view normalization,
                                           std::string_view count_type, bool join,
                                           std::string_view query_type);

[[nodiscard]] pybind11::dict attributes(const hictk::File &f);

}  // namespace hictkpy::file
