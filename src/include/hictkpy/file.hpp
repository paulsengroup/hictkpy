// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// clang-format off
#include "hictkpy/suppress_warnings.hpp"
HICTKPY_DISABLE_WARNING_PUSH
HICTKPY_DISABLE_WARNING_OLD_STYLE_CAST
HICTKPY_DISABLE_WARNING_PEDANTIC
HICTKPY_DISABLE_WARNING_SHADOW
HICTKPY_DISABLE_WARNING_SIGN_CONVERSION
HICTKPY_DISABLE_WARNING_USELESS_CAST

#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
HICTKPY_DISABLE_WARNING_POP
// clang-format on

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

#include "hictk/file.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/pixel_selector.hpp"

namespace hictkpy::file {
void ctor(hictk::File *fp, std::string_view path, std::optional<std::int32_t> resolution,
          std::string_view matrix_type, std::string_view matrix_unit);

[[nodiscard]] std::string repr(const hictk::File &f);

[[nodiscard]] bool is_cooler(std::string_view uri);
[[nodiscard]] bool is_mcool_file(std::string_view uri);
[[nodiscard]] bool is_scool_file(std::string_view uri);
[[nodiscard]] bool is_hic(std::string_view uri);

[[nodiscard]] hictkpy::PixelSelector fetch(const hictk::File &f, std::string_view range1,
                                           std::string_view range2, std::string_view normalization,
                                           std::string_view count_type, bool join,
                                           std::string_view query_type);

[[nodiscard]] nanobind::dict attributes(const hictk::File &f);

[[nodiscard]] std::vector<std::string> avail_normalizations(const hictk::File &f);
[[nodiscard]] std::vector<double> weights(const hictk::File &f, std::string_view normalization,
                                          bool divisive = true);

}  // namespace hictkpy::file
