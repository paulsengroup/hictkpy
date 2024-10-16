// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <hictk/file.hpp>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel_selector.hpp"

namespace hictkpy::file {
void ctor(hictk::File *fp, std::string_view path, std::optional<std::int32_t> resolution,
          std::string_view matrix_type, std::string_view matrix_unit);

[[nodiscard]] std::string repr(const hictk::File &f);

[[nodiscard]] bool is_cooler(std::string_view uri);
[[nodiscard]] bool is_hic(std::string_view uri);

[[nodiscard]] hictkpy::PixelSelector fetch(const hictk::File &f, std::string_view range1,
                                           std::string_view range2, std::string_view normalization,
                                           std::string_view count_type, bool join,
                                           std::string_view query_type);

[[nodiscard]] nanobind::dict attributes(const hictk::File &f);

[[nodiscard]] std::vector<std::string> avail_normalizations(const hictk::File &f);
[[nodiscard]] std::vector<double> weights(const hictk::File &f, std::string_view normalization,
                                          bool divisive = true);

void declare_file_class(nanobind::module_ &m);

}  // namespace hictkpy::file
