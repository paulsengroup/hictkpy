// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cassert>
#include <cstdint>
#include <string>
#include <string_view>
#include <variant>

#include "hictk/balancing/methods.hpp"
#include "hictk/file.hpp"
#include "hictkpy/file.hpp"
#include "hictkpy/pixel_selector.hpp"

namespace py = pybind11;

namespace hictkpy::file {
hictk::File ctor(std::string_view path, std::int32_t resolution, std::string_view matrix_type,
                 std::string_view matrix_unit) {
  return hictk::File{std::string{path}, static_cast<std::uint32_t>(resolution),
                     hictk::hic::ParseMatrixTypeStr(std::string{matrix_type}),
                     hictk::hic::ParseUnitStr(std::string{matrix_unit})};
}

bool is_cooler(std::string_view uri) { return bool(hictk::cooler::utils::is_cooler(uri)); }

hictkpy::PixelSelector fetch(const hictk::File &f, std::string_view range1, std::string_view range2,
                             std::string_view normalization, std::string_view count_type, bool join,
                             std::string_view query_type) {
  if (normalization != "NONE") {
    count_type = "float";
  }

  if (range1.empty()) {
    assert(range2.empty());
    return std::visit(
        [&](const auto &ff) {
          auto sel = ff.fetch(hictk::balancing::Method{normalization});
          using SelT = decltype(sel);
          return hictkpy::PixelSelector(std::make_shared<const SelT>(std::move(sel)), count_type,
                                        join);
        },
        f.get());
  }

  const auto qt =
      query_type == "UCSC" ? hictk::GenomicInterval::Type::UCSC : hictk::GenomicInterval::Type::BED;

  return std::visit(
      [&](const auto &ff) {
        auto sel = range2.empty() || range1 == range2
                       ? ff.fetch(range1, hictk::balancing::Method(normalization), qt)
                       : ff.fetch(range1, range2, hictk::balancing::Method(normalization), qt);
        using SelT = decltype(sel);
        return hictkpy::PixelSelector(std::make_shared<const SelT>(std::move(sel)), count_type,
                                      join);
      },
      f.get());
}

}  // namespace hictkpy::file
