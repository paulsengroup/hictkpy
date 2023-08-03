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

[[nodiscard]] inline py::dict get_cooler_attrs(const hictk::cooler::File &clr) {
  py::dict py_attrs;
  const auto &attrs = clr.attributes();

  py_attrs["bin_size"] = attrs.bin_size;
  py_attrs["bin_type"] = attrs.bin_type;
  py_attrs["format"] = attrs.format;
  py_attrs["format_version"] = attrs.format_version;

  for (const auto &key : {"storage-mode", "creation-date", "generated-by", "assembly", "metadata",
                          "format-url", "nbins", "nchroms", "nnz", "sum", "cis"}) {
    py_attrs[key] = pybind11::none();
  }

  if (attrs.storage_mode.has_value()) {
    py_attrs["storage-mode"] = *attrs.storage_mode;
  }

  if (attrs.creation_date.has_value()) {
    py_attrs["creation-date"] = *attrs.creation_date;
  }
  if (attrs.generated_by.has_value()) {
    py_attrs["generated-by"] = *attrs.generated_by;
  }
  if (attrs.assembly.has_value()) {
    py_attrs["assembly"] = *attrs.assembly;
  }
  if (attrs.metadata.has_value()) {
    py_attrs["metadata"] = *attrs.metadata;
  }
  if (attrs.format_url.has_value()) {
    py_attrs["format-url"] = *attrs.format_url;
  }
  if (attrs.nbins.has_value()) {
    py_attrs["nbins"] = *attrs.nbins;
  }
  if (attrs.nchroms.has_value()) {
    py_attrs["nchroms"] = *attrs.nchroms;
  }
  if (attrs.nnz.has_value()) {
    py_attrs["nnz"] = *attrs.nnz;
  }
  if (attrs.sum.has_value()) {
    std::visit([&](const auto &sum) { py_attrs["sum"] = sum; }, *attrs.sum);
  }
  if (attrs.cis.has_value()) {
    std::visit([&](const auto &cis) { py_attrs["cis"] = cis; }, *attrs.cis);
  }

  return py_attrs;
}

[[nodiscard]] inline py::dict get_hic_attrs(const hictk::hic::File &hf) {
  py::dict py_attrs;

  py_attrs["bin_size"] = hf.bin_size();
  py_attrs["format"] = "HIC";
  py_attrs["format_version"] = hf.version();
  py_attrs["assembly"] = hf.assembly();
  py_attrs["format-url"] = "https://github.com/aidenlab/hic-format";
  py_attrs["nbins"] = hf.bins().size();
  py_attrs["nchroms"] = hf.nchroms();

  return py_attrs;
}

pybind11::dict attributes(const hictk::File &f) {
  if (f.is_cooler()) {
    return get_cooler_attrs(f.get<hictk::cooler::File>());
  }
  return get_hic_attrs(f.get<hictk::hic::File>());
}

}  // namespace hictkpy::file
