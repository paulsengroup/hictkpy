// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <pybind11/pybind11.h>

#include <cstdint>
#include <string>
#include <vector>

#include "hictk/cooler/cooler.hpp"

// This is fine, this header is only supposed to be included in hictkpy.cpp
namespace py = pybind11;

namespace hictkpy::cooler {

inline hictk::cooler::File file_ctor(std::string_view uri) { return hictk::cooler::File(uri); }

inline hictk::cooler::File file_ctor(std::string_view uri, const py::dict& py_chroms,
                                     std::uint32_t bin_size, bool overwrite_if_exists = false) {
  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> chrom_sizes{};

  for (auto it : py_chroms) {
    chrom_names.push_back(py::cast<std::string>(it.first));
    chrom_sizes.push_back(py::cast<std::uint32_t>(it.second));
  }
  const hictk::Reference chroms(chrom_names.begin(), chrom_names.end(), chrom_sizes.begin());
  return hictk::cooler::File::create(uri, chroms, bin_size, overwrite_if_exists);
}

inline bool is_cooler(std::string_view uri) { return bool(hictk::cooler::utils::is_cooler(uri)); }

inline hictk::cooler::File cooler_ctor(std::string_view uri, const py::dict& py_chroms,
                                       std::uint32_t bin_size, bool overwrite_if_exists = false,
                                       bool float_pixels = false) {
  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> chrom_sizes{};

  for (auto it : py_chroms) {
    chrom_names.push_back(py::cast<std::string>(it.first));
    chrom_sizes.push_back(py::cast<std::uint32_t>(it.second));
  }
  const hictk::Reference chroms(chrom_names.begin(), chrom_names.end(), chrom_sizes.begin());
  if (float_pixels) {
    return hictk::cooler::File::create<double>(uri, chroms, bin_size, overwrite_if_exists);
  }
  return hictk::cooler::File::create(uri, chroms, bin_size, overwrite_if_exists);
}

[[nodiscard]] inline py::dict get_cooler_attrs(const hictk::cooler::File& clr) {
  py::dict py_attrs;
  const auto& attrs = clr.attributes();

  py_attrs["bin_size"] = attrs.bin_size;
  py_attrs["bin_type"] = attrs.bin_type;
  py_attrs["format"] = attrs.format;
  py_attrs["format_version"] = attrs.format_version;

  for (const auto& key : {"storage-mode", "creation-date", "generated-by", "assembly", "metadata",
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
    std::visit([&](const auto& sum) { py_attrs["sum"] = sum; }, *attrs.sum);
  }
  if (attrs.cis.has_value()) {
    std::visit([&](const auto& cis) { py_attrs["cis"] = cis; }, *attrs.cis);
  }

  return py_attrs;
}

inline py::object fetch(const hictk::cooler::File& f, std::string_view range1,
                        std::string_view range2, std::string_view normalization,
                        std::string_view count_type, bool join, std::string_view query_type) {
  return file_fetch(f, range1, range2, normalization, count_type, join, query_type);
}

inline py::object fetch_sparse(const hictk::cooler::File& f, std::string_view range1,
                               std::string_view range2, std::string_view normalization,
                               std::string_view count_type, std::string_view query_type) {
  return file_fetch_sparse(f, range1, range2, normalization, count_type, query_type);
}

inline py::object fetch_dense(const hictk::cooler::File& f, std::string_view range1,
                              std::string_view range2, std::string_view normalization,
                              std::string_view count_type, std::string_view query_type) {
  return file_fetch_dense(f, range1, range2, normalization, count_type, query_type);
}
}  // namespace hictkpy::cooler
