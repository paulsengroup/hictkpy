// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

#include "hictkpy/singlecell_file.hpp"

namespace py = pybind11;

namespace hictkpy::singlecell_file {

hictk::cooler::SingleCellFile ctor(std::string_view path) {
  return hictk::cooler::SingleCellFile(std::string{path});
}

std::string repr(const hictk::cooler::SingleCellFile& sclr) {
  return fmt::format(FMT_STRING("SingleCellFile({})"), sclr.path());
}

py::dict get_attrs(const hictk::cooler::SingleCellFile& sclr) {
  py::dict py_attrs;

  const auto& attrs = sclr.attributes();

  py_attrs["bin_size"] = attrs.bin_type;
  py_attrs["bin_type"] = attrs.bin_type;
  py_attrs["format"] = attrs.format;
  py_attrs["format_version"] = attrs.format_version;

  for (const auto& key : {"storage-mode", "creation-date", "generated-by", "assembly", "metadata",
                          "format-url", "nbins", "nchroms", "ncells"}) {
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
  if (attrs.ncells.has_value()) {
    py_attrs["ncells"] = *attrs.ncells;
  }

  return py_attrs;
}

hictk::File getitem(const hictk::cooler::SingleCellFile& sclr, std::string_view cell_id) {
  return hictk::File(sclr.open(cell_id));
}

std::vector<std::string> get_cells(const hictk::cooler::SingleCellFile& sclr) {
  std::vector<std::string> cells{sclr.cells().size()};
  std::copy(sclr.cells().begin(), sclr.cells().end(), cells.begin());
  std::sort(cells.begin(), cells.end());

  return cells;
}

}  // namespace hictkpy::singlecell_file
