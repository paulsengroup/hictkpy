// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/singlecell_file.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <hictk/bin_table.hpp>
#include <hictk/cooler/singlecell_cooler.hpp>
#include <hictk/cooler/validation.hpp>
#include <hictk/file.hpp>
#include <string>
#include <vector>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/reference.hpp"

namespace nb = nanobind;

namespace hictkpy::singlecell_file {

bool is_scool_file(std::string_view path) {
  return bool(hictk::cooler::utils::is_scool_file(path));
}

static void ctor(hictk::cooler::SingleCellFile* fp, std::string_view path) {
  new (fp) hictk::cooler::SingleCellFile(std::string{path});
}

static std::string repr(const hictk::cooler::SingleCellFile& sclr) {
  return fmt::format(FMT_STRING("SingleCellFile({})"), sclr.path());
}

static nb::dict get_attrs(const hictk::cooler::SingleCellFile& sclr) {
  nb::dict py_attrs;

  const auto& attrs = sclr.attributes();

  py_attrs["bin-size"] = attrs.bin_size;
  py_attrs["bin-type"] = attrs.bin_type == hictk::BinTable::Type::fixed ? "fixed" : "variable";
  py_attrs["format"] = attrs.format;
  py_attrs["format-version"] = attrs.format_version;

  for (const auto& key : {"storage-mode", "creation-date", "generated-by", "assembly", "metadata",
                          "format-url", "nbins", "nchroms", "ncells"}) {
    py_attrs[key] = nb::none();
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

static hictk::File getitem(const hictk::cooler::SingleCellFile& sclr, std::string_view cell_id) {
  return hictk::File(sclr.open(cell_id));
}

static std::vector<std::string> get_cells(const hictk::cooler::SingleCellFile& sclr) {
  std::vector<std::string> cells{sclr.cells().size()};
  std::copy(sclr.cells().begin(), sclr.cells().end(), cells.begin());
  std::sort(cells.begin(), cells.end());

  return cells;
}

void declare_singlecell_file_class(nb::module_& m) {
  auto cooler = m.def_submodule("cooler");

  auto scell_file = nb::class_<hictk::cooler::SingleCellFile>(
      cooler, "SingleCellFile", "Class representing a file handle to a .scool file.");

  scell_file.def("__init__", &singlecell_file::ctor, nb::arg("path"),
                 "Open a single-cell Cooler file (.scool).");

  scell_file.def("__repr__", &singlecell_file::repr);

  scell_file.def("path", &hictk::cooler::SingleCellFile::path, "Get the file path.");
  scell_file.def("resolution", &hictk::cooler::SingleCellFile::resolution,
                 "Get the bin size in bp.");
  scell_file.def("chromosomes", &get_chromosomes_from_file<hictk::cooler::SingleCellFile>,
                 nb::arg("include_all") = false,
                 "Get chromosomes sizes as a dictionary mapping names to sizes.");
  scell_file.def("bins", &get_bins_from_file<hictk::cooler::SingleCellFile>,
                 nb::sig("def bins(self) -> pandas.DataFrame"), "Get bins as a pandas DataFrame.");
  scell_file.def("attributes", &singlecell_file::get_attrs, "Get file attributes as a dictionary.");
  scell_file.def("cells", &singlecell_file::get_cells, "Get the list of available cells.");
  scell_file.def("__getitem__", &singlecell_file::getitem,
                 "Open the Cooler file corresponding to the cell ID given as input.");
}

}  // namespace hictkpy::singlecell_file