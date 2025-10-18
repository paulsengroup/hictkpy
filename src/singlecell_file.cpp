// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/singlecell_file.hpp"

#include <fmt/format.h>

#include <filesystem>
#include <hictk/bin_table.hpp>
#include <hictk/cooler/singlecell_cooler.hpp>
#include <hictk/cooler/validation.hpp>
#include <hictk/file.hpp>
#include <string>
#include <vector>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/file.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/reference.hpp"

namespace nb = nanobind;

namespace hictkpy {

static std::string repr(const SingleCellFile& sclr) {
  return fmt::format(FMT_STRING("SingleCellFile({})"), sclr->path());
}

static SingleCellFile& ctx_enter(SingleCellFile& sclr) { return sclr; }

static void ctx_exit(SingleCellFile& sclr, [[maybe_unused]] nb::handle exc_type,
                     [[maybe_unused]] nb::handle exc_value, [[maybe_unused]] nb::handle traceback) {
  sclr.try_close();
}

static File getitem(const SingleCellFile& sclr, std::string_view cell_id) {
  return File{sclr->open(cell_id)};
}

static nb::dict get_attrs(const SingleCellFile& sclr) {
  nb::dict py_attrs;

  const auto& attrs = sclr->attributes();

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

static auto get_bins(const SingleCellFile& sclr) { return get_bins_from_object(*sclr); }

static nb::list get_cells(const SingleCellFile& sclr) {
  nb::list cells;

  for (const auto& cell : sclr->cells()) {
    cells.append(cell);
  }

  return cells;
}

static auto get_chromosomes(const SingleCellFile& sclr, bool include_ALL) {
  return get_chromosomes_from_object(*sclr, include_ALL);
}

static std::filesystem::path get_path(const SingleCellFile& sclr) {
  return std::filesystem::path{sclr->path()};
}

static std::int64_t get_resolution(const SingleCellFile& sclr) {
  return static_cast<std::int64_t>(sclr->resolution());
}

[[noreturn]] static void throw_closed_file_exc(std::string_view path) {
  throw std::runtime_error(fmt::format(
      FMT_STRING("caught an attempt to access file \"{}\", which has already been closed"), path));
}

SingleCellFile::SingleCellFile(const std::filesystem::path& path)
    : _fp(path), _uri(path.string()) {}

hictk::cooler::SingleCellFile* SingleCellFile::operator->() { return &**this; }

const hictk::cooler::SingleCellFile* SingleCellFile::operator->() const { return &**this; }

hictk::cooler::SingleCellFile& SingleCellFile::operator*() {
  if (_fp.has_value()) {
    return *_fp;
  }
  throw_closed_file_exc(_uri);
}

const hictk::cooler::SingleCellFile& SingleCellFile::operator*() const {
  if (_fp.has_value()) {
    return *_fp;
  }
  throw_closed_file_exc(_uri);
}

void SingleCellFile::close() { _fp.reset(); }

bool SingleCellFile::try_close() noexcept {
  if (!_fp) {
    return true;
  }

  try {
    try {
      close();
      return true;
    } catch (const std::exception& e) {
      nanobind::module_::import_("warnings").attr("warn")(e.what());
    } catch (...) {
      nanobind::module_::import_("warnings")
          .attr("warn")(fmt::format(
              FMT_STRING("an error occurred while closing file \"{}\": unknown error"), _uri));
    }
  } catch (...) {  // NOLINT
  }

  return false;
}

bool SingleCellFile::is_scool(const std::filesystem::path& path) {
  return bool(hictk::cooler::utils::is_scool_file(path.string()));
}

void SingleCellFile::bind(nb::module_& m) {
  auto cooler = m.def_submodule("cooler");

  auto scell_file = nb::class_<SingleCellFile>(
      cooler, "SingleCellFile", "Class representing a file handle to a .scool file.");

  scell_file.def(nb::init<const std::filesystem::path&>(), nb::arg("path"),
                 "Open a single-cell Cooler file (.scool).");

  scell_file.def("__repr__", &repr, nb::rv_policy::move);

  scell_file.def("__enter__", &ctx_enter, nb::rv_policy::reference_internal);

  scell_file.def("__exit__", &ctx_exit,
                 // clang-format off
                 nb::arg("exc_type") = nb::none(),
                 nb::arg("exc_value") = nb::none(),
                 nb::arg("traceback") = nb::none()
                 // clang-format on
  );

  scell_file.def("path", &get_path, "Get the file path.", nb::rv_policy::move);

  scell_file.def("close", &SingleCellFile::close, "Manually close the file handle.");

  scell_file.def("resolution", &get_resolution, "Get the bin size in bp.");
  scell_file.def("chromosomes", &get_chromosomes, nb::arg("include_ALL") = false,
                 "Get the chromosome sizes as a dictionary mapping names to sizes.",
                 nb::rv_policy::take_ownership);
  scell_file.def("bins", &get_bins, nb::sig("def bins(self) -> hictkpy.BinTable"),
                 "Get table of bins.", nb::rv_policy::move);
  scell_file.def("attributes", &get_attrs, "Get file attributes as a dictionary.",
                 nb::rv_policy::take_ownership);
  scell_file.def("cells", &get_cells, nb::sig("def cells(self) -> list[str]"),
                 "Get the list of available cells.", nb::rv_policy::move);
  scell_file.def("__getitem__", &getitem, nb::arg("cell_id"),
                 "Open the Cooler file corresponding to the cell ID given as input.",
                 nb::rv_policy::move);
}

}  // namespace hictkpy
