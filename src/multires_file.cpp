// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/multires_file.hpp"

#include <fmt/format.h>

#include <hictk/cooler/validation.hpp>
#include <hictk/multires_file.hpp>
#include <string>
#include <string_view>

#include "hictkpy/nanobind.hpp"
#include "hictkpy/reference.hpp"

namespace nb = nanobind;

namespace hictkpy::multires_file {

static void ctor(hictk::MultiResFile* fp, std::string_view path) {
  new (fp) hictk::MultiResFile(std::string{path});
}

static std::string repr(const hictk::MultiResFile& mrf) {
  return fmt::format(FMT_STRING("MultiResFile({})"), mrf.path());
}

bool is_mcool_file(std::string_view path) {
  return bool(hictk::cooler::utils::is_multires_file(path));
}

void declare_multires_file_class(nb::module_& m) {
  auto mres_file = nb::class_<hictk::MultiResFile>(
      m, "MultiResFile", "Class representing a file handle to a .hic or .mcool file");
  mres_file.def("__init__", &multires_file::ctor, nb::arg("path"),
                "Open a multi-resolution Cooler file (.mcool) or .hic file.");

  mres_file.def("__repr__", &multires_file::repr);

  mres_file.def("path", &hictk::MultiResFile::path, "Get the file path.");
  mres_file.def("chromosomes", &get_chromosomes_from_file<hictk::MultiResFile>,
                nb::arg("include_all") = false,
                "Get chromosomes sizes as a dictionary mapping names to sizes.");
  mres_file.def("resolutions", &hictk::MultiResFile::resolutions,
                "Get the list of available resolutions.");
  mres_file.def("__getitem__", &hictk::MultiResFile::open,
                "Open the Cooler or .hic file corresponding to the resolution given as input.");
}

}  // namespace hictkpy::multires_file
