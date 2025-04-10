// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <spdlog/spdlog.h>

#include <cstdint>
#include <hictk/version.hpp>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/cooler_file_writer.hpp"
#include "hictkpy/file.hpp"
#include "hictkpy/hic_file_writer.hpp"
#include "hictkpy/logger.hpp"
#include "hictkpy/multires_file.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel.hpp"
#include "hictkpy/pixel_selector.hpp"
#include "hictkpy/singlecell_file.hpp"

namespace nb = nanobind;
namespace hictkpy {

[[nodiscard]] static hictkpy::Logger init_logger() {
  hictkpy::Logger logger{spdlog::level::debug};
#ifndef _WIN32
  spdlog::set_default_logger(logger.get_logger());
#endif
  return logger;
}

NB_MODULE(_hictkpy, m) {
  // Leaks appear to only occur when the interpreter shuts down abruptly
  nb::set_leak_warnings(false);
  [[maybe_unused]] const auto logger = init_logger();

  m.attr("__hictk_version__") = hictk::config::version::str();

  m.doc() = "Blazing fast toolkit to work with .hic and .cool files.";

  m.def("is_cooler", &file::is_cooler, nb::arg("path"),
        "Test whether path points to a cooler file.");
  m.def("is_mcool_file", &multires_file::is_mcool_file, nb::arg("path"),
        "Test whether path points to a .mcool file.");
  m.def("is_scool_file", &singlecell_file::is_scool_file, nb::arg("path"),
        "Test whether path points to a .scool file.");
  m.def("is_hic", &file::is_hic, nb::arg("path"), "Test whether path points to a .hic file.");

  declare_thin_pixel_class<std::int64_t>(m, "Int");
  declare_thin_pixel_class<double>(m, "FP");
  declare_pixel_class<std::int64_t>(m, "Int");
  declare_pixel_class<double>(m, "FP");

  BinTable::bind(m);

  PixelSelector::bind(m);

  file::declare_file_class(m);

  multires_file::declare_multires_file_class(m);
  singlecell_file::declare_singlecell_file_class(m);

  CoolerFileWriter::bind(m);
  HiCFileWriter::bind(m);
}

}  // namespace hictkpy
