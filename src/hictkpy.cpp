// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <arrow/config.h>
#include <arrow/python/api.h>
#include <spdlog/spdlog.h>

#include <cstdint>
#include <hictk/version.hpp>
#include <stdexcept>

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

static void check_arrow_abi_compat() {
  if (arrow::GetBuildInfo().version_major != ARROW_VERSION_MAJOR ||
      arrow::GetBuildInfo().version_minor != ARROW_VERSION_MINOR) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "Detected Arrow ABI version mismatch!\n"
            "hictkpy was compiled with Arrow v{}, which is not ABI compatible with the currently "
            "installed version of Arrow (v{}).\n"
            "Please install a compatible version of Arrow with e.g. \"pip install "
            "pyarrow=={}.{}\"."),
        ARROW_VERSION_STRING, arrow::GetBuildInfo().version_string, ARROW_VERSION_MAJOR,
        ARROW_VERSION_MINOR));
  }
}

NB_MODULE(_hictkpy, m) {
  [[maybe_unused]] const auto logger = init_logger();

  check_arrow_abi_compat();

  if (arrow::py::import_pyarrow() == -1) {
    throw std::runtime_error("failed to initialize pyarrow runtime");
  }

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

  PixelSelector::bind(m);

  file::declare_file_class(m);

  multires_file::declare_multires_file_class(m);
  singlecell_file::declare_singlecell_file_class(m);

  CoolerFileWriter::bind(m);
  HiCFileWriter::bind(m);
}

}  // namespace hictkpy
