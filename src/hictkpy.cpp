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
#include "hictkpy/locking.hpp"
#include "hictkpy/logger.hpp"
#include "hictkpy/multires_file.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel.hpp"
#include "hictkpy/pixel_selector.hpp"
#include "hictkpy/singlecell_file.hpp"

namespace nb = nanobind;
namespace hictkpy {

[[nodiscard]] static std::unique_ptr<Logger> init_logger() {
  try {
    auto logger = std::make_unique<Logger>(spdlog::level::trace);
    [[maybe_unused]] const GilScopedAcquire gil{true};
    nb::module_::import_("atexit").attr("register")(
        nb::cpp_function([logger_ptr = logger.get()]() { logger_ptr->shutdown(); }));
    return logger;
  } catch (const std::exception& e) {
    raise_python_runtime_warning(FMT_STRING("failed to configure hictkpy's logger: {}"), e.what());
  } catch (...) {
    raise_python_runtime_warning(FMT_STRING("failed to configure hictkpy's logger: unknown error"));
  }
  return {};
}

NB_MODULE(_hictkpy, m) {
  // Leaks appear to only occur when the interpreter shuts down abruptly
  nb::set_leak_warnings(false);
  static const auto tsan_proxy_mutex = GilScopedAcquire::try_register_with_tsan();
  static const auto logger = init_logger();

  m.attr("__hictk_version__") = hictk::config::version::str();

  m.doc() = "Blazing fast toolkit to work with .hic and .cool files.";

  m.def("is_cooler", &File::is_cooler, nb::arg("path"),
        "Test whether path points to a cooler file.");
  m.def("is_mcool_file", &MultiResFile::is_mcool, nb::arg("path"),
        "Test whether path points to a .mcool file.");
  m.def("is_scool_file", &SingleCellFile::is_scool, nb::arg("path"),
        "Test whether path points to a .scool file.");
  m.def("is_hic", &File::is_hic, nb::arg("path"), "Test whether path points to a .hic file.");

  auto logging = m.def_submodule("logging");
  logging.def(
      "setLevel",
      [&](const std::variant<std::int64_t, std::string>& level) {
        if (logger) {
          std::visit([&](const auto& level_) { logger->set_level(level_); }, level);
        }
      },
      nb::arg("level"), "Test the log level for hictkpy's logger.");

  BinTable::bind(m);
  Pixel::bind(m);
  PixelSelector::bind(m);

  File::bind(m);
  MultiResFile::bind(m);
  SingleCellFile::bind(m);

  CoolerFileWriter::bind(m);
  HiCFileWriter::bind(m);
}

}  // namespace hictkpy
