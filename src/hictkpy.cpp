// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cstdint>
#include <exception>
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

static void set_nanobind_leak_warnings() {
#ifndef NDEBUG
  const auto x = true;
#else
  const auto x = false;
#endif
  // Leaks appear to only occur when the interpreter shuts down abruptly
  nb::set_leak_warnings(x);
}

static void handle_proc_forking(Logger* logger_ptr, bool print_warning = true) {
  if (logger_ptr) {
    if (print_warning) {
      raise_python_user_warning(FMT_STRING(
          "hictkpy: detected a call to fork():\n"
          "hictkpy's logger does not support multiprocessing when using fork() as start method.\n"
          "Please change process start method to spawn or forkserver.\n"
          "For more details, refer to Python's documentation:\n"
          "https://docs.python.org/3/library/"
          "multiprocessing.html#multiprocessing.set_start_method"));
    }
    logger_ptr->reset_after_fork();
  }
}

[[nodiscard]] static std::unique_ptr<Logger> init_logger() {
  try {
    auto logger = std::make_unique<Logger>(spdlog::level::trace);
    [[maybe_unused]] const GilScopedAcquire gil{true};
    nb::module_::import_("atexit").attr("register")(
        nb::cpp_function([logger_ptr = logger.get()]() { logger_ptr->shutdown(); }));

    auto os = nb::module_::import_("os");
    if (nb::hasattr(os, "register_at_fork")) {
      os.attr("register_at_fork")(
          nb::arg("after_in_parent") = nb::cpp_function(
              [logger_ptr = logger.get()]() { handle_proc_forking(logger_ptr, true); }),
          nb::arg("after_in_child") = nb::cpp_function(
              [logger_ptr = logger.get()]() { handle_proc_forking(logger_ptr, false); }));
    }

    return logger;
  } catch (const std::exception& e) {
    raise_python_runtime_warning(FMT_STRING("failed to configure hictkpy's logger: {}"), e.what());
  } catch (...) {
    raise_python_runtime_warning(FMT_STRING("failed to configure hictkpy's logger: unknown error"));
  }
  return {};
}

NB_MODULE(_hictkpy, m) {
  set_nanobind_leak_warnings();
  [[maybe_unused]] static const auto& tsan_proxy_mutex = GilScopedAcquire::try_register_with_tsan();
  static const auto logger = init_logger();

  cooler::init_global_state();

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
      nb::arg("level"),
      "Test the log level for hictkpy's logger.\n"
      "Accepts the predefined levels defined by the logging module.");

  logging.def(
      "flush",
      []() {
        if (logger) {
          logger->flush();
        }
      },
      "Flush all log messages.");

  logging.def(
      "_log",
      [&](std::string_view level, std::string_view msg) {
        if (level == "trace") {
          SPDLOG_TRACE(FMT_STRING("{}"), msg);
        } else if (level == "debug") {
          SPDLOG_DEBUG(FMT_STRING("{}"), msg);
        } else if (level == "info") {
          SPDLOG_INFO(FMT_STRING("{}"), msg);
        } else if (level == "warn") {
          SPDLOG_WARN(FMT_STRING("{}"), msg);
        } else if (level == "err") {
          SPDLOG_ERROR(FMT_STRING("{}"), msg);
        } else if (level == "critical") {
          SPDLOG_CRITICAL(FMT_STRING("{}"), msg);
        } else {
          SPDLOG_INFO(FMT_STRING("{}"), msg);
        }
      },
      nb::arg("level"), nb::arg("msg"));

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
