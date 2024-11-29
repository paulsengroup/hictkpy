// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/logger.hpp"

#include <spdlog/sinks/callback_sink.h>
#include <spdlog/spdlog.h>

#include <cstdint>
#include <memory>
#include <string>
#include <utility>

#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"

namespace nb = nanobind;

namespace hictkpy {

[[nodiscard]] static std::int32_t to_py_lvl(spdlog::level::level_enum level) {
  // https://docs.python.org/3/library/logging.html#logging-levels
  // NOLINTBEGIN(*-avoid-magic-numbers)
  switch (level) {
    case spdlog::level::trace:
      [[fallthrough]];
    case spdlog::level::debug:
      return 10;
    case spdlog::level::info:
      return 20;
    case spdlog::level::warn:
      return 30;
    case spdlog::level::err:
      return 40;
    case spdlog::level::critical:
      [[fallthrough]];
    case spdlog::level::off:
      return 50;
    default:
      unreachable_code();
  }
  // NOLINTEND(*-avoid-magic-numbers)
}

Logger::Logger(spdlog::level::level_enum level) : _logger(init_cpp_logger(level)) {}

Logger::Logger(std::string_view level) : Logger(spdlog::level::from_str(std::string{level})) {}

[[nodiscard]] static nb::object get_py_logger() {
  const auto logging = nb::module_::import_("logging");
  return logging.attr("getLogger")("hictkpy");
}

std::shared_ptr<spdlog::logger> Logger::get_logger() { return _logger; }

std::shared_ptr<spdlog::logger> Logger::init_cpp_logger(
    [[maybe_unused]] spdlog::level::level_enum level_) {
#ifndef _WIN32
  auto sink = std::make_shared<spdlog::sinks::callback_sink_mt>(
      [logger = get_py_logger()](const spdlog::details::log_msg& msg) mutable {
        [[maybe_unused]] const nb::gil_scoped_acquire gil{};
        auto msg_py = nb::str(msg.payload.data(), msg.payload.size());
        logger.attr("log")(to_py_lvl(msg.level), msg_py);
      });

  sink->set_pattern("%v");
  sink->set_level(level_);

  auto logger = std::make_shared<spdlog::logger>("hictkpy", std::move(sink));
  logger->set_level(level_);

  return logger;
#else
  return nullptr;
#endif
}

}  // namespace hictkpy
