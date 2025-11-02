// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/logger.hpp"

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <spdlog/sinks/callback_sink.h>
#include <spdlog/spdlog.h>

#include <cstdint>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

#include "hictkpy/common.hpp"
#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"

namespace nb = nanobind;

namespace hictkpy {

[[nodiscard]] static std::int32_t to_py_lvl(spdlog::level::level_enum level) {
  using T = std::underlying_type_t<spdlog::level::level_enum>;
  level = spdlog::level::level_enum{
      std::max(conditional_static_cast<T>(SPDLOG_ACTIVE_LEVEL), conditional_static_cast<T>(level))};

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
  auto sink = std::make_shared<spdlog::sinks::callback_sink_mt>(
      [logger = get_py_logger()](const spdlog::details::log_msg& msg) mutable {
        if constexpr (SPDLOG_ACTIVE_LEVEL < SPDLOG_LEVEL_DEBUG) {
          // useful for debugging and required to avoid deadlocks
          if (msg.level == SPDLOG_LEVEL_TRACE) {
            fmt::println(stderr, FMT_STRING("[{:%Y-%m-%d %T}] [trace]: {}"), msg.time, msg.payload);
            return;
          }
        }

        HICTKPY_GIL_SCOPED_ACQUIRE
        auto msg_py = nb::str(msg.payload.data(), msg.payload.size());
        logger.attr("log")(to_py_lvl(msg.level), msg_py);
      });

  sink->set_pattern("%v");
  sink->set_level(level_);

  auto logger = std::make_shared<spdlog::logger>("hictkpy", std::move(sink));
  logger->set_level(level_);

  return logger;
}

}  // namespace hictkpy
