// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/logger.hpp"

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <spdlog/sinks/callback_sink.h>
#include <spdlog/spdlog.h>

#include <cassert>
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

template <typename... T>
static void println_noexcept(fmt::format_string<T...> fmt, T&&... args) noexcept {
  try {
    fmt::println(stderr, fmt, std::forward<T>(args)...);
  } catch (...) {  // NOLINT
  }
}

[[nodiscard]] static std::int32_t spdlog_lvl_to_py(spdlog::level::level_enum level) {
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

[[nodiscard]] static double spdlog_to_py_timestamp(const spdlog::details::log_msg& msg) noexcept {
  try {
    const auto t =
        std::chrono::duration_cast<std::chrono::microseconds>(msg.time.time_since_epoch()).count();
    return static_cast<double>(t) / 1.0e6;
  } catch (...) {
    return 0;
  }
}

Logger::Message::Message(const spdlog::details::log_msg& msg)
    : payload(msg.payload.data(), msg.payload.size()),
      timestamp(spdlog_to_py_timestamp(msg)),
      level(msg.level) {}

nb::object Logger::Message::to_py_logrecord() const {
  [[maybe_unused]] const auto gil = GilScopedAcquire<true>::create();
  auto LogRecord = nb::module_::import_("logging").attr("LogRecord");

  auto record = LogRecord("hictkpy",                // name
                          spdlog_lvl_to_py(level),  // log level
                          "",                       // pathname
                          0,                        // line number
                          payload,                  // log message
                          nb::tuple(),              // args
                          nb::none()                // exc_info
  );

  record.attr("created") = timestamp;
  return record;
}

Logger::Logger(spdlog::level::level_enum level) {
  init_cpp_logger(level);
  _logger_thread = spawn_logger_thread();
}

Logger::Logger(std::string_view level) : Logger(spdlog::level::from_str(std::string{level})) {}

Logger::~Logger() noexcept {
  _early_exit = true;
  if (!_logger_thread.joinable()) {
    return;
  }

  try {
    _logger_thread.join();
  } catch (const std::exception& e) {
    println_noexcept(FMT_STRING("Logger: failed to join logger thread: {}"), e.what());
  } catch (...) {  // NOLINT
    println_noexcept(FMT_STRING("Logger: failed to join logger thread: unknown error"));
  }
}

[[nodiscard]] static nb::object get_py_logger() {
  [[maybe_unused]] const auto gil = GilScopedAcquire<true>::create();
  const auto logging = nb::module_::import_("logging");
  return logging.attr("getLogger")("hictkpy");
}

void Logger::init_cpp_logger(spdlog::level::level_enum level_) {
  auto sink = std::make_shared<spdlog::sinks::callback_sink_mt>(
      [this](const spdlog::details::log_msg& msg) {
        try {
          if (HICTKPY_UNLIKELY(!_logger)) {
            return;
          }
          enqueue(msg);
        } catch (const std::exception& e) {  // NOLINT
          println_noexcept(FMT_STRING("Logger: error in log handler: {}"), e.what());
        } catch (...) {  // NOLINT
          println_noexcept(FMT_STRING("Logger: error in log handler: unknown error"));
        }
      });

  sink->set_pattern("%v");
  sink->set_level(level_);

  _logger = std::make_shared<spdlog::logger>("hictkpy", std::move(sink));
  _logger->set_level(level_);
  spdlog::set_default_logger(_logger);
}

std::thread Logger::spawn_logger_thread() {
  return std::thread{[this]() {
    try {
      auto logger = get_py_logger();
      while (!_early_exit) {
        const auto msg = try_dequeue();
        if (msg.has_value()) {
          [[maybe_unused]] const auto gil = GilScopedAcquire<true>::create();
          logger.attr("handle")(msg->to_py_logrecord());
        }
      }
    } catch (const std::exception& e) {
      println_noexcept(FMT_STRING("Logger: an error occurred in logger thread: {}"), e.what());
    } catch (...) {  // NOLINT
      println_noexcept(FMT_STRING("Logger: an error occurred in logger thread: unknown error"));
    }
  }};
}

std::unique_lock<std::timed_mutex> Logger::lock(const std::chrono::milliseconds& timeout) {
  return {_msg_queue_mtx, timeout};
}

void Logger::enqueue(Message msg) noexcept {
  try {
    static const std::chrono::milliseconds timeout{100};
    while (HICTKPY_LIKELY(!_early_exit)) {
      [[maybe_unused]] const auto lck = lock(timeout);
      if (HICTKPY_LIKELY(lck.owns_lock())) {
        _msg_queue.emplace(std::move(msg));
        return;
      }
    }
  } catch (...) {  // NOLINT
  }
}

std::optional<Logger::Message> Logger::try_dequeue() noexcept {
  try {
    static const std::chrono::milliseconds timeout{100};
    [[maybe_unused]] const auto lck = lock(timeout);
    if (HICTKPY_UNLIKELY(lck.owns_lock() && !_msg_queue.empty())) {
      auto msg = std::make_optional(_msg_queue.front());
      _msg_queue.pop();
      return msg;
    }
  } catch (...) {  // NOLINT
  }

  return {};
}

}  // namespace hictkpy
