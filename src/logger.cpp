// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/logger.hpp"

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <spdlog/sinks/callback_sink.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <exception>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <thread>
#include <type_traits>
#include <utility>

#include "hictkpy/common.hpp"
#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"

namespace nb = nanobind;

namespace hictkpy {

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

[[nodiscard]] static spdlog::level::level_enum py_to_spdlog_lvl(std::int64_t level) {
  // NOLINTBEGIN(*-avoid-magic-numbers)
  if (level >= 50) {
    return spdlog::level::off;
  }
  if (level >= 40) {
    return spdlog::level::err;
  }
  if (level >= 30) {
    return spdlog::level::warn;
  }
  if (level >= 20) {
    return spdlog::level::info;
  }
  if (level >= 10) {
    return spdlog::level::debug;
  }
  return spdlog::level::trace;
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

[[nodiscard]] static nb::object get_py_logger() {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  const auto logging = nb::module_::import_("logging");
  return logging.attr("getLogger")("hictkpy");
}

static void log_message_py(nb::object& logger, const std::string& level, std::string_view message) {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  std::ignore = logger.attr(level.c_str())(message);
}

static void log_message_py(const std::string& level, std::string_view message) {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  auto logger = get_py_logger();
  log_message_py(logger, level, message);
}

Logger::Message::Message(const spdlog::details::log_msg& msg)
    : payload(msg.payload.data(), msg.payload.size()),
      timestamp(spdlog_to_py_timestamp(msg)),
      level(msg.level) {}

nb::object Logger::Message::to_py_logrecord() const {
  [[maybe_unused]] const GilScopedAcquire gil{true};
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

Logger::Message Logger::Message::EOQ() {
  Message msg{};
  msg.timestamp = -1;
  return msg;
}

bool Logger::Message::is_eoq_signal() const noexcept { return timestamp == -1; }

Logger::Logger() : Logger(spdlog::level::level_enum::warn) {}

Logger::Logger(spdlog::level::level_enum level) {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  log_message_py("debug", "setting up hictkpy's logger...");

  init_cpp_logger(level);
  _logger_thread = spawn_logger_thread();
  log_message_py("debug", "successfully set up the hictkpy's logger!");
}

Logger::Logger(std::string_view level) : Logger(spdlog::level::from_str(std::string{level})) {}

Logger::~Logger() noexcept { shutdown(); }

void Logger::init_cpp_logger(spdlog::level::level_enum level_) {
  auto sink = std::make_shared<spdlog::sinks::callback_sink_mt>(
      [this](const spdlog::details::log_msg& msg) {
        try {
          if (HICTKPY_UNLIKELY(!_logger)) {
            return;
          }
          enqueue(msg);
        } catch (const std::exception& e) {  // NOLINT
          raise_python_user_warning(FMT_STRING("hictkpy::Logger: error in log handler: {}"),
                                    e.what());
        } catch (...) {  // NOLINT
          raise_python_user_warning(
              FMT_STRING("hictkpy::Logger: error in log handler: unknown error"));
        }
      });

  sink->set_pattern("%v");
  sink->set_level(level_);

  _logger = std::make_shared<spdlog::logger>("hictkpy", std::move(sink));
  _logger->set_level(level_);
  spdlog::set_default_logger(_logger);
}

std::thread Logger::spawn_logger_thread() {
  log_message_py("debug", "about to spawn hictkpy's logger thread...");
  return std::thread{[this]() {
    auto logger = []() -> std::optional<nb::object> {
      try {
        return std::make_optional(get_py_logger());
      } catch (const std::exception& e) {
        raise_python_user_warning(
            FMT_STRING("hictkpy::Logger: an error occurred in logger thread: {}"), e.what());
      } catch (...) {  // NOLINT
        raise_python_user_warning(
            FMT_STRING("hictkpy::Logger: an error occurred in logger thread: unknown error"));
      }
      return {};
    }();

    if (!logger.has_value()) {
      return;
    }

    try {
      log_message_py(*logger, "debug", "hictkpy's logger thread successfully started!");
      while (true) {
        const auto msg = try_dequeue();
        if (msg.has_value()) {
          if (msg->is_eoq_signal()) {
            log_message_py("debug",
                           "EOQ signal received: hictkpy's logger thread is shutting down...");

            [[maybe_unused]] const GilScopedAcquire gil{true};
            logger.reset();
            return;
          }
          [[maybe_unused]] const GilScopedAcquire gil{true};
          std::ignore = logger->attr("handle")(msg->to_py_logrecord());
        }
      }
    } catch (const std::exception& e) {
      raise_python_user_warning(
          FMT_STRING("hictkpy::Logger: an error occurred in logger thread: {}"), e.what());
    } catch (...) {  // NOLINT
      raise_python_user_warning(
          FMT_STRING("hictkpy::Logger: an error occurred in logger thread: unknown error"));
    }

    [[maybe_unused]] const GilScopedAcquire gil{true};
    logger.reset();
  }};
}

std::unique_lock<std::timed_mutex> Logger::lock(const std::chrono::milliseconds& timeout) {
  return {_msg_queue_mtx, timeout};
}

void Logger::close_msg_queue() noexcept {
  [[maybe_unused]] const std::scoped_lock lck(_msg_queue_mtx);
  while (true) {
    try {
      _msg_queue.emplace(Message::EOQ());
      return;
    } catch (...) {  // NOLINT
    }
  }
}

void Logger::enqueue(Message msg) noexcept {
  try {
    const std::chrono::milliseconds timeout{100};
    while (true) {
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
    const std::chrono::milliseconds timeout{100};
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

void Logger::set_level(const std::string& py_level) noexcept {
  if (!_logger) {
    return;
  }

  try {
    [[maybe_unused]] const GilScopedAcquire gil{true};
    const auto logging = nb::module_::import_("logging");
    const auto py_level_int = nb::cast<std::int64_t>(logging.attr(py_level.c_str()));

    set_level(py_level_int);
  } catch (const std::exception& e) {
    raise_python_user_warning(FMT_STRING("hictkpy::Logger: failed to change log level: {}"),
                              e.what());
  } catch (...) {
    raise_python_user_warning(
        FMT_STRING("hictkpy::Logger: failed to change log level: unknown error"));
  }
}

void Logger::set_level(std::int64_t py_level) noexcept {
  if (!_logger) {
    return;
  }

  std::optional<nb::object> previous_level{};
  auto reset_log_level = [&]() {
    try {
      if (previous_level.has_value()) {
        get_py_logger().attr("setLevel")(*previous_level);
      }
    } catch (...) {  // NOLINT
    }
  };

  try {
    const auto level = py_to_spdlog_lvl(py_level);
    [[maybe_unused]] const GilScopedAcquire gil{true};
    auto logger = get_py_logger();
    previous_level = logger.attr("level");
    std::ignore = get_py_logger().attr("setLevel")(py_level);

    _logger->set_level(level);
    for (auto& sink : _logger->sinks()) {
      sink->set_level(level);
    }

    return;

  } catch (const std::exception& e) {
    raise_python_user_warning(FMT_STRING("hictkpy::Logger: failed to change log level: {}"),
                              e.what());
  } catch (...) {
    raise_python_user_warning(
        FMT_STRING("hictkpy::Logger: failed to change log level: unknown error"));
  }
  reset_log_level();
}

void Logger::shutdown() noexcept {
  if (!_logger_thread.joinable()) {
    return;
  }

  try {
    log_message_py("debug", "shutting down hictkpy's logger...");
    close_msg_queue();
    {
      [[maybe_unused]] const nb::gil_scoped_release gil{};
      _logger_thread.join();
    }
    log_message_py("debug", "hictkpy's logger was successfully shutdown!");
  } catch (const std::exception& e) {
    raise_python_user_warning(FMT_STRING("hictkpy::Logger: failed to join logger thread: {}"),
                              e.what());
  } catch (...) {  // NOLINT
    raise_python_user_warning(
        FMT_STRING("hictkpy::Logger: failed to join logger thread: unknown error"));
  }
}

}  // namespace hictkpy
