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

static constexpr bool PRINTDEBUG_LOGGING{false};

template <typename... T>
static void printdebug(fmt::format_string<T...> fmt, T&&... args) noexcept {
  if constexpr (PRINTDEBUG_LOGGING) {
    println_stderr_noexcept(fmt, std::forward<T>(args)...);
  }
}

static void printdebug(std::string_view msg) noexcept { printdebug(FMT_STRING("{}"), msg); }

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
      return 50;
    case spdlog::level::off:
      return 99;
    default:
      unreachable_code();
  }
  // NOLINTEND(*-avoid-magic-numbers)
}

[[nodiscard]] static spdlog::level::level_enum py_to_spdlog_lvl(std::int64_t level) {
  // NOLINTBEGIN(*-avoid-magic-numbers)
  if (level > 50) {
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
  // NOLINTEND(*-avoid-magic-numbers)
}

[[nodiscard]] static double spdlog_to_py_timestamp(const spdlog::details::log_msg& msg) noexcept {
  try {
    const auto t =
        std::chrono::duration_cast<std::chrono::microseconds>(msg.time.time_since_epoch()).count();
    return static_cast<double>(t) / 1.0e6;  // NOLINT(*-avoid-magic-numbers)
  } catch (...) {
    return 0;
  }
}

[[nodiscard]] static nb::object get_py_logger() {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  const auto logging = nb::module_::import_("logging");
  return logging.attr("getLogger")("hictkpy");
}

LogMessage::LogMessage(const spdlog::details::log_msg& msg)
    : payload(msg.payload.data(), msg.payload.size()),
      timestamp(spdlog_to_py_timestamp(msg)),
      level(msg.level) {}

LogMessage::LogMessage(std::string payload_, double timestamp_, spdlog::level::level_enum level_)
    : payload(std::move(payload_)), timestamp(timestamp_), level(level_) {}

nb::object LogMessage::to_py_logrecord() const {
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

LogMessage LogMessage::EOQ() { return {"", -1, spdlog::level::off}; }

void LogMessage::serialize(std::string& buff) const {
  // NOLINTBEGIN(*-pro-type-reinterpret-cast)
  const std::string_view timestamp_view{reinterpret_cast<const char*>(&timestamp),
                                        sizeof(timestamp)};
  const std::string_view level_view{reinterpret_cast<const char*>(&level), sizeof(level)};
  // NOLINTEND(*-pro-type-reinterpret-cast)

  buff.clear();
  buff.reserve(timestamp_view.size() + level_view.size() + payload.size());

  buff.append(timestamp_view);
  buff.append(level_view);
  buff.append(payload);
}

void LogMessage::deserialize(std::string_view buff, LogMessage& msg) {
  if (buff.size() < sizeof(timestamp) + sizeof(level)) {
    // empty or corrupted message
    msg.payload.clear();
    msg.timestamp = 0;
    msg.level = spdlog::level::off;
    return;
  }

  std::memcpy(&msg.timestamp, buff.data(), sizeof(timestamp));
  buff.remove_prefix(sizeof(timestamp));

  std::memcpy(&msg.level, buff.data(), sizeof(level));
  buff.remove_prefix(sizeof(level));

  msg.payload.assign(buff);
}

LogMessage LogMessage::deserialize(std::string_view buff) {
  LogMessage msg;
  deserialize(buff, msg);
  return msg;
}

bool LogMessage::is_eoq_signal() const noexcept { return timestamp == -1; }

MessageQueue::operator bool() const noexcept { return !_closed; }

void MessageQueue::enqueue(Message msg, std::chrono::milliseconds wait_time) noexcept {
  try {
    while (!_closed) {
      [[maybe_unused]] const std::unique_lock lck{_mtx, wait_time};
      if (lck.owns_lock()) {
        _msg_queue.emplace(std::move(msg));
        return;
      }
    }
  } catch (...) {  // NOLINT
  }
}

std::optional<LogMessage> MessageQueue::try_dequeue() noexcept {
  try {
    while (!_closed) {
      [[maybe_unused]] const std::unique_lock lck{_mtx};
      if (lck.owns_lock() && !_msg_queue.empty()) {
        auto msg = _msg_queue.front();
        _msg_queue.pop();
        if (msg.is_eoq_signal()) {
          _closed = true;
        }
        return msg;
      }
    }
  } catch (...) {  // NOLINT
  }
  return {};
}

[[nodiscard]] std::optional<LogMessage> MessageQueue::try_dequeue_timed() noexcept {
  return try_dequeue_timed(std::chrono::milliseconds{100});  // NOLINT(*-avoid-magic-numbers)
}

template <typename Duration>
[[nodiscard]] std::optional<LogMessage> MessageQueue::try_dequeue_timed(
    Duration duration) noexcept {
  try {
    while (!_closed) {
      [[maybe_unused]] const std::unique_lock lck{_mtx, duration};
      if (lck.owns_lock() && !_msg_queue.empty()) {
        auto msg = _msg_queue.front();
        _msg_queue.pop();
        if (msg.is_eoq_signal()) {
          _closed = true;
        }
        return msg;
      }
    }
  } catch (...) {  // NOLINT
  }
  return {};
}

void MessageQueue::dequeue_all(std::vector<Message>& msgs) noexcept {
  try {
    msgs.clear();
    [[maybe_unused]] const std::unique_lock lck{_mtx};
    msgs.reserve(_msg_queue.size());
    while (!_msg_queue.empty()) {
      msgs.emplace_back(_msg_queue.front());
      _msg_queue.pop();
    }
  } catch (...) {  // NOLINT
  }
}

void MessageQueue::send_eoq() noexcept { enqueue(Message::EOQ()); }

std::shared_ptr<spdlog::logger> Logger::init_cpp_logger(spdlog::level::level_enum level_) noexcept {
  printdebug("hictkpy::Logger::Logger(): setting up hictkpy's logger...");

  try {
    auto logger = std::make_shared<spdlog::logger>("hictkpy");
    if (!logger) {
      raise_python_user_warning("hictkpy::Logger: setup failed: logging is disabled");
      return nullptr;
    }

    logger->set_level(level_);

    auto sink = std::make_shared<spdlog::sinks::callback_sink_mt>(
        [this, logger](const spdlog::details::log_msg& msg) {
          try {
            if (HICTKPY_UNLIKELY(!logger)) {
              return;
            }
            enqueue(msg);
          } catch (const std::exception& e) {
            raise_python_user_warning(FMT_STRING("hictkpy::Logger: error in log handler: {}"),
                                      e.what());
          } catch (...) {  // NOLINT
            raise_python_user_warning(
                FMT_STRING("hictkpy::Logger: error in log handler: unknown error"));
          }
        });

    sink->set_pattern("%v");
    sink->set_level(level_);
    logger->sinks().emplace_back(std::move(sink));

    printdebug("hictkpy::Logger::Logger(): successfully set up the hictkpy's logger!");
    return logger;

  } catch (const std::exception& e) {
    raise_python_user_warning(FMT_STRING("hictkpy::Logger: setup failed: {}\nlogging is disabled"),
                              e.what());
  } catch (...) {
    raise_python_user_warning("hictkpy::Logger: setup failed: logging is disabled");
  }
  return nullptr;
}

[[nodiscard]] static std::optional<nb::object> get_logger_noexcept() noexcept {
  try {
    return get_py_logger();
  } catch (const std::exception& e) {
    raise_python_user_warning(FMT_STRING("hictkpy::Logger: an error occurred in logger thread: {}"),
                              e.what());
  } catch (...) {  // NOLINT
    raise_python_user_warning(
        FMT_STRING("hictkpy::Logger: an error occurred in logger thread: unknown error"));
  }
  return {};
}

// NOLINTNEXTLINE(*-function-cognitive-complexity)
std::thread Logger::start_logger_thread() noexcept {
  auto thr_started = std::make_shared<std::atomic<bool>>(false);
  std::unique_ptr<std::thread> thr{};
  try {
    thr = std::make_unique<std::thread>([this, thr_started]() mutable {
      assert(thr_started);
      try {
        *thr_started = true;
        while (true) {
          const auto msg = _msg_queue.try_dequeue_timed();
          if (_return_immediately) {
            printdebug("hictkpy::Logger: logger thread is returning immediately");
            return;
          }
          if (msg.has_value()) {
            if (msg->is_eoq_signal()) {
              printdebug("hictkpy::Logger: EOQ signal received: logger thread has been shutdown");
              return;
            }
            [[maybe_unused]] const GilScopedAcquire gil{true};
            if (auto logger = get_logger_noexcept(); logger.has_value()) {
              std::ignore = logger->attr("handle")(msg->to_py_logrecord());
            }
          }
        }
      } catch (const nb::python_error& e) {
        printdebug(
            FMT_STRING(
                "hictkpy::Logger: logger thread returning due to the following exception: {}"),
            e.what());
      } catch (const std::exception& e) {
        raise_python_user_warning(
            FMT_STRING("hictkpy::Logger: an error occurred in logger thread: {}"), e.what());
      } catch (...) {  // NOLINT
        raise_python_user_warning(
            FMT_STRING("hictkpy::Logger: an error occurred in logger thread: unknown error"));
      }
    });
  } catch (...) {
    raise_python_user_warning(FMT_STRING("hictkpy::Logger: failed to start logger thread!"));
    if (!!thr && thr->joinable()) {
      thr->join();
    }
  }

  while (!*thr_started) {  // NOLINTNEXTLINE(*-avoid-magic-numbers)
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
  }
  printdebug("hictkpy::Logger: logger thread successfully started!");

  return std::move(*thr);
}

Logger::Logger() : Logger(spdlog::level::level_enum::warn) {}

Logger::Logger(spdlog::level::level_enum level)
    : _logger(init_cpp_logger(level)), _logger_thread(start_logger_thread()) {
  if (_logger_thread.joinable()) {
    spdlog::set_default_logger(_logger);
  } else if (auto logger = spdlog::default_logger(); logger) {
    logger->set_level(spdlog::level::off);
  }
}

Logger::Logger(std::string_view level) : Logger(spdlog::level::from_str(std::string{level})) {}

Logger::~Logger() noexcept {
  _return_immediately = true;
  if (_logger_thread.joinable()) {
    _logger_thread.join();
  }
}

void Logger::enqueue(Message msg) noexcept { _msg_queue.enqueue(std::move(msg)); }

std::optional<Logger::Message> Logger::try_dequeue() noexcept { return _msg_queue.try_dequeue(); }

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
  // NOLINTNEXTLINE(*-deadcode.DeadStores)
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

void Logger::flush() noexcept {
  try {
    [[maybe_unused]] const GilScopedAcquire gil{true};
    std::vector<Message> msgs;
    _msg_queue.dequeue_all(msgs);
    auto logger = get_logger_noexcept();
    if (!logger.has_value()) {
      return;
    }

    for (const auto& msg : msgs) {
      if (msg.is_eoq_signal()) {
        _return_immediately = true;
        return;
      }
      std::ignore = logger->attr("handle")(msg.to_py_logrecord());
    }
  } catch (...) {  // NOLINT
  }
}

void Logger::shutdown() noexcept {
  try {
    [[maybe_unused]] const nb::gil_scoped_release gil{};
    printdebug("hictkpy::Logger::shutdown() called!");
    if (_logger_thread.joinable()) {
      _msg_queue.send_eoq();
      printdebug("hictkpy::Logger::shutdown(): sending EOQ signal!");
      _logger_thread.join();
    }
    if (auto logger = spdlog::default_logger(); !!logger) {
      printdebug("hictkpy::Logger::shutdown(): disabling spdlog logger...");
      logger->set_level(spdlog::level::off);
    }
  } catch ([[maybe_unused]] const std::exception& e) {
    printdebug(FMT_STRING("hictkpy::Logger::shutdown(): failed to join logger thread: {}"),
               e.what());
  } catch (...) {  // NOLINT
    printdebug("hictkpy::Logger::shutdown(): failed to join logger thread: unknown error");
  }
}

void Logger::reset_after_fork() noexcept {
  printdebug("hictkpy::Logger::reset_after_fork() called!");
  [[maybe_unused]] const GilScopedAcquire gil{true};
  shutdown();
}

}  // namespace hictkpy
