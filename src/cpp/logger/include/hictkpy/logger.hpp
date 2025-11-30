// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog-inl.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <queue>
#include <string>
#include <string_view>
#include <thread>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

struct LogMessage {
  std::string payload{};
  double timestamp{};
  spdlog::level::level_enum level{};

  LogMessage() = default;
  LogMessage(const spdlog::details::log_msg &msg);  // NOLINT(*-explicit-conversions)
  LogMessage(std::string payload_, double timestamp_, spdlog::level::level_enum level_);

  [[nodiscard]] nanobind::object to_py_logrecord() const;
  [[nodiscard]] bool is_eoq_signal() const noexcept;
  [[nodiscard]] static LogMessage EOQ();
  void serialize(std::string &buff) const;
  static void deserialize(std::string_view buff, LogMessage &msg);
  [[nodiscard]] static LogMessage deserialize(std::string_view buff);
};

class MessageQueue {
  using Message = LogMessage;
  std::string _buff;
  std::queue<Message> _msg_queue;
  std::timed_mutex _mtx;
  std::atomic<bool> _closed{false};

 public:
  MessageQueue() = default;

  MessageQueue(const MessageQueue &) = delete;
  MessageQueue(MessageQueue &&) noexcept = delete;

  ~MessageQueue() noexcept = default;

  MessageQueue &operator=(const MessageQueue &) = delete;
  MessageQueue &operator=(MessageQueue &&) noexcept = delete;

  explicit operator bool() const noexcept;

  void enqueue(Message msg,  // NOLINTNEXTLINE(*-avoid-magic-numbers)
               std::chrono::milliseconds wait_time = std::chrono::milliseconds{50}) noexcept;
  [[nodiscard]] std::optional<Message> try_dequeue() noexcept;
  [[nodiscard]] std::optional<Message> try_dequeue_timed() noexcept;
  template <typename Duration>
  [[nodiscard]] std::optional<Message> try_dequeue_timed(Duration duration) noexcept;
  void dequeue_all(std::vector<Message> &msgs) noexcept;
  void send_eoq() noexcept;
};

class Logger {
  using Message = LogMessage;

  std::shared_ptr<spdlog::logger> _logger{};
  MessageQueue _msg_queue{};
  std::thread _logger_thread{};
  std::atomic<bool> _return_immediately{false};

 public:
  Logger();
  explicit Logger(spdlog::level::level_enum level);
  explicit Logger(std::string_view level);

  Logger(const Logger &) = delete;
  Logger(Logger &&) noexcept = delete;

  ~Logger() noexcept;

  Logger &operator=(const Logger &) = delete;
  Logger &operator=(Logger &&) noexcept = delete;

  void enqueue(Message msg) noexcept;
  [[nodiscard]] std::optional<Message> try_dequeue() noexcept;
  void shutdown() noexcept;
  void reset_after_fork() noexcept;
  void set_level(std::int64_t level) noexcept;
  void set_level(const std::string &level) noexcept;
  void flush() noexcept;

 private:
  [[nodiscard]] std::shared_ptr<spdlog::logger> init_cpp_logger(
      spdlog::level::level_enum level) noexcept;
  [[nodiscard]] std::thread start_logger_thread() noexcept;
};

}  // namespace hictkpy
