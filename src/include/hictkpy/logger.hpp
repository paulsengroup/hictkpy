// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog-inl.h>

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

class Logger {
  struct Message {
    std::string payload{};
    double timestamp{};
    spdlog::level::level_enum level{};

    Message() = default;
    Message(const spdlog::details::log_msg &msg);  // NOLINT(*-explicit-conversions)

    [[nodiscard]] nanobind::object to_py_logrecord() const;
    [[nodiscard]] bool is_eoq_signal() const noexcept;
    [[nodiscard]] static Message EOQ();
  };

  std::shared_ptr<spdlog::logger> _logger{};
  std::queue<Message> _msg_queue;
  std::timed_mutex _msg_queue_mtx;
  std::thread _logger_thread{};

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
  void set_level(std::int64_t level) noexcept;
  void set_level(const std::string &level) noexcept;

 private:
  void init_cpp_logger(spdlog::level::level_enum level);
  [[nodiscard]] std::thread spawn_logger_thread();
  std::unique_lock<std::timed_mutex> lock(const std::chrono::milliseconds &timeout);
  void close_msg_queue() noexcept;
};

}  // namespace hictkpy
