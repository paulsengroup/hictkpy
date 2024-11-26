// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog-inl.h>

#include <memory>
#include <string_view>

namespace hictkpy {

class Logger {
  std::shared_ptr<spdlog::logger> _logger{};

 public:
  explicit Logger(spdlog::level::level_enum level = spdlog::level::warn);
  explicit Logger(std::string_view level);

  [[nodiscard]] std::shared_ptr<spdlog::logger> get_logger();

 private:
  [[nodiscard]] static std::shared_ptr<spdlog::logger> init_cpp_logger(
      spdlog::level::level_enum level);
};

}  // namespace hictkpy
