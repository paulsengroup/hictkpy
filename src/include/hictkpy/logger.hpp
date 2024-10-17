// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog-inl.h>

#include <memory>
#include <string_view>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

class Logger {
  nanobind::object _py_logger{};
  std::shared_ptr<spdlog::logger> _logger{};

 public:
  explicit Logger(spdlog::level::level_enum level_ = spdlog::level::warn);
  explicit Logger(std::string_view level_);

  [[nodiscard]] std::shared_ptr<spdlog::logger> get_logger();

 private:
  [[nodiscard]] static nanobind::object init_py_logger();
  [[nodiscard]] static std::shared_ptr<spdlog::logger> init_cpp_logger(
      spdlog::level::level_enum level_, nanobind::object py_logger);
};

}  // namespace hictkpy
