// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>
#include <hictk/multires_file.hpp>
#include <optional>
#include <string>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

class MultiResFile {
  std::optional<hictk::MultiResFile> _fp{};
  std::string _uri;

 public:
  MultiResFile() = delete;
  explicit MultiResFile(const std::filesystem::path& path);

  MultiResFile(const MultiResFile&) = delete;
  MultiResFile(MultiResFile&&) noexcept = default;

  ~MultiResFile() noexcept = default;

  MultiResFile& operator=(const MultiResFile&) = delete;
  MultiResFile& operator=(MultiResFile&&) noexcept = default;

  [[nodiscard]] const hictk::MultiResFile& operator*() const;
  [[nodiscard]] hictk::MultiResFile& operator*();

  [[nodiscard]] const hictk::MultiResFile* operator->() const;
  [[nodiscard]] hictk::MultiResFile* operator->();

  void close();
  bool try_close() noexcept;

  static void bind(nanobind::module_& m);

  [[nodiscard]] static bool is_mcool(const std::filesystem::path& path);
};

}  // namespace hictkpy
