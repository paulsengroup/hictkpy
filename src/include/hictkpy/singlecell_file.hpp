// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>
#include <hictk/cooler/singlecell_cooler.hpp>
#include <optional>
#include <string>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

class SingleCellFile {
  std::optional<hictk::cooler::SingleCellFile> _fp{};
  std::string _uri;

 public:
  SingleCellFile() = delete;
  explicit SingleCellFile(const std::filesystem::path& path);

  SingleCellFile(const SingleCellFile&) = delete;
  SingleCellFile(SingleCellFile&&) noexcept = default;

  ~SingleCellFile() noexcept = default;

  SingleCellFile& operator=(const SingleCellFile&) = delete;
  SingleCellFile& operator=(SingleCellFile&&) noexcept = default;

  [[nodiscard]] const hictk::cooler::SingleCellFile& operator*() const;
  [[nodiscard]] hictk::cooler::SingleCellFile& operator*();

  [[nodiscard]] const hictk::cooler::SingleCellFile* operator->() const;
  [[nodiscard]] hictk::cooler::SingleCellFile* operator->();

  void close();
  bool try_close() noexcept;

  static void bind(nanobind::module_& m);

  [[nodiscard]] static bool is_scool(const std::filesystem::path& path);
};

}  // namespace hictkpy
