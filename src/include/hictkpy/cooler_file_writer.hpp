// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filesystem>
#include <hictk/cooler/singlecell_cooler.hpp>
#include <hictk/reference.hpp>
#include <hictk/tmpdir.hpp>
#include <optional>
#include <string>
#include <string_view>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

class CoolerFileWriter {
  std::string _path{};
  hictk::internal::TmpDir _tmpdir{};
  std::optional<hictk::cooler::SingleCellFile> _w{};
  std::uint32_t _compression_lvl{};
  bool _finalized{false};

 public:
  CoolerFileWriter(std::string_view path_, nanobind::dict chromosomes_, std::uint32_t resolution_,
                   std::string_view assembly = "unknown",
                   const std::filesystem::path& tmpdir = std::filesystem::temp_directory_path(),
                   std::uint32_t compression_lvl = 6);

  [[nodiscard]] std::string path() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;

  [[nodiscard]] const hictk::Reference& chromosomes() const;

  void add_pixels(nanobind::object df);

  void serialize(const std::string& log_lvl_str = "warn");

  [[nodiscard]] std::string repr() const;
  static void bind(nanobind::module_& m);

 private:
  [[nodiscard]] static hictk::cooler::SingleCellFile create_file(
      std::string_view path, hictk::Reference chromosomes, std::uint32_t resolution,
      std::string_view assembly, const std::filesystem::path& tmpdir);
};

}  // namespace hictkpy
