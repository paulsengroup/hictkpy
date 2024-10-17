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
#include "hictkpy/reference.hpp"

namespace hictkpy {

class CoolerFileWriter {
  std::filesystem::path _path{};
  hictk::internal::TmpDir _tmpdir{};
  std::optional<hictk::cooler::SingleCellFile> _w{};
  std::uint32_t _compression_lvl{};
  bool _finalized{false};

 public:
  CoolerFileWriter() = delete;
  CoolerFileWriter(std::filesystem::path path_, const ChromosomeDict& chromosomes_,
                   std::uint32_t resolution_, std::string_view assembly,
                   const std::filesystem::path& tmpdir, std::uint32_t compression_lvl);

  [[nodiscard]] const std::filesystem::path& path() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;

  [[nodiscard]] const hictk::Reference& chromosomes() const;

  void add_pixels(const nanobind::object& df);

  void serialize(const std::string& log_lvl_str = "warn");

  [[nodiscard]] std::string repr() const;
  static void bind(nanobind::module_& m);

 private:
  [[nodiscard]] static hictk::cooler::SingleCellFile create_file(
      std::string_view path, const hictk::Reference& chromosomes, std::uint32_t resolution,
      std::string_view assembly, const std::filesystem::path& tmpdir);
};

}  // namespace hictkpy
