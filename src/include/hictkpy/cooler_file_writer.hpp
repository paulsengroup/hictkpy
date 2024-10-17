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
  CoolerFileWriter() = delete;
  CoolerFileWriter(std::string_view path_, const nanobind::dict& chromosomes_,
                   std::uint32_t resolution_, std::string_view assembly,
                   const std::filesystem::path& tmpdir, std::uint32_t);

  [[nodiscard]] std::string_view path() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;

  [[nodiscard]] const hictk::Reference& chromosomes() const;

  void add_pixels(const nanobind::object& df);

  void finalize(std::string_view log_lvl_str, std::size_t chunk_size, std::size_t update_frequency);

  [[nodiscard]] std::string repr() const;
  static void bind(nanobind::module_& m);

 private:
  [[nodiscard]] static hictk::cooler::SingleCellFile create_file(
      std::string_view path, const hictk::Reference& chromosomes, std::uint32_t resolution,
      std::string_view assembly, const std::filesystem::path& tmpdir);
};

}  // namespace hictkpy
