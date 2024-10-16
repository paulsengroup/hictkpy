// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/hic/file_writer.hpp>
#include <hictk/reference.hpp>
#include <hictk/tmpdir.hpp>
#include <string>
#include <string_view>
#include <vector>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

class HiCFileWriter {
  hictk::internal::TmpDir _tmpdir{};
  hictk::hic::internal::HiCFileWriter _w{};
  bool _finalized{false};

 public:
  HiCFileWriter(std::string_view path, const nanobind::dict& chromosomes,
                const std::vector<std::uint32_t>& resolutions, std::string_view assembly,
                std::size_t n_threads, std::size_t chunk_size, const std::filesystem::path& tmpdir,
                std::uint32_t compression_lvl, bool skip_all_vs_all_matrix);
  HiCFileWriter(std::string_view path, const nanobind::dict& chromosomes, std::uint32_t resolution,
                std::string_view assembly, std::size_t n_threads, std::size_t chunk_size,
                const std::filesystem::path& tmpdir, std::uint32_t compression_lvl,
                bool skip_all_vs_all_matrix);

  [[nodiscard]] std::string_view path() const noexcept;
  [[nodiscard]] const std::vector<std::uint32_t>& resolutions() const noexcept;

  [[nodiscard]] const hictk::Reference& chromosomes() const;

  void add_pixels(const nanobind::object& df);

  void serialize(const std::string& log_lvl_str = "warn");

  [[nodiscard]] std::string repr() const;

  static void bind(nanobind::module_& m);
};

}  // namespace hictkpy
