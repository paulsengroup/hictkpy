// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/hic/file_writer.hpp>
#include <hictk/reference.hpp>
#include <string>
#include <string_view>
#include <vector>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

class HiCFileWriter {
  hictk::hic::internal::HiCFileWriter _w{};
  bool _finalized{false};

 public:
  HiCFileWriter(std::string_view path, nanobind::dict chromosomes,
                const std::vector<std::uint32_t>& resolutions,
                std::string_view assembly = "unknown", std::size_t n_threads = 1,
                std::size_t chunk_size = 10'000'000,
                const std::filesystem::path& tmpdir = std::filesystem::temp_directory_path(),
                std::uint32_t compression_lvl = 10, bool skip_all_vs_all_matrix = false);
  HiCFileWriter(std::string_view path, nanobind::dict chromosomes, std::uint32_t resolution,
                std::string_view assembly = "unknown", std::size_t n_threads = 1,
                std::size_t chunk_size = 10'000'000,
                const std::filesystem::path& tmpdir = std::filesystem::temp_directory_path(),
                std::uint32_t compression_lvl = 10, bool skip_all_vs_all_matrix = false);

  [[nodiscard]] std::string_view path() const noexcept;
  [[nodiscard]] const std::vector<std::uint32_t>& resolutions() const noexcept;

  [[nodiscard]] const hictk::Reference& chromosomes() const;

  void add_pixels(nanobind::object df);

  void serialize(const std::string& log_lvl_str = "warn");

  [[nodiscard]] std::string repr() const;

  static void bind(nanobind::module_& m);
};

}  // namespace hictkpy
