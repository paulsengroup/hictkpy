// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <nanobind/nanobind.h>
#include <spdlog/spdlog.h>

#include <hictk/cooler/cooler.hpp>
#include <hictk/cooler/singlecell_cooler.hpp>
#include <hictk/hic/file_writer.hpp>
#include <hictk/reference.hpp>
#include <optional>
#include <variant>
#include <vector>

namespace hictkpy {

class HiCFileWriter {
  hictk::hic::internal::HiCFileWriter _w{};
  bool _finalized{false};

 public:
  HiCFileWriter(std::string_view path, hictk::Reference chromosomes,
                const std::vector<std::uint32_t>& resolutions,
                std::string_view assembly = "unknown", std::size_t n_threads = 1,
                std::size_t chunk_size = 10'000'000,
                const std::filesystem::path& tmpdir = std::filesystem::temp_directory_path(),
                std::uint32_t compression_lvl = 10, bool skip_all_vs_all_matrix = false);

  [[nodiscard]] std::string_view path() const noexcept;
  [[nodiscard]] const std::vector<std::uint32_t>& resolutions() const noexcept;

  [[nodiscard]] const hictk::Reference& chromosomes() const;

  void add_pixels(nanobind::object df);

  void serialize(const std::string& log_lvl_str = "warn");
};

void hic_file_writer_ctor(hictkpy::HiCFileWriter* fp, std::string_view path,
                          nanobind::dict chromosomes, const std::vector<std::uint32_t>& resolutions,
                          std::string_view assembly, std::size_t n_threads, std::size_t chunk_size,
                          std::string_view tmpdir, std::uint32_t compression_lvl,
                          bool skip_all_vs_all_matrix);

void hic_file_writer_ctor_single_res(hictkpy::HiCFileWriter* fp, std::string_view path,
                                     nanobind::dict chromosomes, std::uint32_t resolution,
                                     std::string_view assembly, std::size_t n_threads,
                                     std::size_t chunk_size, std::string_view tmpdir,
                                     std::uint32_t compression_lvl, bool skip_all_vs_all_matrix);

[[nodiscard]] std::string hic_file_writer_repr(hictkpy::HiCFileWriter& w);

class CoolFileWriter {
  std::string _path{};
  hictk::internal::TmpDir _tmpdir{};
  std::optional<hictk::cooler::SingleCellFile> _w{};
  std::uint32_t _compression_lvl{};
  bool _finalized{false};

 public:
  CoolFileWriter(std::string_view path_, hictk::Reference chromosomes_, std::uint32_t resolution_,
                 std::string_view assembly = "unknown",
                 const std::filesystem::path& tmpdir = std::filesystem::temp_directory_path(),
                 std::uint32_t compression_lvl = 6);

  [[nodiscard]] std::string path() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;

  [[nodiscard]] const hictk::Reference& chromosomes() const;

  void add_pixels(nanobind::object df);

  void serialize(const std::string& log_lvl_str = "warn");

 private:
  [[nodiscard]] static hictk::cooler::SingleCellFile create_file(
      std::string_view path, hictk::Reference chromosomes, std::uint32_t resolution,
      const std::filesystem::path& tmpdir);
};

void cool_file_writer_ctor(hictkpy::CoolFileWriter* fp, std::string_view path,
                           nanobind::dict chromosomes, std::uint32_t resolution,
                           std::string_view assembly, std::string_view tmpdir,
                           std::uint32_t compression_lvl);

[[nodiscard]] std::string cool_file_writer_repr(hictkpy::CoolFileWriter& w);

}  // namespace hictkpy
