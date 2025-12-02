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
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/file.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/reference.hpp"
#include "hictkpy/table.hpp"

namespace hictkpy {

class HiCFileWriter {
  std::string _path;
  std::uint32_t _base_resolution{};
  std::optional<hictk::internal::TmpDir> _tmpdir{};
  std::optional<hictk::hic::internal::HiCFileWriter> _w{};

 public:
  HiCFileWriter() = delete;
  HiCFileWriter(const std::filesystem::path& path_, const ChromosomeDict& chromosomes,
                const std::vector<std::uint32_t>& resolutions_, std::string_view assembly,
                std::size_t n_threads, std::size_t chunk_size, const std::filesystem::path& tmpdir_,
                std::uint32_t compression_lvl, bool skip_all_vs_all_matrix);
  HiCFileWriter(const std::filesystem::path& path_, const ChromosomeDict& chromosomes,
                std::uint32_t resolution, std::string_view assembly, std::size_t n_threads,
                std::size_t chunk_size, const std::filesystem::path& tmpdir_,
                std::uint32_t compression_lvl, bool skip_all_vs_all_matrix);
  HiCFileWriter(const std::filesystem::path& path_, const hictkpy::BinTable& bins_,
                std::string_view assembly, std::size_t n_threads, std::size_t chunk_size,
                const std::filesystem::path& tmpdir_, std::uint32_t compression_lvl,
                bool skip_all_vs_all_matrix);

  [[nodiscard]] std::filesystem::path path() const noexcept;
  [[nodiscard]] auto resolutions() const;

  [[nodiscard]] const hictk::Reference& chromosomes() const;
  [[nodiscard]] hictkpy::BinTable bins(std::uint32_t resolution) const;

  void add_pixels_from_dict(const nanobind::dict& columns, bool validate);
  void add_pixels_from_table(const nanobind::object& df, bool validate);

  [[nodiscard]] File finalize(std::optional<std::string_view> log_lvl_str = {});
  [[nodiscard]] bool finalized() const noexcept;

  void try_cleanup() noexcept;

  [[nodiscard]] std::string repr() const;

  static void bind(nanobind::module_& m);

 private:
  void add_pixels(const PyArrowTable& table, bool validate);

  const std::filesystem::path& tmpdir() const;
  const hictk::hic::internal::HiCFileWriter& get() const;
  hictk::hic::internal::HiCFileWriter& get();
};

}  // namespace hictkpy
