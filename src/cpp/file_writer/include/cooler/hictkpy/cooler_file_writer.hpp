// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/bin_table.hpp>
#include <hictk/cooler/singlecell_cooler.hpp>
#include <hictk/reference.hpp>
#include <hictk/tmpdir.hpp>
#include <memory>
#include <optional>
#include <string>
#include <string_view>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/file.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/reference.hpp"
#include "hictkpy/table.hpp"

namespace hictkpy {

class CoolerFileWriter {
  std::filesystem::path _path{};
  std::optional<hictk::internal::TmpDir> _tmpdir{};
  std::optional<hictk::cooler::SingleCellFile> _w{};
  std::uint32_t _compression_lvl{};

 public:
  CoolerFileWriter() = delete;
  CoolerFileWriter(std::filesystem::path path_, const ChromosomeDict& chromosomes_,
                   std::uint32_t resolution_, std::string_view assembly,
                   const std::filesystem::path& tmpdir_, std::uint32_t compression_lvl);
  CoolerFileWriter(std::filesystem::path path_, const hictkpy::BinTable& bins_,
                   std::string_view assembly, const std::filesystem::path& tmpdir,
                   std::uint32_t compression_lvl);

  [[nodiscard]] const std::filesystem::path& path() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const;

  [[nodiscard]] const hictk::Reference& chromosomes() const;
  [[nodiscard]] std::shared_ptr<const hictk::BinTable> bins_ptr() const;

  void add_pixels_from_dict(const nanobind::dict& columns, bool sorted, bool validate);
  void add_pixels_from_table(const nanobind::object& df, bool sorted, bool validate);

  [[nodiscard]] File finalize(std::optional<std::string_view> log_lvl_str, std::size_t chunk_size,
                              std::size_t update_frequency);
  [[nodiscard]] bool finalized() const noexcept;

  void try_cleanup() noexcept;

  [[nodiscard]] std::string repr() const;
  static void bind(nanobind::module_& m);

 private:
  [[nodiscard]] static std::optional<hictk::cooler::SingleCellFile> create_file(
      std::string_view path, const hictk::BinTable& bins, std::string_view assembly,
      const std::filesystem::path& tmpdir_);

  void add_pixels(const PyArrowTable& table, bool sorted, bool validate);

  void reset();

  [[nodiscard]] const std::filesystem::path& tmpdir() const;
  [[nodiscard]] const hictk::cooler::SingleCellFile& get() const;
  [[nodiscard]] hictk::cooler::SingleCellFile& get();
};

}  // namespace hictkpy
