// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filesystem>
#include <hictk/cooler/cooler.hpp>
#include <hictk/file.hpp>
#include <hictk/hic.hpp>
#include <optional>
#include <string>
#include <string_view>

#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"

namespace hictkpy {

class File {
  std::optional<hictk::File> _fp{};
  std::string _uri{};

 public:
  File() = delete;
  explicit File(hictk::File f);
  explicit File(hictk::cooler::File f);
  explicit File(hictk::hic::File f);
  File(const std::filesystem::path& path, std::optional<std::int32_t> resolution,
       std::string_view matrix_type, std::string_view matrix_unit);

  File(const File&) = delete;
  File(File&&) noexcept = default;

  ~File() noexcept;

  File& operator=(const File&) = delete;
  File& operator=(File&&) noexcept = default;

  [[nodiscard]] const hictk::File& operator*() const;
  [[nodiscard]] hictk::File& operator*();

  [[nodiscard]] const hictk::File* operator->() const;
  [[nodiscard]] hictk::File* operator->();

  void close();
  bool try_close() noexcept;

  static void bind(nanobind::module_& m);

  [[nodiscard]] static bool is_cooler(const std::filesystem::path& uri);
  [[nodiscard]] static bool is_hic(const std::filesystem::path& uri);
  [[nodiscard]] CoolerGlobalLock::UniqueLock lock() const;
};

namespace cooler {
void init_global_state();
}

}  // namespace hictkpy
