// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <arrow/type_fwd.h>

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "hictkpy/nanobind.hpp"
#include "hictkpy/variant.hpp"

namespace hictkpy {

[[nodiscard]] nanobind::object export_pyarrow_table(std::shared_ptr<arrow::Table> arrow_table);

// clang-format off
inline constexpr std::array<std::string_view, 3> coo_columns{
  "bin1_id", "bin2_id", "count"
};
inline constexpr std::array<std::string_view, 3> bed3_columns{
  "chrom", "start", "end"
};
inline constexpr std::array<std::string_view, 7> bg2_columns{
  "chrom1", "start1", "end1",
  "chrom2", "start2", "end2",
  "count"
};
// clang-format on

class PyArrowTable {
 public:
  enum class Type : std::uint_fast8_t { BED3, COO, BG2, UNKNOWN };

 private:
  std::optional<nanobind::object> _owner{};
  std::shared_ptr<arrow::Table> _table{};
  Type _table_type{Type::UNKNOWN};

 public:
  PyArrowTable() = default;
  explicit PyArrowTable(std::shared_ptr<arrow::Table> table,
                        std::optional<nanobind::object> owner = {});
  PyArrowTable(std::shared_ptr<arrow::Table> table, Type table_type,
               std::optional<nanobind::object> owner = {});
  PyArrowTable(const PyArrowTable& other);
  PyArrowTable(PyArrowTable&& other) noexcept;
  ~PyArrowTable() noexcept;

  explicit operator bool() const noexcept;
  PyArrowTable& operator=(const PyArrowTable& other);
  PyArrowTable& operator=(PyArrowTable&& other) noexcept;

  [[nodiscard]] std::shared_ptr<arrow::Table> get() const noexcept;
  [[nodiscard]] auto type() const noexcept -> Type;
  [[nodiscard]] bool has_owner() const noexcept;
  void set_owner(nanobind::object owner);
};

[[nodiscard]] PyArrowTable import_pyarrow_table(const nanobind::object& df,
                                                const std::vector<std::string>& column_names = {});

[[nodiscard]] Dtype infer_column_dtype(const std::shared_ptr<const arrow::Table>& df,
                                       const std::string& column_name);

[[nodiscard]] constexpr bool is_dictionary_dtype(arrow::Type::type type) noexcept;
[[nodiscard]] constexpr bool is_string_dtype(arrow::Type::type type) noexcept;
[[nodiscard]] constexpr bool is_integral_dtype(arrow::Type::type type) noexcept;
[[nodiscard]] constexpr bool is_floating_point_dtype(arrow::Type::type type) noexcept;
[[nodiscard]] constexpr bool is_numeric_dtype(arrow::Type::type type) noexcept;

}  // namespace hictkpy

#include "../../table_impl.hpp"
