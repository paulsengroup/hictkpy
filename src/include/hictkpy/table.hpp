// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <hictk/generic_variant.hpp>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "hictkpy/nanobind.hpp"

namespace arrow {
class Table;

}  // namespace arrow

namespace hictkpy {

[[nodiscard]] nanobind::object export_pyarrow_table(std::shared_ptr<arrow::Table> arrow_table);

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

  PyArrowTable& operator=(const PyArrowTable& other);
  PyArrowTable& operator=(PyArrowTable&& other) noexcept;

  [[nodiscard]] std::shared_ptr<arrow::Table> get() const noexcept;
  [[nodiscard]] auto type() const noexcept -> Type;
  [[nodiscard]] bool has_owner() const noexcept;
  void set_owner(nanobind::object owner);
};

[[nodiscard]] PyArrowTable import_pyarrow_table(const nanobind::object& df,
                                                const std::vector<std::string>& column_names = {});

[[nodiscard]] std::optional<hictk::internal::GenericVariant> infer_column_dtype(
    const std::shared_ptr<const arrow::Table>& df, const std::string& column_name);

}  // namespace hictkpy
