// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/table.hpp"

#include <arrow/c/abi.h>
#include <arrow/c/bridge.h>
#include <arrow/chunked_array.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <fmt/format.h>

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictkpy/common.hpp"
#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/type.hpp"
#include "hictkpy/variant.hpp"

namespace nb = nanobind;

template <>
struct std::default_delete<ArrowSchema> {
  void operator()(ArrowSchema* schema) const noexcept {
    HICTKPY_GIL_SCOPED_ACQUIRE
    if (schema->release) {
      schema->release(schema);
    }
    delete schema;  // NOLINT(*-owning-memory)
  }
};

template <>
struct std::default_delete<ArrowArrayStream> {
  void operator()(ArrowArrayStream* array) const noexcept {
    HICTKPY_GIL_SCOPED_ACQUIRE
    if (array->release) {
      array->release(array);
    }
    delete array;  // NOLINT(*-owning-memory)
  }
};

namespace hictkpy {

template <typename Container>
[[nodiscard]] static bool str_contains(const Container& buff, std::string_view query) {
  const auto match =
      std::find_if(buff.begin(), buff.end(), [&](const auto& x) { return x == query; });
  return match != buff.end();
}

[[nodiscard]] static bool starts_with(std::string_view str, std::string_view prefix) noexcept {
  return str.find(prefix) == 0;
}

[[nodiscard]] static bool is_valid_chrom_col(std::string_view name,
                                             const arrow::ChunkedArray& col) noexcept {
  if (!starts_with(name, "chrom")) {
    return false;
  }

  const auto tid = col.type()->id();

  if (is_string_dtype(tid)) {
    return true;
  }

  if (!is_dictionary_dtype(tid)) {
    return false;
  }

  const auto dtype = std::static_pointer_cast<arrow::DictionaryType>(col.type());
  return is_string_dtype(dtype->value_type()->id());
}

[[nodiscard]] static bool is_valid_pos_col(std::string_view name,
                                           const arrow::ChunkedArray& col) noexcept {
  if (!starts_with(name, "start") && !starts_with(name, "end")) {
    return false;
  }

  const auto tid = col.type()->id();
  return is_integral_dtype(tid);
}

[[nodiscard]] static bool is_valid_bin_id_col(std::string_view name,
                                              const arrow::ChunkedArray& col) noexcept {
  if (!starts_with(name, "bin")) {
    return false;
  }

  const auto tid = col.type()->id();
  return is_integral_dtype(tid);
}

[[nodiscard]] static bool is_valid_count_col(std::string_view name,
                                             const arrow::ChunkedArray& col) noexcept {
  if (name != "count") {
    return false;
  }

  const auto tid = col.type()->id();
  return is_numeric_dtype(tid);
}

template <std::size_t N>
[[nodiscard]] static PyArrowTable::Type infer_table_type_helper(
    const std::shared_ptr<const arrow::Table>& df,
    const std::array<std::string_view, N>& column_names) {
  using T = PyArrowTable::Type;
  for (const auto& col_name : column_names) {
    if (const auto status = df->schema()->CanReferenceFieldByName(col_name); !status.ok()) {
      throw std::invalid_argument(
          fmt::format(FMT_STRING("unable to uniquely reference column \"{}\" by name: {}"),
                      col_name, status.message()));
    }

    const auto col = df->GetColumnByName(std::string{col_name});
    // clang-format off
    if (!is_valid_chrom_col(col_name, *col)  &&
        !is_valid_bin_id_col(col_name, *col) &&
        !is_valid_pos_col(col_name, *col)    &&
        !is_valid_count_col(col_name, *col)) {
      return T::UNKNOWN;
    }
    // clang-format on
  }

  if (bg2_columns.data() == column_names.data()) {
    return T::BG2;
  }
  if (coo_columns.data() == column_names.data()) {
    return T::COO;
  }
  if (bed3_columns.data() == column_names.data()) {
    return T::BED3;
  }

  unreachable_code();
}

[[nodiscard]] static PyArrowTable::Type infer_table_type(
    const std::shared_ptr<const arrow::Table>& df) {
  if (!df) {
    return PyArrowTable::Type::UNKNOWN;
  }

  std::uint8_t coo_cols_found{};
  std::uint8_t bed3_cols_found{};
  std::uint8_t bg2_cols_found{};

  for (const auto& col : df->schema()->field_names()) {
    if (col == "count") {
      ++bg2_cols_found;
      ++coo_cols_found;
      continue;
    }
    if (str_contains(bg2_columns, col)) {
      ++bg2_cols_found;
      continue;
    }
    if (str_contains(coo_columns, col)) {
      ++coo_cols_found;
      continue;
    }
    if (str_contains(bed3_columns, col)) {
      ++bed3_cols_found;
    }
  }

  if (bg2_cols_found == bg2_columns.size()) {
    return infer_table_type_helper(df, bg2_columns);
  }
  if (coo_cols_found == coo_columns.size()) {
    return infer_table_type_helper(df, coo_columns);
  }
  if (bed3_cols_found == bed3_columns.size()) {
    return infer_table_type_helper(df, bed3_columns);
  }
  return PyArrowTable::Type::UNKNOWN;
}

// https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#pycapsule-standard
// https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html

[[nodiscard]] static nb::object export_arrow_schema(const arrow::Schema& schema_in) {
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#arrowschema-export
  auto schema = std::make_unique<ArrowSchema>();

  const auto status = arrow::ExportSchema(schema_in, schema.get());
  if (!status.ok()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to export arrow::Schema as ArrowSchema: {}"), status.message()));
  }

  HICTKPY_GIL_SCOPED_ACQUIRE
  auto capsule = make_capsule(std::move(schema), "arrow_schema");

  auto obj = nb::module_::import_("types").attr("SimpleNamespace")();
  obj.attr("__setattr__")("__arrow_c_schema__",
                          nb::cpp_function([capsule_ = std::move(capsule)]() { return capsule_; }));
  return obj;
}

[[nodiscard]] static nb::object export_arrow_array_stream(
    std::shared_ptr<arrow::ChunkedArray> column) {
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#arrowstream-export
  auto array_stream = std::make_unique<ArrowArrayStream>();
  const auto status = arrow::ExportChunkedArray(std::move(column), array_stream.get());
  if (!status.ok()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to export arrow::ChunkedArray as ArrowArrayStream: {}"),
                    status.message()));
  }

  HICTKPY_GIL_SCOPED_ACQUIRE
  auto capsule = make_capsule(std::move(array_stream), "arrow_array_stream");
  auto obj = nb::module_::import_("types").attr("SimpleNamespace")();
  obj.attr("__setattr__")(
      "__arrow_c_stream__",
      nb::cpp_function(
          [capsule_ = std::move(capsule)]([[maybe_unused]] const nb::any& _) { return capsule_; },
          nb::arg("requested_schema") = nb::none()));
  return obj;
}

nb::object export_pyarrow_table(std::shared_ptr<arrow::Table> arrow_table) {
  assert(arrow_table);

  HICTKPY_GIL_SCOPED_ACQUIRE
  const auto pa = import_pyarrow_checked();

  std::vector columns(arrow_table->columns().begin(), arrow_table->columns().end());
  std::vector<nb::object> columns_py(columns.size());

  auto schema = pa.attr("schema")(export_arrow_schema(*arrow_table->schema()));
  arrow_table.reset();

  for (std::size_t i = 0; i < columns.size(); ++i) {
    // NOLINTNEXTLINE(*-pro-bounds-pointer-arithmetic)
    columns_py[i] = pa.attr("chunked_array")(export_arrow_array_stream(std::move(columns[i])));
  }

  return nb::cast(pa.attr("Table").attr("from_arrays")(columns_py, nb::arg("schema") = schema),
                  nb::rv_policy::take_ownership);
}

PyArrowTable::PyArrowTable(std::shared_ptr<arrow::Table> table, std::optional<nb::object> owner)
    : _owner(std::move(owner)), _table(std::move(table)), _table_type(infer_table_type(_table)) {}
PyArrowTable::PyArrowTable(std::shared_ptr<arrow::Table> table, Type table_type,
                           std::optional<nb::object> owner)
    : _owner(std::move(owner)), _table(std::move(table)), _table_type(table_type) {}
PyArrowTable::PyArrowTable(const PyArrowTable& other)
    : _table(other._table), _table_type(other._table_type) {
  if (other.has_owner()) {
    HICTKPY_GIL_SCOPED_ACQUIRE
    _owner = other._owner;
  }
}
PyArrowTable::PyArrowTable(PyArrowTable&& other) noexcept
    : _table(std::move(other._table)), _table_type(other._table_type) {
  if (other.has_owner()) {
    HICTKPY_GIL_SCOPED_ACQUIRE
    _owner = std::move(other._owner);
  }
}

PyArrowTable::~PyArrowTable() noexcept {
  if (has_owner()) {
    HICTKPY_GIL_SCOPED_ACQUIRE
    _owner.reset();
  }
}

PyArrowTable::operator bool() const noexcept { return !!_table; }

PyArrowTable& PyArrowTable::operator=(const PyArrowTable& other) {
  if (this == &other) {
    return *this;
  }

  {
    HICTKPY_GIL_SCOPED_ACQUIRE
    _owner = other._owner;
  }
  _table = other._table;
  _table_type = other._table_type;

  return *this;
}
PyArrowTable& PyArrowTable::operator=(PyArrowTable&& other) noexcept {
  if (this == &other) {
    return *this;
  }

  {
    HICTKPY_GIL_SCOPED_ACQUIRE
    _owner = std::move(other._owner);
  }
  _table = std::move(other._table);
  _table_type = other._table_type;

  return *this;
}

std::shared_ptr<arrow::Table> PyArrowTable::get() const noexcept {
  assert(_table);
  return _table;
}
auto PyArrowTable::type() const noexcept -> Type { return _table_type; }

bool PyArrowTable::has_owner() const noexcept { return _owner.has_value(); }
void PyArrowTable::set_owner(nb::object owner) { _owner = std::move(owner); }

[[nodiscard]] static std::shared_ptr<arrow::Schema> import_arrow_schema(const nb::object& df) {
  HICTKPY_GIL_SCOPED_ACQUIRE
  if (!nb::hasattr(df, "schema")) {
    throw nb::attribute_error("object does not have attribute \"schema\"");
  }

  if (!nb::hasattr(df.attr("schema"), "__arrow_c_schema__")) {
    throw nb::attribute_error("schema object does not have attribute \"__arrow_c_schema__\"");
  }

  auto py_schema = nb::cast<nb::capsule>(df.attr("schema").attr("__arrow_c_schema__")());
  if (py_schema.name() != std::string_view{"arrow_schema"}) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("expected capsule of with name arrow_schema, found {}"), py_schema.name()));
  }

  auto schema = arrow::ImportSchema(static_cast<ArrowSchema*>(py_schema.data()));
  if (!schema.ok()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to import arrow::Schema from pyarrow.Table.schema: {}"),
                    schema.status().message()));
  }
  return schema.MoveValueUnsafe();
}

[[nodiscard]] static std::shared_ptr<arrow::ChunkedArray> import_arrow_column(
    const nb::object& df, std::string_view column_name) {
  HICTKPY_GIL_SCOPED_ACQUIRE
  if (!nb::hasattr(df, "column")) {
    throw nb::attribute_error("object does not have attribute \"column\"");
  }

  auto py_column = df.attr("column")(column_name);

  if (!nb::hasattr(py_column, "__arrow_c_stream__")) {
    const auto msg = fmt::format(
        FMT_STRING("column {} does not have attribute \"__arrow_c_schema__\""), column_name);
    throw nb::attribute_error(msg.c_str());
  }

  auto py_array = nb::cast<nb::capsule>(py_column.attr("__arrow_c_stream__")());
  if (py_array.name() != std::string_view{"arrow_array_stream"}) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("expected capsule of with name arrow_schema, found {}"), py_array.name()));
  }

  auto array = arrow::ImportChunkedArray(static_cast<ArrowArrayStream*>(py_array.data()));
  if (!array.ok()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to import column \"{}\" from pyarrow.Table: {}"),
                    column_name, array.status().message()));
  }
  return array.MoveValueUnsafe();
}

[[nodiscard]] static PyArrowTable import_arrow_table(const nb::object& df,
                                                     const std::vector<std::string>& column_names) {
  auto schema = import_arrow_schema(df);

  if (column_names.empty()) {
    const auto field_names = schema->field_names();
    if (field_names.empty()) {
      return PyArrowTable{nullptr, PyArrowTable::Type::UNKNOWN};
    }
    return import_arrow_table(df, field_names);
  }

  arrow::FieldVector fields{};
  fields.reserve(column_names.size());

  for (std::string_view column_name : column_names) {
    auto field = schema->GetFieldByName(column_name);
    if (!field) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("table does not have a column named \"{}\""), column_name));
    }
    fields.emplace_back(std::move(field));
  }

  arrow::ChunkedArrayVector columns{};
  columns.reserve(column_names.size());

  for (const auto& column_name : column_names) {
    columns.emplace_back(import_arrow_column(df, column_name));
  }

  schema = std::make_shared<arrow::Schema>(std::move(fields), schema->metadata());

  auto table = arrow::Table::Make(std::move(schema), std::move(columns));
  HICTKPY_GIL_SCOPED_ACQUIRE
  return PyArrowTable{std::move(table), df};
}

[[nodiscard]] static bool is_pyarrow_table(const nb::object& df) noexcept {
  try {
    HICTKPY_GIL_SCOPED_ACQUIRE
    auto pa = import_pyarrow_checked();
    return nb::isinstance(df, pa.attr("Table"));
  } catch (...) {
    return false;
  }
}

[[nodiscard]] static bool is_pandas_dataframe(const nb::object& df) noexcept {
  try {
    HICTKPY_GIL_SCOPED_ACQUIRE
    auto pd = import_module_checked("pandas");
    return nb::isinstance(df, pd.attr("DataFrame"));
  } catch (...) {
    return false;
  }
}

[[nodiscard]] static PyArrowTable import_pandas_dataframe(
    const nb::object& df, const std::vector<std::string>& column_names) {
  auto pyarrow_table = [&]() {
    HICTKPY_GIL_SCOPED_ACQUIRE
    auto pa = import_pyarrow_checked();
    if (column_names.empty()) {
      return std::make_optional(pa.attr("Table").attr("from_pandas")(df));
    }
    return std::make_optional(
        pa.attr("Table").attr("from_pandas")(df, nb::arg("columns") = column_names));
  }();

  // NOLINTNEXTLINE(*-unchecked-optional-access)
  auto table = import_arrow_table(*pyarrow_table, column_names);
  {
    HICTKPY_GIL_SCOPED_ACQUIRE
    pyarrow_table.reset();
  }
  return table;
}

PyArrowTable import_pyarrow_table(const nb::object& df,
                                  const std::vector<std::string>& column_names) {
  try {
    check_pyarrow_is_importable();
  } catch (nb::python_error& e) {  // NOLINTNEXTLINE(*-vararg)
    nb::raise_from(e, PyExc_ModuleNotFoundError,
                   "Loading interactions from a DataFrame requires pyarrow");
  }

  if (is_pyarrow_table(df)) {
    return import_arrow_table(df, column_names);
  }

  if (is_pandas_dataframe(df)) {
    return import_pandas_dataframe(df, column_names);
  }

  HICTKPY_GIL_SCOPED_ACQUIRE
  throw std::invalid_argument(fmt::format(
      FMT_STRING("expected table to be of type pandas.DataFrame or pyarrow.Table, found {}"),
      format_py_type(df.type())));
}

Dtype infer_column_dtype(const std::shared_ptr<const arrow::Table>& df,
                         const std::string& column_name) {
  const auto col = df->GetColumnByName(column_name);
  if (!col) {
    return {};
  }

  const auto dtype = col->type();

  using T = arrow::Type::type;
  switch (dtype->id()) {
    case T::UINT8:
      return std::uint8_t{};
    case T::UINT16:
      return std::uint16_t{};
    case T::UINT32:
      return std::uint32_t{};
    case T::UINT64:
      return std::uint64_t{};
    case T::INT8:
      return std::int8_t{};
    case T::INT16:
      return std::int16_t{};
    case T::INT32:
      return std::int32_t{};
    case T::INT64:
      return std::int64_t{};
    case T::FLOAT:
      return float{};
    case T::DOUBLE:
      return double{};
    case T::STRING:
      [[fallthrough]];
    case T::STRING_VIEW:
      [[fallthrough]];
    case T::LARGE_STRING:
      return std::string{};
    default: {
      assert(dtype->id() != T::DICTIONARY);
      return {};
    }
  }
}

}  // namespace hictkpy
