// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/to_pyarrow.hpp"

#include <Python.h>
#include <arrow/array/array_base.h>
#include <arrow/buffer.h>
#include <arrow/c/abi.h>
#include <arrow/c/bridge.h>
#include <arrow/table.h>
#include <fmt/format.h>
#include <fmt/ranges.h>  // TODO remove

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"

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

[[nodiscard]] static std::shared_ptr<arrow::Schema> import_arrow_schema(const nb::object& df) {
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

[[nodiscard]] static std::shared_ptr<arrow::Table> import_arrow_table(
    const nb::object& df, const std::vector<std::string>& column_names) {
  auto schema = import_arrow_schema(df);

  if (column_names.empty()) {
    return import_arrow_table(df, schema->field_names());
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
  return arrow::Table::Make(std::move(schema), std::move(columns));
}

[[nodiscard]] static std::shared_ptr<arrow::Table> import_pandas_dataframe_slow(
    const nb::object& df, const std::vector<std::string>& column_names) {
  // TODO, how can we deal with columns with unknown type? Can we do that from arrow?
  throw std::runtime_error("not implemented!");
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

[[nodiscard]] static std::shared_ptr<arrow::Table> import_pandas_dataframe(
    const nb::object& df, const std::vector<std::string>& column_names) {
  auto log_failure = [](const char* msg) {
    SPDLOG_WARN(FMT_STRING("failed to import pandas.DataFrame using pyarrow ({}): falling back to "
                           "a slower method"),
                msg);
  };

  try {
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
  } catch (const nb::builtin_exception& e) {
    if (e.type() != nb::exception_type::import_error) {
      log_failure(e.what());
    }
  } catch (const std::exception& e) {
    log_failure(e.what());
  } catch (...) {
    log_failure("unknown error");
  }
  return import_pandas_dataframe_slow(df, column_names);
}

std::shared_ptr<arrow::Table> import_pyarrow_table(const nb::object& df,
                                                   const std::vector<std::string>& column_names) {
  assert(!column_names.empty());
  if (is_pyarrow_table(df)) {
    return import_arrow_table(df, column_names);
  }

  if (is_pandas_dataframe(df)) {
    return import_pandas_dataframe(df, column_names);
  }

  auto type_name = df.type();
  if (nb::hasattr(type_name, "__name__")) {
    type_name = type_name.attr("__name__");
  }

  throw std::invalid_argument(fmt::format(
      FMT_STRING("expected table to be of type pandas.DataFrame or pyarrow.Table, found {}"),
      nb::str(type_name).c_str()));
}

void test_import_table(const nb::object& df, const std::vector<std::string>& column_names) {
  fmt::println(stderr, FMT_STRING("{}"), import_pyarrow_table(df, column_names)->ToString());
}

}  // namespace hictkpy
