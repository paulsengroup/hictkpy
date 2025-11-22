// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/to_pyarrow.hpp"

#include <Python.h>
#include <arrow/c/abi.h>
#include <arrow/c/bridge.h>
#include <arrow/table.h>
#include <fmt/format.h>

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
// NOLINTEND(*-owning-memory,*-no-malloc)

}  // namespace hictkpy
